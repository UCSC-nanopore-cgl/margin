﻿#!/usr/bin/env python2.7
from __future__ import print_function

import argparse
import os
import multiprocessing
import subprocess
import sys
import textwrap
import tarfile
from urlparse import urlparse
import math
import shutil
import glob
import time
import datetime
import logging
from contextlib import closing

import yaml
from bd2k.util.files import mkdir_p
from bd2k.util.processes import which
from toil.job import Job
from toil.lib.docker import dockerCall, dockerCheckOutput, _fixPermissions
from toil_lib import require, UserError
from toil_lib.files import tarball_files, copy_files
from toil_lib.jobs import map_job
from toil_lib.urls import download_url, s3am_upload

# formats
SCHEMES = ('http', 'file', 's3', 'ftp')

# filenames
DEFAULT_CONFIG_NAME = 'config-toil-marginphase.yaml'
DEFAULT_MANIFEST_NAME = 'manifest-toil-marginphase.tsv'

# docker images
DOCKER_SAMTOOLS = "quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c"
DOCKER_MARGIN_PHASE = "quay.io/ucsc_cgl/margin_phase:latest"

# resource
MP_CPU = 2
MP_MEM_BAM_FACTOR = 1024 #todo account for learning iterations
MP_MEM_REF_FACTOR = 2
MP_DSK_BAM_FACTOR = 5 #input bam chunk, output (in sam fmt), vcf etc
MP_DSK_REF_FACTOR = 2

# for debugging (shouldn't be committed as true)
DEBUG = False
DOCKER_LOGGING = False

#chunk_info_keys
CI_CHUNK_BOUNDARY_START = "chunk_bounary_start_pos" #position where chunk split occured
CI_CHUNK_BOUNDARY_END = "chunk_bounary_end_pos" #position where chunk split occured
CI_CHUNK_START = "chunk_start_pos" #chunk boundary modified by the margin
CI_CHUNK_END = "chunk_end_pos" #chunk boundary modified by the margin
CI_READ_COUNT = "read_count" #how many reads were in the chunk
CI_CHUNK_SIZE = "chunk_size" #chunk size in bytes
CI_CHUNK_INDEX = "chunk_index" #index of the chunk
CI_OUTPUT_FILE_ID = "output_file_id"

# output naming conventions
SAM_HAP_1_SUFFIX = "out.1.sam"
SAM_HAP_2_SUFFIX = "out.2.sam"
VCF_SUFFIX = "out.vcf"

def parse_samples(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """

    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if line.isspace() or line.startswith('#'):
                continue
            sample = line.strip().split('\t')

            # validate structure
            if len(sample) != 4:
                raise UserError('Bad manifest format! Expected 4 tab-separated columns, got: {}'.format(sample))

            # extract sample parts
            uuid, url, contig_name, reference_url = sample

            # validation?

            sample = [uuid, url, contig_name, reference_url]
            samples.append(sample)
    return samples


def prepare_input(job, sample, config):

    # job prep
    config = argparse.Namespace(**vars(config))
    uuid, url, contig_name, reference_url = sample
    config.uuid = uuid
    config.contig_name = contig_name
    config.reference_url = reference_url
    work_dir = job.fileStore.getLocalTempDir()
    start = time.time()
    job.fileStore.logToMaster("START:{}:{}".format(config.uuid, datetime.datetime.now()))

    # todo global resource estimation
    config.maxCores = min(config.maxCores, multiprocessing.cpu_count())
    if DEBUG:
        config.maxMemory = min(config.maxMemory, 210 * 1024 * 1024 * 1024)
    #config.disk

    # download references
    #ref fasta
    download_url(job, url=reference_url, work_dir=work_dir)
    ref_genome_filename = os.path.basename(reference_url)
    ref_genome_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, ref_genome_filename))
    config.reference_genome_fileid = ref_genome_fileid
    ref_genome_size = os.stat(os.path.join(work_dir, ref_genome_filename)).st_size
    #ref vcf
    download_url(job, url=config.reference_vcf, work_dir=work_dir)
    ref_vcf_filename = os.path.basename(config.reference_vcf)
    ref_vcf_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, ref_vcf_filename))
    config.reference_vcf_fileid = ref_vcf_fileid
    #params
    download_url(job, url=config.params, work_dir=work_dir)
    params_filename = os.path.basename(config.params)
    params_fileid = job.fileStore.writeGlobalFile(os.path.join(work_dir, params_filename))
    config.params_fileid = params_fileid


    # download bam
    download_url(job, url=url, work_dir=work_dir)
    bam_filename = os.path.basename(url)
    data_bam_location = os.path.join("/data", bam_filename)

    # index the bam
    docker_params = ["index", data_bam_location]
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, docker_params))
    dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=docker_params)

    # sanity check
    workdir_bai_location = os.path.join(work_dir, bam_filename + ".bai")
    if not os.path.isfile(workdir_bai_location):
        raise UserError("BAM index file not created for {}: {}".format(bam_filename, workdir_bai_location))

    # get start and end location
    get_idx_cmd = [
        ["samtools", "view", data_bam_location],
        ["head", "-n", "1"],
        [os.path.join("/data", _write_select_column_script(work_dir))]
    ]
    start_idx_str = _dockerCheckOutput_except_141(job, tool=DOCKER_SAMTOOLS, work_dir=work_dir, parameters=get_idx_cmd).strip()
    get_idx_cmd[1][0] = "tail"
    end_idx_str = _dockerCheckOutput_except_141(job, tool=DOCKER_SAMTOOLS, work_dir=work_dir, parameters=get_idx_cmd).strip()
    job.fileStore.logToMaster("{}: start_pos:{}, end_pos:{}".format(config.uuid, start_idx_str, end_idx_str))
    start_idx = int(start_idx_str) - 1
    end_idx = int(end_idx_str) + 1

    # get reads from positions
    chunk_infos = list()
    idx = start_idx
    while idx < end_idx:
        ci = dict()
        ci[CI_CHUNK_BOUNDARY_START] = idx
        chunk_start = idx - (config.partition_margin if idx != start_idx else 0)
        ci[CI_CHUNK_START] = chunk_start
        idx += config.partition_size
        ci[CI_CHUNK_BOUNDARY_END] = idx
        chunk_end = idx + config.partition_margin
        ci[CI_CHUNK_END] = chunk_end
        chunk_infos.append(ci)

    # enqueue jobs
    job.fileStore.logToMaster("{}: Enqueueing {} jobs".format(config.uuid, len(chunk_infos)))
    idx = 0
    enqueued_jobs = 0
    returned_tarballs = list()
    total_mem = 0
    for ci in chunk_infos:
        #prep
        ci[CI_CHUNK_INDEX] = idx
        chunk_start = ci[CI_CHUNK_START]
        chunk_end = ci[CI_CHUNK_END]
        chunk_position_description = "{}:{}-{}".format(config.contig_name, chunk_start, chunk_end)
        bam_split_command = ["view", "-b", data_bam_location, chunk_position_description]
        chunk_name = "{}.{}.bam".format(config.uuid, idx)
        #write chunk
        chunk_location = os.path.join(work_dir, chunk_name)
        with open(chunk_location, 'w') as out:
            if DOCKER_LOGGING:
                job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, bam_split_command))
            dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=bam_split_command, outfile=out)
        #document read count
        chunk_size = os.stat(chunk_location).st_size
        ci[CI_CHUNK_SIZE] = chunk_size
        read_count= _get_bam_read_count(job, work_dir, chunk_name)
        ci[CI_READ_COUNT] = read_count
        job.fileStore.logToMaster("{}: chunk from {} for idx {} is {}b ({}mb) and has {} reads"
                                  .format(config.uuid, chunk_position_description, idx, chunk_size,
                                          int(chunk_size / 1024 / 1024), read_count))
        if DEBUG:
            job.fileStore.logToMaster("Copying {} to {}".format(chunk_name, config.output_dir))
            copy_files(file_paths=[chunk_location], output_dir=config.output_dir)

        # enqueue marginPhase job
        if read_count > 0:
            chunk_fileid = job.fileStore.writeGlobalFile(chunk_location)
            mp_cores = int(min(MP_CPU, config.maxCores))
            mp_mem = int(min(int(chunk_size * MP_MEM_BAM_FACTOR + ref_genome_size * MP_MEM_REF_FACTOR), config.maxMemory))
            mp_disk = int(min(int(chunk_size * MP_DSK_BAM_FACTOR + ref_genome_size * MP_DSK_REF_FACTOR), config.maxDisk))
            job.fileStore.logToMaster("{}:{} requesting {} cores, {}b ({}mb) disk, {}b ({}gb) mem"
                                      .format(config.uuid, idx, mp_cores, mp_disk, int(mp_disk / 1024 / 1024 ),
                                              mp_mem, int(mp_mem / 1024 / 1024 / 1024)))
            total_mem += mp_mem
            mp_mem = str(int(mp_mem / 1024)) + "K"
            mp_disk = str(int(mp_disk) / 1024) + "K"
            margin_phase_job = job.addChildJobFn(run_margin_phase, config, chunk_fileid, ci,
                                                 memory=mp_mem, cores=mp_cores, disk=mp_disk)
            returned_tarballs.append(margin_phase_job.rv())
            enqueued_jobs += 1
        idx += 1

    job.fileStore.logToMaster("{}: Enqueued {} jobs, requested total of {}gb ({}b) mem"
                              .format(config.uuid, enqueued_jobs, int(total_mem/1024/1024/1024), total_mem))

    # enqueue merging and consolidation job
    merge_job = job.addFollowOnJobFn(merge_chunks, config, returned_tarballs)
    merge_job.addFollowOnJobFn(consolidate_output, config, merge_job.rv())

    # log
    _log_time(job, "prepare_input", start, config.uuid)


def _write_select_column_script(work_dir, column=4):
    # I feel bad for doing this, but I can't send single quotes into a toil command without them becoming escaped
    # so I'm just creating a script which does that
    filename = "select_column_{}.sh".format(column)
    file_location = os.path.join(work_dir, filename)
    with open(file_location, 'w') as out:
        print("#!/usr/bin/awk -f", file=out)
        print("{print $%d}" % column, file=out)
    os.chmod(file_location, 1023) #rwxrwxrwx
    return filename


def _dockerCheckOutput_except_141(job, tool, work_dir, parameters):
    # there's something strange with the return code for commands which stop reading from stdin (like "head")
    # and so we need to ignore the returncode
    try:
        return dockerCheckOutput(job, tool, parameters=parameters, workDir=work_dir)
    except subprocess.CalledProcessError, e:
        if e.returncode == 141:
            return e.output
        else:
            raise e

def _get_bam_read_count(job, work_dir, bam_name):
    params = [
        ["samtools", "view", os.path.join("/data", bam_name)],
        ["wc", "-l"]
    ]
    line_count_str = dockerCheckOutput(job, DOCKER_SAMTOOLS, params, work_dir)
    return int(line_count_str)


def run_margin_phase(job, config, chunk_file_id, chunk_info):
    # prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    chunk_idx = chunk_info[CI_CHUNK_INDEX]
    chunk_identifier = "{}.{}".format(config.uuid, chunk_idx)
    chunk_name = "{}.in.bam".format(chunk_identifier)
    chunk_location = os.path.join(work_dir, chunk_name)
    job.fileStore.logToMaster("run_margin_phase:{}:{}:{}".format(config.uuid, chunk_idx, datetime.datetime.now()))

    # download bam chunk
    job.fileStore.readGlobalFile(chunk_file_id, chunk_location)
    if not os.path.isfile(chunk_location):
        raise UserError("Failed to download chunk {} from {}".format(chunk_name, chunk_file_id))

    # download references
    #ref genome
    genome_reference_name = "reference.fa"
    genome_reference_location = os.path.join(work_dir, genome_reference_name)
    job.fileStore.readGlobalFile(config.reference_genome_fileid, genome_reference_location)
    if not os.path.isfile(genome_reference_location):
        raise UserError("Failed to download genome reference {} from {}"
                        .format(os.path.basename(config.reference_genome), config.reference_genome_fileid))
    #ref vcf
    vcf_reference_name = "reference.vcf"
    vcf_reference_location = os.path.join(work_dir, vcf_reference_name)
    job.fileStore.readGlobalFile(config.reference_vcf_fileid, vcf_reference_location)
    if not os.path.isfile(genome_reference_location):
        raise UserError("Failed to download vcf reference {} from {}"
                        .format(os.path.basename(config.reference_vcf), config.reference_vcf_fileid))
    # params
    params_name = "params.json"
    params_location = os.path.join(work_dir, params_name)
    job.fileStore.readGlobalFile(config.params_fileid, params_location)
    if not os.path.isfile(params_location):
        raise UserError("Failed to download params {} from {}"
                        .format(os.path.basename(config.params), config.params_fileid))

    # run marginPhase
    params = [os.path.join("/data", chunk_name), os.path.join("/data", genome_reference_name),
              "-p", os.path.join("/data", params_name), "-o", os.path.join("/data","{}.out".format(chunk_identifier)),
              "-r",  os.path.join("/data", vcf_reference_name)]
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_MARGIN_PHASE, params))
    dockerCall(job, tool=DOCKER_MARGIN_PHASE, workDir=work_dir, parameters=params)
    os.rename(os.path.join(work_dir, "marginPhase.log"), os.path.join(work_dir, "{}.log".format(chunk_identifier)))

    # document output
    job.fileStore.logToMaster("Output files for {}:".format(chunk_identifier))
    output_file_locations = glob.glob(os.path.join(work_dir, "{}*".format(chunk_identifier)))
    for f in output_file_locations:
        job.fileStore.logToMaster("\t{}".format(os.path.basename(f)))

    # tarball the output and save
    tarball_name = "{}.tar.gz".format(chunk_identifier)
    tarball_files(tar_name=tarball_name, file_paths=output_file_locations, output_dir=work_dir)
    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
    chunk_info[CI_OUTPUT_FILE_ID] = output_file_id

    # log
    _log_time(job, "run_margin_phase", start, chunk_identifier)
    return chunk_info


def merge_chunks(job, config, chunk_infos):
    # prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    job.fileStore.logToMaster("merging_chunks:{}:{}".format(config.uuid, datetime.datetime.now()))
    job.fileStore.logToMaster("merging {} chunks".format(len(chunk_infos)))
    # work directory for tar management
    tar_work_dir = os.path.join(work_dir, "tmp")
    # output files
    merged_chunks_directory = os.path.join(work_dir, "merged")
    os.mkdir(merged_chunks_directory)
    merged_chunk_idx = 0  # for areas where merging can't happen
    merged_hap1_name, merged_hap2_name, merged_vcf_name = None, None, None
    merged_hap1_file, merged_hap2_file, merged_vcf_file = None, None, None

    # sort by chunk index and validate
    chunk_infos.sort(key=(lambda x: x[CI_CHUNK_INDEX]))
    idx = 0
    missing_indices = []
    for ci in chunk_infos:
        while ci[CI_CHUNK_INDEX] > idx:
            missing_indices.append(idx)
            idx += 1
        idx += 1
    if len(missing_indices) > 0:
        job.fileStore.logToMaster("{}: Found {} missing indices: {}"
                                  .format(config.uuid, len(missing_indices), missing_indices))

    # prep for iteration
    prev_hap1_read_ids, prev_hap2_read_ids = set(), set()
    prev_chunk = {CI_CHUNK_INDEX: "start"}

    # iterate over all chunks
    for chunk in chunk_infos:
        # get current chunk
        if os.path.isdir(tar_work_dir):
            shutil.rmtree(tar_work_dir)
        sam_hap1_file, sam_hap2_file, vcf_file = _extract_chunk_tarball(job, config, tar_work_dir, chunk)

        # get reads
        read_start_pos = chunk[CI_CHUNK_START]
        read_end_pos = chunk[CI_CHUNK_BOUNDARY_START] + config.partition_margin
        curr_hap1_read_ids = _get_read_ids_in_range(job, tar_work_dir, os.path.basename(sam_hap1_file),
                                           config.contig_name, read_start_pos, read_end_pos)
        curr_hap2_read_ids = _get_read_ids_in_range(job, tar_work_dir, os.path.basename(sam_hap2_file),
                                           config.contig_name, read_start_pos, read_end_pos)
        job.fileStore.logToMaster("{}: found {} reads for the start of chunk {} with read boundaries {} - {}"
                                  .format(config.uuid, (len(curr_hap1_read_ids) + len(curr_hap2_read_ids)),
                                          chunk[CI_CHUNK_INDEX], read_start_pos, read_end_pos))

        # log the matching info
        same_haplotype_ordering = _should_same_haplotype_ordering_be_maintained(job, config, prev_chunk, chunk,
                                                                                prev_hap1_read_ids, prev_hap2_read_ids,
                                                                                curr_hap1_read_ids, curr_hap2_read_ids)
        # exclude reads already in a haplotype
        read_ids_to_exclude = set()
        for id in prev_hap1_read_ids:
            read_ids_to_exclude.add(id)
        for id in prev_hap2_read_ids:
            read_ids_to_exclude.add(id)

        # this indicates there was no (or equal) read overlap.  probably it means we've just started the process
        if same_haplotype_ordering is None:
            job.fileStore.logToMaster("{}: starting new merged chunk idx {} from chunk {}"
                                      .format(config.uuid, merged_chunk_idx, chunk[CI_CHUNK_INDEX]))

            # get merged haplotype names and files
            merged_hap1_name = "{}.merged.{}.hap1.sam".format(config.uuid, merged_chunk_idx)
            merged_hap2_name = "{}.merged.{}.hap2.sam".format(config.uuid, merged_chunk_idx)
            merged_vcf_name  = "{}.merged.{}.vcf".format(config.uuid, merged_chunk_idx)
            merged_hap1_file = os.path.join(merged_chunks_directory, merged_hap1_name)
            merged_hap2_file = os.path.join(merged_chunks_directory, merged_hap2_name)
            merged_vcf_file  = os.path.join(merged_chunks_directory, merged_vcf_name)

            # prep the hap1 and hap2 bams - the headers should be the same for all chunks, so we can just start with
            #   the extracted haplotype sam files
            subprocess.check_call(["cp", sam_hap1_file, merged_hap1_file])
            subprocess.check_call(["cp", sam_hap2_file, merged_hap2_file])
            subprocess.check_call(["cp", vcf_file, merged_vcf_file])

            # increment merged chunk idx
            merged_chunk_idx += 1
        elif same_haplotype_ordering:
            job.fileStore.logToMaster("{}:chunk{}: writing same ordering"
                                      .format(config.uuid, chunk[CI_CHUNK_INDEX]))
            #append reads
            excl_ids_hap1 = _append_sam_reads_to_file(job, config, sam_hap1_file, merged_hap1_file, read_ids_to_exclude)
            excl_ids_hap2 = _append_sam_reads_to_file(job, config, sam_hap2_file, merged_hap2_file, read_ids_to_exclude)
            #document excluded reads
            excl_ids_hap1_cnt = len(excl_ids_hap1)
            excl_ids_hap2_cnt = len(excl_ids_hap2)
            job.fileStore.logToMaster("{}:chunk{}:hap1: excluded {} reads ({}% of overlap) during merge"
                                      .format(config.uuid, chunk[CI_CHUNK_INDEX], excl_ids_hap1_cnt,
                                              int(100.0 * excl_ids_hap1_cnt / len(curr_hap1_read_ids))))
            job.fileStore.logToMaster("{}:chunk{}:hap2: excluded {} reads ({}% of overlap) during merge"
                                      .format(config.uuid, chunk[CI_CHUNK_INDEX], excl_ids_hap2_cnt,
                                              int(100.0 * excl_ids_hap2_cnt / len(curr_hap2_read_ids))))
            # append vcf calls
            _append_vcf_calls_to_file(job, config, vcf_file, merged_vcf_file, False)
        else:
            job.fileStore.logToMaster("{}:chunk{}: writing different ordering"
                                      .format(config.uuid, chunk[CI_CHUNK_INDEX]))
            #append reads
            excl_ids_hap1 = _append_sam_reads_to_file(job, config, sam_hap1_file, merged_hap2_file, read_ids_to_exclude)
            excl_ids_hap2 = _append_sam_reads_to_file(job, config, sam_hap2_file, merged_hap1_file, read_ids_to_exclude)
            #document excluded reads
            excl_ids_hap1_cnt = len(excl_ids_hap1)
            excl_ids_hap2_cnt = len(excl_ids_hap2)
            job.fileStore.logToMaster("{}:chunk{}:hap1: excluded {} reads ({}% of overlap) during merge"
                                      .format(config.uuid, chunk[CI_CHUNK_INDEX], excl_ids_hap1_cnt,
                                              int(100.0 * excl_ids_hap1_cnt / len(curr_hap1_read_ids))))
            job.fileStore.logToMaster("{}:chunk{}:hap2: excluded {} reads ({}% of overlap) during merge"
                                      .format(config.uuid, chunk[CI_CHUNK_INDEX], excl_ids_hap2_cnt,
                                              int(100.0 * excl_ids_hap2_cnt / len(curr_hap2_read_ids))))
            # append vcf calls
            _append_vcf_calls_to_file(job, config, vcf_file, merged_vcf_file, True)

        # prep for iteration / cleanup
        read_start_pos = chunk[CI_CHUNK_BOUNDARY_END] - config.partition_margin
        read_end_pos = chunk[CI_CHUNK_END]
        prev_hap1_read_ids = _get_read_ids_in_range(job, tar_work_dir, os.path.basename(sam_hap1_file),
                                                config.contig_name, read_start_pos, read_end_pos)
        prev_hap2_read_ids = _get_read_ids_in_range(job, tar_work_dir, os.path.basename(sam_hap2_file),
                                                config.contig_name, read_start_pos, read_end_pos)
        prev_chunk = chunk
        job.fileStore.logToMaster("{}: found {} reads for the end of chunk {} with read boundaries {} - {}"
                                  .format(config.uuid, (len(curr_hap1_read_ids) + len(curr_hap2_read_ids)),
                                          prev_chunk[CI_CHUNK_INDEX], read_start_pos, read_end_pos))

    # post-processing
    for file_name in os.listdir(merged_chunks_directory):
        if file_name.endswith(".sam"): _sort_sam_file(job, config, merged_chunks_directory, file_name)
        if file_name.endswith(".vcf"): _sort_vcf_file(job, config, merged_chunks_directory, file_name)

    # tarball the output and save
    job.fileStore.logToMaster("{}: Output files for merge:".format(config.uuid))
    output_file_locations = glob.glob(os.path.join(merged_chunks_directory, "*"))
    output_file_locations.sort()
    for f in output_file_locations:
        job.fileStore.logToMaster("\t{}".format(os.path.basename(f)))
    tarball_name = "{}.merged.tar.gz".format(config.uuid)
    tarball_files(tar_name=tarball_name, file_paths=output_file_locations, output_dir=work_dir)
    output_file_id = job.fileStore.writeGlobalFile(os.path.join(work_dir, tarball_name))
    # we need to return the input list of chunk infos for consolidation
    chunk_infos.append({CI_OUTPUT_FILE_ID: output_file_id, CI_CHUNK_INDEX: "merged"})

    _log_time(job, "merge_chunks", start, config.uuid)
    return chunk_infos


def _should_same_haplotype_ordering_be_maintained(job, config, prev_chunk, curr_chunk,
                                                  prev_hap1_read_ids, prev_hap2_read_ids,
                                                  curr_hap1_read_ids, curr_hap2_read_ids):
    # prep
    match_identifier = "{}:read_matching {}-{}".format(config.uuid, prev_chunk[CI_CHUNK_INDEX], curr_chunk[CI_CHUNK_INDEX])

    # total counts
    prev_hap1_read_count = len(prev_hap1_read_ids)
    prev_hap2_read_count = len(prev_hap2_read_ids)
    curr_hap1_read_count = len(curr_hap1_read_ids)
    curr_hap2_read_count = len(curr_hap2_read_ids)
    job.fileStore.logToMaster("{}: \tprev_hap1_cnt:{} \tprev_hap2_cnt:{} \tcurr_hap1_cnt:{} \tcurr_hap2_cnt:{}"
                              .format(match_identifier, prev_hap1_read_count, prev_hap2_read_count,
                                      curr_hap1_read_count, curr_hap2_read_count))
    if len(curr_hap1_read_ids.intersection(curr_hap2_read_ids)) > 0:
        job.fileStore.logToMaster("{}: chunk {} had {} reads in both haplotypes!"
                                  .format(config.uuid, curr_chunk[CI_CHUNK_INDEX],
                                          len(curr_hap1_read_ids.intersection(curr_hap2_read_ids))))

    # intersect counts
    reads_in_curr1_and_prev1 = 0
    reads_in_curr1_and_prev2 = 0
    reads_in_curr1_and_neither_prev = 0
    reads_in_curr2_and_prev1 = 0
    reads_in_curr2_and_prev2 = 0
    reads_in_curr2_and_neither_prev = 0

    # calculate intersection
    for id in curr_hap1_read_ids:
        if id in prev_hap1_read_ids: reads_in_curr1_and_prev1 += 1
        if id in prev_hap2_read_ids: reads_in_curr1_and_prev2 += 1
        if id not in prev_hap1_read_ids and id not in prev_hap2_read_ids: reads_in_curr1_and_neither_prev += 1
    for id in curr_hap2_read_ids:
        if id in prev_hap1_read_ids: reads_in_curr2_and_prev1 += 1
        if id in prev_hap2_read_ids: reads_in_curr2_and_prev2 += 1
        if id not in prev_hap1_read_ids and id not in prev_hap2_read_ids: reads_in_curr2_and_neither_prev += 1

    reads_supporting_same_ordering = reads_in_curr1_and_prev1 + reads_in_curr2_and_prev2
    reads_supporting_different_ordering = reads_in_curr1_and_prev2 + reads_in_curr2_and_prev1
    reads_in_currs_not_in_prevs = curr_hap1_read_count + curr_hap2_read_count \
                                  - reads_in_curr1_and_neither_prev - reads_in_curr2_and_neither_prev

    ratio_supporting_same_ordering, ratio_supporting_different_ordering = -1, -1
    if reads_in_currs_not_in_prevs != 0:
        ratio_supporting_same_ordering = 1.0 * reads_supporting_same_ordering / reads_in_currs_not_in_prevs
        ratio_supporting_different_ordering = 1.0 * reads_supporting_different_ordering / reads_in_currs_not_in_prevs

    # log stuff (maybe this can be removed later)
    job.fileStore.logToMaster("{}: \tcur1_prev1:{} \tcur1_prev2:{} \tcur1_only:{} \tcur2_prev1:{} \tcur2_prev2:{} \tcur2_only:{}"
                              .format(match_identifier, reads_in_curr1_and_prev1, reads_in_curr1_and_prev2,
                                      reads_in_curr1_and_neither_prev, reads_in_curr2_and_prev1,
                                      reads_in_curr2_and_prev2, reads_in_curr2_and_neither_prev))
    job.fileStore.logToMaster("{}: \treads_supporting_current_order:{} ({}) \treads_supporting_different_order:{} ({})"
                              .format(match_identifier, reads_supporting_same_ordering, ratio_supporting_same_ordering,
                                      reads_supporting_different_ordering, ratio_supporting_different_ordering))

    #return recommendation
    # None if no recommendation, else returns whether data indicates same ordering (T or F)
    if ratio_supporting_same_ordering == ratio_supporting_different_ordering: return None
    return ratio_supporting_same_ordering > ratio_supporting_different_ordering


def _sort_sam_file(job, config, work_dir, sam_file_name):
    # prep
    job.fileStore.logToMaster("{}: sorting {}".format(config.uuid, sam_file_name))
    sorted_file_name = "{}.sorted.sam".format(sam_file_name)
    # sort
    sort_cmd = ["sort", "-o", os.path.join("/data/", sorted_file_name),
                os.path.join("/data", sam_file_name)]
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, sort_cmd))
    dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=sort_cmd)
    # replace
    subprocess.check_call(["mv", os.path.join(work_dir, sorted_file_name), os.path.join(work_dir, sam_file_name)])


def _sort_vcf_file(job, config, work_dir, vcf_file_name):
    #prep
    job.fileStore.logToMaster("{}: sorting {}".format(config.uuid, vcf_file_name))
    vcf_file = os.path.join(work_dir, vcf_file_name)
    sorted_vcf_file = os.path.join(work_dir, "{}.sorted.vcf".format(vcf_file_name))
    header = list()
    content = list()
    # read input into memeory
    with open(vcf_file, 'r') as input:
        for line in input:
            if line.startswith("#"):
                header.append(line)
            else:
                content.append(line)
    #sort
    content.sort(key=lambda x: int(x.split("\t")[1]))
    # write to file
    with open(sorted_vcf_file, 'w') as output:
        for line in header:
            output.write(line)
        for line in content:
            output.write(line)
    # replace
    subprocess.check_call(["mv", sorted_vcf_file, vcf_file])


def _append_sam_reads_to_file(job, config, input_sam_file, output_sam_file, excluded_read_ids=set()):
    written_lines = 0
    read_ids_not_written = set()
    with open(output_sam_file, 'a') as output, open(input_sam_file, 'r') as input:
        for line in input:
            if line.startswith("@"): continue  # header
            read_id = line.split("\t")[0]
            if read_id in excluded_read_ids: # already written in a previous chunk
                read_ids_not_written.add(read_id)
                continue
            output.write(line)
            written_lines += 1
    job.fileStore.logToMaster("{}: wrote {} lines ({} excluded) from {} to {}"
                              .format(config.uuid, written_lines, len(read_ids_not_written),
                                      os.path.basename(input_sam_file), os.path.basename(output_sam_file)))
    return read_ids_not_written


def _append_vcf_calls_to_file(job, config, input_vcf_file, output_vcf_file, reverse_phasing=False):
    written_lines = 0
    with open(output_vcf_file, 'a') as output, open(input_vcf_file, 'r') as input:
        for line in input:
            if line.startswith("#"): continue  # header
            if reverse_phasing:
                line = line.rstrip().split("\t")
                phase = line[-1]
                has_bar = "|" in phase
                has_slash = "/" in phase
                if (not has_bar and not has_slash) or (has_bar and has_slash):
                    raise UserError("{}: Malformed vcf {} phasing line: {}"
                                    .format(config.uuid, os.path.basename(input_vcf_file), "\\t".join(line)))
                phase = phase.split("|") if has_bar else phase.split("/")
                phase.reverse()
                phase = ("|" if has_bar else "/").join(phase)
                line[-1] = phase
                line = "\t".join(line) + "\n"
            output.write(line)
            written_lines += 1
    job.fileStore.logToMaster("{}: wrote {} lines from {} to {}"
                              .format(config.uuid, written_lines,
                                      os.path.basename(input_vcf_file), os.path.basename(output_vcf_file)))
    return written_lines


def _extract_chunk_tarball(job, config, tar_work_dir, chunk):
    # prep
    os.mkdir(tar_work_dir)
    tar_file = os.path.join(tar_work_dir, "chunk.tar.gz")

    # get file
    job.fileStore.readGlobalFile(chunk[CI_OUTPUT_FILE_ID], tar_file, mutable=True)
    with tarfile.open(tar_file, 'r') as tar:
        tar.extractall(tar_work_dir)

    # find desired files
    sam_hap1, sam_hap2, vcf = None, None, None
    for name in os.listdir(tar_work_dir):
        if name.endswith(VCF_SUFFIX): vcf = name
        elif name.endswith(SAM_HAP_1_SUFFIX): sam_hap1 = name
        elif name.endswith(SAM_HAP_2_SUFFIX): sam_hap2 = name
    if sam_hap1 is None or sam_hap2 is None or vcf is None:
        raise UserError("{}: Missing expected output file, sam_hap1:{} sam_hap2:{} vcf:{} chunk_info:{}"
                        .format(config.uuid, sam_hap1, sam_hap2, vcf, chunk))
    sam_hap1_file = os.path.join(tar_work_dir, sam_hap1)
    sam_hap2_file = os.path.join(tar_work_dir, sam_hap2)
    vcf_file = os.path.join(tar_work_dir, vcf)

    # return file locations
    return sam_hap1_file, sam_hap2_file, vcf_file


def _get_read_ids_in_range(job, work_dir, file_name, contig_name, start_pos, end_pos):
    # samtools can't get random locations from sam files, so we convert to bam first :/
    bam_name = "{}.bam".format(file_name)
    bai_name = "{}.bai".format(bam_name)

    if not os.path.isfile(os.path.join(work_dir, bai_name)):
        # convert to bam
        convert_cmd = ["view", "-b", os.path.join(file_name), "-o", os.path.join(bam_name)]
        if DOCKER_LOGGING:
            job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, convert_cmd))
        dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=convert_cmd)

        # index
        index_cmd = ["index", os.path.join("/data", bam_name)]
        if DOCKER_LOGGING:
            job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, index_cmd))
        dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=index_cmd)

    # read_ids prep
    reads_filename = "{}.reads.txt".format(file_name)
    samtools_cmd = ["samtools", "view", os.path.join("/data", bam_name), "{}:{}-{}".format(contig_name, start_pos, end_pos)]
    column_script = [os.path.join("/data", _write_select_column_script(work_dir, 1))]
    tee_script = ["tee", os.path.join("/data", reads_filename)]

    # call docker
    params = [samtools_cmd, column_script, tee_script]
    if DOCKER_LOGGING:
        job.fileStore.logToMaster("Running {} with parameters: {}".format(DOCKER_SAMTOOLS, params))
    dockerCall(job, tool=DOCKER_SAMTOOLS, workDir=work_dir, parameters=params)

    # get output
    read_ids = set()
    with open(os.path.join(work_dir, reads_filename), 'r') as reads:
        for id in reads:
            read_ids.add(id.strip())
    return read_ids


def consolidate_output(job, config, chunk_infos):
    #prep
    start = time.time()
    work_dir = job.fileStore.getLocalTempDir()
    out_tar = os.path.join(work_dir, '{}.tar.gz'.format(config.uuid))
    job.fileStore.logToMaster("consolidate_output:{}:{}".format(config.uuid, datetime.datetime.now()))
    job.fileStore.logToMaster("consolidating {} files".format(len(chunk_infos)))

    # build tarball
    out_tars = [out_tar]
    with tarfile.open(out_tar, 'w:gz') as f_out:
        for ci in chunk_infos:
            file_id = ci[CI_OUTPUT_FILE_ID]
            tar_file = os.path.join(work_dir, "{}.tar.gz".format(ci[CI_CHUNK_INDEX]))
            job.fileStore.readGlobalFile(file_id, tar_file)
            out_tars.append(tar_file)
            with tarfile.open(tar_file, 'r') as f_in:
                for tarinfo in f_in:
                    with closing(f_in.extractfile(tarinfo)) as f_in_file:
                        f_out.addfile(tarinfo, fileobj=f_in_file)
    job.fileStore.logToMaster("{}: Consolidated {} tarballs".format(config.uuid, len(out_tars)))

    # Move to output location
    if urlparse(config.output_dir).scheme == 's3':
        job.fileStore.logToMaster('{}: Uploading {} to S3: {}'.format(config.uuid, out_tar, config.output_dir))
        s3am_upload(fpath=out_tar, s3_dir=config.output_dir, num_cores=config.cores)
    else:
        job.fileStore.logToMaster('{}: Moving {} to output dir: {}'.format(config.uuid, out_tar, config.output_dir))
        mkdir_p(config.output_dir)
        copy_files(file_paths=[out_tar], output_dir=config.output_dir)

    # log
    _log_time(job, "consolidate_output", start, config.uuid)
    job.fileStore.logToMaster("END:{}:{}".format(config.uuid, datetime.datetime.now()))



def _log_time(job, function_name, start_time, sample_identifier=''):
    job.fileStore.logToMaster("TIME:{}:{}:{}".format(sample_identifier, function_name, int(time.time() - start_time)))


def _get_default_docker_params(work_dir):
    return ['--rm','--log-driver','none', '-v', '{}:/data'.format(work_dir)]

def generate_config():
    return textwrap.dedent("""
        # MarginPhase Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in this configuration file and then rerun the pipeline: "toil-marginphase run"
        #
        # URLs can take the form: http://, ftp://, file://, s3://
        # Local inputs follow the URL convention: file:///full/path/to/input
        # S3 URLs follow the convention: s3://bucket/directory/file.txt
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false.
        ##############################################################################################################

        # Required: Output location of sample. Can be full path to a directory or an s3:// URL
        # Warning: S3 buckets must exist prior to upload or it will fail.
        # Warning: Do not use "file://" syntax if output directory is local location
        output-dir:

        # Required: URL {scheme} to reference VCF
        reference-vcf: s3://your/reference/here.vcf

        # Required: URL {scheme} to marginPhase params JSON file
        params: ftp://your/params/here.json

        # Required: Size of each bam partition
        partition-size: 2000000

        # Required: Margin to apply on each partition
        partition-margin: 5000

    """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample to be run.
        #
        #   There are 4 tab-separated columns: filetype, paired/single, UUID, URL(s) to sample
        #
        #   UUID            This should be a unique identifier for the sample to be processed.  Must belong to a single contig
        #   URL             A URL {scheme} pointing to the sample or a full path to a directory
        #   CONTIG_NAME     Contig name
        #   REFERENCE_URL   A URL {scheme} pointing to reference fasta file
        #
        #   If a full path to a directory is provided for a sample, every file inside needs to be a fastq(.gz).
        #   Do not have other superfluous files / directories inside or the pipeline will complain.
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/file.bam    chr3    file:///path/to/chr3.reference.fa
        #   UUID_2  s3://path/to/file.bam   chrX    s3://path/to/chrX.reference.fa
        #   UUID_3  s3://path/to/file.bam   4   file:///path/to/chr4.reference.fa
        #
        #   Place your samples below, one per line.
        """.format(scheme=[x + '://' for x in SCHEMES])[1:])


def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message

    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    MarginPhase pipeline

    =======================================
    Dependencies
    Curl:       apt-get install curl
    Docker:     wget -qO- https://get.docker.com/ | sh
    Toil:       pip install toil
    Boto:       pip install boto (OPTIONAL)
    """

    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')

    # Generate subparsers
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')

    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the MarginPhase pipeline')
    group = parser_run.add_mutually_exclusive_group()
    parser_run.add_argument('--config', default=DEFAULT_CONFIG_NAME, type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    group.add_argument('--manifest', default=DEFAULT_MANIFEST_NAME, type=str,
                       help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                            '\nDefault value: "%(default)s"')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Add Toil options
    Job.Runner.addToilOptions(parser_run)
    args = parser.parse_args()

    # Parse subparsers related to generation of config and manifest
    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_CONFIG_NAME), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, DEFAULT_MANIFEST_NAME), generate_manifest)

    # Pipeline execution
    elif args.command == 'run':
        # sanity check
        require(os.path.exists(args.config), '{} not found. Please run '
                                             '"toil-marginphase generate-config"'.format(args.config))
        require(os.path.exists(args.manifest), '{} not found and no samples provided. Please '
                                               'run "toil-marginphase generate-manifest"'.format(args.manifest))

        # get samples
        samples = parse_samples(path_to_manifest=args.manifest)

        # Parse config
        parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
        config = argparse.Namespace(**parsed_config)
        config.maxCores = int(args.maxCores) if args.maxCores else sys.maxsize
        config.maxDisk = int(args.maxDisk) if args.maxDisk else sys.maxint
        config.maxMemory = args.maxMemory if args.maxMemory else str(sys.maxint)

        # Config sanity checks
        require(config.output_dir, 'No output location specified')
        if urlparse(config.output_dir).scheme != "s3":
            mkdir_p(config.output_dir)
        if not config.output_dir.endswith('/'):
            config.output_dir += '/'
        # require(config.reference_contig, 'No reference contig specified')
        # require(config.contig_name, 'No contig name specified')
        require(config.reference_vcf, 'No reference vcf specified')
        require(config.params, 'No params specified')
        require(config.partition_size, "Configuration parameter partition_size is required")
        require(config.partition_margin, "Configuration parameter partition_margin is required")

        # Program checks
        for program in ['docker']:
            require(next(which(program), None), program + ' must be installed on every node.'.format(program))

        # Start the workflow
        Job.Runner.startToil(Job.wrapJobFn(map_job, prepare_input, samples, config), args)


if __name__ == '__main__':
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
