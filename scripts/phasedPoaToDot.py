#!/usr/bin/env python3
# from __future__ import print_function

""" Visualize POA file by converting to a dot file for graphviz """

import argparse
import math
import sys

import pandas as pd
from graphviz import Digraph

# Get user arguments
parser = argparse.ArgumentParser("Visualize a phased POA file from margin")
parser.add_argument("poaFile")
parser.add_argument('--start', '-s', dest='start', default=0, required=False, type=int,
                    help='First position in the ref to include (inclusive)')
parser.add_argument('--length', '-l', dest='length', default=sys.maxsize, required=False, type=int,
                    help='Length of POA subgraph backbone to visualize')
parser.add_argument('--rleCoordinate', '-r', dest='rleCoordinate', default=False, required=False, action="store_true",
                    help='Use RLE coordinates for specifying positions in POA backbone')
parser.add_argument('--output', '-o', dest='output', default='./phased_poa', required=False, type=str,
                    help='File to place PDF format POA graph viz in')
parser.add_argument('--show', '-w', dest='show', default=False, required=False, action="store_true",
                    help='Open the PDF output automatically')
args = parser.parse_args()

# Read poa into a data frame
df = pd.read_csv(args.poaFile, converters={'DELETES': str, 'INSERTS': str})

# Build graph
dot = Digraph(comment='Partial Order Alignment')
dot.attr(rankdir='LR')  # Force graph to have left to right layout

# Check for out of bounds
if args.start < 0:
    args.start = 0

# Convert start / length to rle coordinates
if not args.rleCoordinate:
    # Make a map of non-rle coordinates to rle coordinates
    toRleCoordinate = {}
    i = 0
    for refIndex, repeatCount in zip(df["REF_INDEX"][1:], df["REPEAT_COUNT"][1:]):
        for j in range(repeatCount):
            toRleCoordinate[i + j] = refIndex - 1
        i += repeatCount
    # Update start coordinate to be in rle space
    args.start = toRleCoordinate[args.start]

# Move to POA coordinate (with the 0 base being the prefix position)
args.start += 1

# Set the length of the POA interval
if args.length == sys.maxsize or args.start + args.length > df.shape[0]:
    args.length = df.shape[0] - args.start
end = args.start + args.length  # Exclusive last node

insertColor = "darkgreen"
backboneColor = "blue"
hap1Color = "red"
hap2Color = "purple"

# Make sink node for connecting deletions that delete the last position
dot.node(str(args.start - 1), "SOURCE")

# Make backbone nodes for each base
for i in range(args.start, end):
    row = df.iloc[i]

    labels = [
        "index:{}, base:{}, rc:{}, weight:{:.2f}".format(i, row["REF_BASE"], row["REPEAT_COUNT"], row["TOTAL_WEIGHT"])]
    labels.append(
        "h1:{:.2f} h2:{:.2f} h1+:{:.2f}, h2+{:.2f}\n".format(row["FRACTION_HAP1_WEIGHT"], row["FRACTION_HAP2_WEIGHT"],
                                                             row["FRACTION_POS_STRAND_HAP1"],
                                                             row["FRACTION_POS_STRAND_HAP2"]))

    ## Add base info
    for base in "ACTG":
        NORM_BASE_WEIGHT, FRACTION_BASE_HAP1, FRACTION_BASE_HAP2, FRACTION_BASE_POS_STRAND_HAP1, FRACTION_BASE_POS_STRAND_HAP2 = \
            row["FRACTION_BASE_%s_WEIGHT" % base], row["FRACTION_BASE_%s_HAP1" % base], row[
                "FRACTION_BASE_%s_HAP2" % base], \
            row["FRACTION_BASE_%s_POS_STRAND_HAP1" % base], row["FRACTION_BASE_%s_POS_STRAND_HAP2" % base]
        if NORM_BASE_WEIGHT >= 0.01:
            labels.append("{}:{:.2f}, h1:{:.2f}, h2:{:.2f},\n\th1+:{:.2f}, h2+:{:.2f}\n".format(base, NORM_BASE_WEIGHT,
                                                                                                FRACTION_BASE_HAP1,
                                                                                                FRACTION_BASE_HAP2,
                                                                                                FRACTION_BASE_POS_STRAND_HAP1,
                                                                                                FRACTION_BASE_POS_STRAND_HAP2))
    # labels.append("\n")

    ## Add repeat counts
    for repeatCount in range(1, 51):
        PROB_REPEAT_COUNT_HAP1 = row["PROB_HAP1_REPEAT_COUNT_%i" % repeatCount]
        PROB_REPEAT_COUNT_HAP2 = row["PROB_HAP2_REPEAT_COUNT_%i" % repeatCount]
        if PROB_REPEAT_COUNT_HAP1 > 0.001 or PROB_REPEAT_COUNT_HAP2 > 0.001:
            labels.append("repeat count:{}, h1:{:.3f}, h2:{:.3f}".format(repeatCount, PROB_REPEAT_COUNT_HAP1,
                                                                         PROB_REPEAT_COUNT_HAP2))

    dot.node(str(i), "\n".join(labels), shape="rectangle", fontcolor=backboneColor, color=backboneColor,
             penwidth=str(math.log(1 + row["TOTAL_WEIGHT"])))

# Make sink node for connecting deletions that delete the last position
dot.node(str(end), "END")

# Connect backbone nodes
for i in range(args.start - 1, end):
    row = df.iloc[i + 1]
    dot.edge(str(i), str(i + 1), penwidth=str(math.log(1 + row["TOTAL_WEIGHT"])))  ## Aggregate edge
    ## Haplotype specific edges
    # h1Weight = row["TOTAL_WEIGHT"] * row["FRACTION_HAP1_WEIGHT"]
    # h2Weight = row["TOTAL_WEIGHT"] * row["FRACTION_HAP2_WEIGHT"]
    # dot.edge(str(i), str(i+1), label="{:.2f}".format(h1Weight), color=hap1Color, penwidth=str(math.log(1 + h1Weight)))
    # dot.edge(str(i), str(i+1), label="{:.2f}".format(h2Weight), color=hap2Color, penwidth=str(math.log(1 + h2Weight)))

# Make deletions
for i in range(args.start - 1, end - 1):
    deletions = [j for j in df.iloc[i]["DELETES"].split("|") if j != '']
    for j in range(0, len(deletions), 6):
        DELETE_LENGTH, TOTAL_WEIGHT, FRACTION_HAP1_WEIGHT, FRACTION_HAP2_WEIGHT, FRACTION_POS_STRAND_HAP1, FRACTION_POS_STRAND_HAP2 = map(
            float, deletions[j:j + 6])

        # Non hap specific weights
        # dot.edge(str(i), str(i+int(DELETE_LENGTH)+1),
        #         label="weight:{:.2f}\nh1:{:.2f}, h2:{:.2f}\nh1+:{:.2f}, h2+:{:.2f}".format(TOTAL_WEIGHT, FRACTION_HAP1_WEIGHT,
        #                                                FRACTION_HAP2_WEIGHT, FRACTION_POS_STRAND_HAP1, FRACTION_POS_STRAND_HAP2),
        #         fontcolor=deleteColor, color=deleteColor, penwidth=str(math.log(1 + TOTAL_WEIGHT)))

        ## Haplotype specific edges
        h1Weight = TOTAL_WEIGHT * FRACTION_HAP1_WEIGHT
        h2Weight = TOTAL_WEIGHT * FRACTION_HAP2_WEIGHT

        if FRACTION_HAP1_WEIGHT > 0.01:
            dot.edge(str(i), str(i + int(DELETE_LENGTH) + 1),
                     label="weight:{:.2f}, %:{:.2f}, \n%+:{:.2f}".format(h1Weight, FRACTION_HAP1_WEIGHT,
                                                                         FRACTION_POS_STRAND_HAP1),
                     color=hap1Color, penwidth=str(math.log(1 + h1Weight)))

        if FRACTION_HAP2_WEIGHT > 0.01:
            dot.edge(str(i), str(i + int(DELETE_LENGTH) + 1),
                     label="weight:{:.2f}, %:{:.2f}, \n%+:{:.2f}".format(h2Weight, FRACTION_HAP2_WEIGHT,
                                                                         FRACTION_POS_STRAND_HAP2),
                     color=hap2Color, penwidth=str(math.log(1 + h2Weight)))

# Make insertions
for i in range(args.start - 1, end - 1):
    insertions = [j for j in df.iloc[i]["INSERTS"].split("|") if j != '']
    for j in range(0, len(insertions), 6):
        INSERT_SEQ = insertions[j]
        TOTAL_WEIGHT, FRACTION_HAP1_WEIGHT, FRACTION_HAP2_WEIGHT, FRACTION_POS_STRAND_HAP1, FRACTION_POS_STRAND_HAP2 = map(
            float, insertions[j + 1:j + 6])
        insertNodeName = "{},{}".format(i, j)
        labels = [INSERT_SEQ]
        labels.append(
            "weight:{:.2f}\nh1:{:.2f}, h2:{:.2f}\nh1+:{:.2f}, h2+:{:.2f}".format(TOTAL_WEIGHT, FRACTION_HAP1_WEIGHT,
                                                                                 FRACTION_HAP2_WEIGHT,
                                                                                 FRACTION_POS_STRAND_HAP1,
                                                                                 FRACTION_POS_STRAND_HAP2))
        penwidth = str(math.log(1 + TOTAL_WEIGHT))
        dot.node(insertNodeName, label="\n".join(labels), fontcolor=insertColor, color=insertColor, penwidth=penwidth)

        ## Haplotype specific edges

        h1Weight = TOTAL_WEIGHT * FRACTION_HAP1_WEIGHT
        h2Weight = TOTAL_WEIGHT * FRACTION_HAP2_WEIGHT

        if FRACTION_HAP1_WEIGHT > 0.01:
            penwidth = str(math.log(1 + h1Weight))
            dot.edge(str(i), insertNodeName, color=hap1Color, penwidth=penwidth,
                     label="weight:{:.2f}, %:{:.2f}, \n%+:{:.2f}".format(h1Weight, FRACTION_HAP1_WEIGHT,
                                                                         FRACTION_POS_STRAND_HAP1))
            dot.edge(insertNodeName, str(i + 1), color=hap1Color, penwidth=penwidth)

        if FRACTION_HAP2_WEIGHT > 0.01:
            penwidth = str(math.log(1 + h2Weight))
            dot.edge(str(i), insertNodeName, color=hap2Color, penwidth=penwidth,
                     label="weight:{:.2f}, %:{:.2f}, \n%+:{:.2f}".format(h2Weight, FRACTION_HAP2_WEIGHT,
                                                                         FRACTION_POS_STRAND_HAP2))
            dot.edge(insertNodeName, str(i + 1), color=hap2Color, penwidth=penwidth)

        # dot.edge(str(i), insertNodeName, color=insertColor, penwidth=penwidth)
        # dot.edge(insertNodeName, str(i+1), color=insertColor, penwidth=penwidth)

# Draw
dot.render(args.output, view=args.show)
