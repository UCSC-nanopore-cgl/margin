#!/usr/bin/env python3
# from __future__ import print_function
import math
import os
import subprocess
import sys


def main():
    trueAssembly = open(sys.argv[1]).readlines()
    assembly = sys.argv[2]
    verbose = True if len(sys.argv) == 4 and sys.argv[-1] == "verbose" else False

    read = []
    for line in trueAssembly[1:]:
        if ">" in line:
            break
        read.append(line[:-1])
    read = "".join(read)

    # This spike in was used to confirm that the x sequence that is printed by the code  is the true assembly 
    # and the y sequence is the assembly
    # so x gaps are insertions with respect to the true assembly
    # and y gaps are deletions with respect to the true assembly
    # read = read[:30000] + "GATTACAGATTACA" + read[30000:]

    print("Read len:", len(read), file=sys.stderr)
    # read = read[:70000]
    fh = open("temp.fa", 'w')
    fh.write(">hello\n%s\n" % read)
    fh.close()

    subprocess.call("lastz temp.fa %s --format=axt > temp.axt" % (assembly,), shell=True)

    axt = open("temp.axt", "r").readlines()
    axt = [i for i in axt if i[0] != '#']

    identities = []
    totalMatches = 0
    totalXLen = 0
    totalYLen = 0
    totalXGaps = 0
    totalYGaps = 0
    totalXGapEvents = 0
    totalYGapEvents = 0
    totalMismatches = 0
    for j in range(len(axt)):
        if "hello" in axt[j]:

            trueAssemblyAlignmentStart = int(axt[j].split()[2])
            predictedAssemblyAlignmentStart = int(axt[j].split()[5])
            x = axt[j + 1].upper()
            y = axt[j + 2].upper()

            def getAlignmentToSequenceCoordinateMap(alignedSeq):
                toSeqCoordinate = {}
                i = 0
                for k in range(len(alignedSeq)):
                    toSeqCoordinate[k] = i
                    if alignedSeq[k] != '-':
                        i += 1
                return toSeqCoordinate

            xAlignmenToSeqCoordinate = getAlignmentToSequenceCoordinateMap(x)
            yAlignmenToSeqCoordinate = getAlignmentToSequenceCoordinateMap(y)

            xLen = len([i for i in x if i != '-'])
            yLen = len([i for i in y if i != '-'])

            if xLen > 50000 and yLen > 50000:
                # x = x[10000:-10000]
                # y = y[10000:-10000]
                # xLen = len([i for i in x if i != '-'])
                # yLen = len([i for i in y if i != '-'])

                xGaps = len(x) - xLen
                yGaps = len(y) - yLen
                xGapEvents = len([i for i in range(1, len(x)) if x[i] == '-' and x[i - 1] != '-'])
                yGapEvents = len([i for i in range(1, len(y)) if y[i] == '-' and y[i - 1] != '-'])
                mismatches = sum([1 if (x[i] != y[i] and x[i] != '-' and y[i] != '-') else 0 for i in range(len(x))])
                mismatchLocations = [(trueAssemblyAlignmentStart + xAlignmenToSeqCoordinate[i] - 1,
                                      predictedAssemblyAlignmentStart + yAlignmenToSeqCoordinate[i] - 1,
                                      x[i - 10:i + 10], y[i - 10:i + 10]) for i in range(len(x)) if
                                     (x[i] != y[i] and x[i] != '-' and y[i] != '-')]
                matches = sum([1 if x[i] == y[i] else 0 for i in range(len(x))])

                identity = matches * 2 / float(xLen + yLen)

                identities.append((xLen, yLen, matches, identity, mismatches, xGaps, yGaps, mismatchLocations))

                totalMatches += matches
                totalXLen += xLen
                totalYLen += yLen
                totalXGaps += xGaps
                totalYGaps += yGaps
                totalXGapEvents += xGapEvents
                totalYGapEvents += yGapEvents
                totalMismatches += mismatches

    identities.sort()

    for xLen, yLen, matches, identity, mismatches, xGaps, yGaps, mismatchLocations in identities:
        print("xLen: %s, yLen: %s, matches: %s, identity: %s, mismatches: %s, xGaps: %s, yGaps: %s" % (
            xLen, yLen, matches, identity, mismatches, xGaps, yGaps), file=sys.stderr)
        if verbose:
            for trueAssemblyIndex, predictedAssemblyIndex, xString, yString in mismatchLocations:
                print(" Mismatch", trueAssemblyIndex, predictedAssemblyIndex, xString, yString)

    qv = lambda x: -10.0 * math.log10(1.0 - x) if x < 1.0 else 1000

    print(
        "total-xLen: %s, total-yLen: %s, total-matches: %s, total-identity: %s, "
        "qv: %s, qv-matches: %s, qv-inserts: %s, qv-deletes: %s, qv-indels: %s, "
        "total-mismatches: %s, total-xGaps (inserts): %s, total-yGaps (deletes): %s, total-xGapEvents (insert events): %s, total-yGapEvents (delete events): %s" %
        (totalXLen, totalYLen, totalMatches,
         totalMatches / (totalMatches + totalXGaps + totalYGaps + totalMismatches + 0.0001),  # identity
         qv(totalMatches / (totalMatches + totalXGaps + totalYGaps + totalMismatches + 0.0001)),  # overall qv
         qv(totalMatches / (totalMatches + totalMismatches + 0.0001)),  # mismatch qv
         qv(totalMatches / (totalXGaps + totalMatches + 0.0001)),  # insert qv
         qv(totalMatches / (totalYGaps + totalMatches + 0.0001)),  # delete qv
         qv(totalMatches / (totalMatches + totalXGaps + totalYGaps + 0.0001)),  # delete qv
         totalMismatches, totalXGaps, totalYGaps, totalXGapEvents, totalYGapEvents), file=sys.stderr)

    # Cleanup
    os.remove("temp.fa")
    os.remove("temp.axt")


if __name__ == '__main__':
    main()
