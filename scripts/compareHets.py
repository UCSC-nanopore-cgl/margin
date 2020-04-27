import sys

if len(sys.argv) != 3:
    print("compareHets.py trueHetsFile predictedHetsFile")
    sys.exit(0)

trueHetsFile = sys.argv[1]
predictedHetsFile = sys.argv[2]


def rle(s):
    # Makes RLE string
    return [s[i] for i in range(len(s)) if i == 0 or s[i] != s[i - 1]]


def compareRLEs(str1, str2):
    return rle(str1[5:-5]) != rle(str2[5:-5])


def rc(s):
    c = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([c[i] for i in s[::-1]])


def getMismatches(mismatchFile):
    mismatches = {}
    mismatchLines = []
    with open(mismatchFile) as fh:
        for line in fh:
            tokens = line.split()
            if len(tokens) > 0 and tokens[0] == "Mismatch":
                mismatchLines.append((tokens[-2], tokens[-1], int(tokens[1]), int(tokens[2])))
                minStr = lambda x: x if x < rc(x) else rc(x)
                # tokens[-1] = minStr(tokens[-1])
                # tokens[-2] = minStr(tokens[-2])
                if tokens[-1] > tokens[-2]:
                    mismatches[tokens[-2]] = (tokens[-1], int(tokens[1]), int(tokens[2]))
                else:
                    mismatches[tokens[-1]] = (tokens[-2], int(tokens[1]), int(tokens[2]))
    return mismatches, mismatchLines

trueMismatches, trueMismatchLines = getMismatches(trueHetsFile)
predictedMismatches, predictedMismatchLines = getMismatches(predictedHetsFile)

totalCommonMismatches = set(trueMismatches.keys()).intersection(set(predictedMismatches.keys()))

print("Total common mismatches %s, total true mismatches %s, total predicted mismatches %s" % \
      (len(totalCommonMismatches), len(trueMismatches), len(predictedMismatches)))

print("Predicted mismatches")
for mismatch1, mismatch2, locationH1, locationH2 in predictedMismatchLines:
    print("Predicted mismatch: %s %s %s %s %s %s" % (
        mismatch1, mismatch2, locationH1, locationH2, mismatch1 in trueMismatches or mismatch2 in trueMismatches,
        compareRLEs(mismatch1, mismatch2)))

print("True mismatches")
for mismatch1, mismatch2, locationH1, locationH2 in trueMismatchLines:
    print("True mismatch: %s %s %s %s %s %s" % (
        mismatch1, mismatch2, locationH1, locationH2,
        mismatch1 in predictedMismatches or mismatch2 in predictedMismatches,
        compareRLEs(mismatch1, mismatch2)))
