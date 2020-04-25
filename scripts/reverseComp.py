import sys

""" Quick and dirty make reverse complement of a fasta file"""


def rc(s):
    c = {"A": "T", "C": "G", "G": "C", "T": "A", "a": "t", "c": "g", "g": "c", "t": "a"}
    return "".join([(c[i] if i in c else i) for i in s[::-1]])


with open(sys.argv[1]) as fh:
    read = []
    for line in fh:
        if ">" in line:
            print(line[:-1])
            if read != []:
                print(rc("".join(read))[1:] + "\n")
                read = []
        elif line != "\n":
            read.append(line)
    if read != []:
        print(rc("".join(read))[1:])
