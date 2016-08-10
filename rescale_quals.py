#!/usr/bin/python
# Rescales phred 64 to phred 33. 
# TODO Change MAPQ to be the basequal somehow
# Adapted from @matted (https://www.biostars.org/p/51005/)
import sys

def process(fname):
    fin = open(fname)
    for line in fin:
        if line[0] == '@': # print header lines:
            print line,
            continue
        line = line.strip().split("\t")
        if len(line[10]) != len(line[9]):
            print >>sys.stderr, "number of quality scores does not match number of bases: [%s] vs. [%s]" % (line[9], line[10])
            line[10] = 'O' * len(line[9])
        line[10] = "".join([chr(ord(c)-31) for c in line[10]]) # ASCII-64 to 33 encoding of Phred vals
        print "\t".join(line)


    fin.close()

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print >>sys.stderr, "usage (converts fastq 1.3/1.5 to sanger): samtools view -h in.bam | bam_rescale_quals.py - | samtools view -bS - > out.bam"
    else:
        if sys.argv[1] == "-":
            sys.argv[1] = "/dev/stdin"
        process(sys.argv[1])

