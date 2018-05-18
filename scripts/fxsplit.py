#!/usr/bin/env python3
"""split FASTA/Q files"""

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>
# Second Author: Scott G Daniel <scottdaniel@email.arizona.edu>
import gzip
import sys
import argparse
import os
from Bio import SeqIO

# --------------------------------------------------
def get_args():
    """get args"""
    parser = argparse.ArgumentParser(
        description="Split FASTA files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-i", "--infile", help="Input file",
                        type=str, metavar="FILE", required=True)

    parser.add_argument("-f", "--format", help="Format (fasta, fastq)",
                        type=str, metavar="FILE", default="fasta")

    parser.add_argument("-n", "--num", help="Number of records per file",
                        type=int, metavar="NUM", default=50)

    parser.add_argument("-o", "--out_dir", help="Output directory",
                        type=str, metavar="DIR", default="split-files")

    return parser.parse_args()

# --------------------------------------------------
def main():
    args = get_args()
    infile = args.infile
    file_format = args.format.lower()
    out_dir = args.out_dir
    max_per = args.num

    if not os.path.isfile(infile):
        print('--infile "{}" is not valid'.format(infile))
        sys.exit(1)

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if max_per < 1:
        print("--num cannot be less than one")
        sys.exit(1)

    if not file_format in set(['fasta', 'fastq']):
        print("--format ({}) must be fasta/q".format(file_format))
        sys.exit(1)

    i = 0
    nseq = 0
    nfile = 0
    out_fh = None
    basename, ext = os.path.splitext(os.path.basename(infile))

    handle = None
    if ext == ".gz":
        handle = gzip.open(infile, "rt")
        basename, ext = os.path.splitext(basename)
    else:
        handle = open(infile, "rt")

    record_iter = SeqIO.parse(handle, file_format)
    for i, batch in enumerate(batch_iterator(record_iter, max_per)):
        filename = os.path.join(out_dir, basename + '.' + str(i+1) + ext)
        with open(filename, "w") as handle:
            count = SeqIO.write(batch, handle, file_format)
        print("Wrote {:d} records to {:s}".format(count, filename))

# from biopython.org/wiki/Split_large_file
# --------------------------------------------------
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


# --------------------------------------------------
if __name__ == '__main__':
    main()
