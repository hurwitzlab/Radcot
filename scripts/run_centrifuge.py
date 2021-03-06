#!/usr/bin/env python3
"""Run Centrifuge"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile as tmp
from pprint import pprint

#WORK env var will be present on TACC
#But may not be set when testing locally
if os.getenv('WORK') is None:
    os.environ['WORK'] = './'

# --------------------------------------------------
def get_args():
    """get args"""
    parser = argparse.ArgumentParser(
        description='Options for run_centrifuge.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-q', '--query',
                        help='File or directory of input',
                        metavar='str',
                        type=str,
                        nargs='+',
                        default='')

    parser.add_argument('-1', '--forward',
                        help='Comma-separated forward reads\n'
                        'that pair with reverse reads',
                        metavar='str',
                        type=str,
                        default='')

    parser.add_argument('-2', '--reverse',
                        help='Comma-separated reverse reads\n'
                        'that pair with forward reads',
                        metavar='str',
                        type=str,
                        default='')

    parser.add_argument('-U', '--unpaired',
                        help='Comma-separated unpaired reads',
                        metavar='str',
                        type=str,
                        default='')

    parser.add_argument('-f', '--format',
                        help='Input file format',
                        metavar='str',
                        type=str,
                        default='fasta')

    parser.add_argument('-i', '--index',
                        help='Centrifuge index name',
                        metavar='str',
                        type=str,
                        default='p_compressed+h+v')

    parser.add_argument('-I', '--index_dir',
                        help='Centrifuge index directory',
                        metavar='str',
                        type=str,
                        default='/work/05066/imicrobe/iplantc.org/data/centrifuge-indexes')

    parser.add_argument('-o', '--out_dir',
                        help='Output directory',
                        metavar='str',
                        type=str,
                        default=os.path.join(os.getcwd(), 'centrifuge-out'))

    parser.add_argument('-g', '--genome-dir',
                        dest='genome_dir', metavar='DIRECTORY',
                        default=os.path.join(os.getenv('WORK'),'genomes'),
                        help="Directory with all the genomes (*.fna's).\n"
                        "This is important because it must be the same in\n"
                        "all steps. Thus, if you keep the default here, you must\n"
                        "keep the default in all steps\n"
                        " [ Default = $WORK/genomes ]")

    parser.add_argument('-x', '--exclude_taxids',
                        help='Exclude tax ids',
                        metavar='str',
                        type=str,
                        default='')

    parser.add_argument('-X', '--max_seqs_per_file',
                        help='Maxium num of sequences per file',
                        metavar='int',
                        type=int,
                        default=1000000)

    parser.add_argument('-t', '--threads',
                        help='Num of threads per instance of centrifuge',
                        metavar='int',
                        type=int,
                        default=1)

    parser.add_argument('-P', '--procs',
                        help='Max number of processes to run',
                        metavar='int',
                        type=int,
                        default=1)
    
    parser.add_argument("-m", "--min_abundance", action="store", \
                        help="Minimum abundance needed to download a species\' genome", \
                        default='.01', type=float)

    parser.add_argument("-a", "--annotation_type", \
                        default="refseq", choices=['refseq','patric'], \
                        help="Which type of annotation to get: refseq or patric. " \
                        "NOTE: Choosing \"refseq\" will discard genomes that only have PATRIC annotaion. " \
                        "NOTE2: patric annotations will generally include refseq genes. " \
                        "Default is \"refseq\" as it is generally better curated.")

    return parser.parse_args()

# --------------------------------------------------
def warn(msg):
    """Print a message to STDERR"""
    print(msg, file=sys.stderr)

# --------------------------------------------------
def die(msg='Something went wrong'):
    """Print a message to STDERR and exit with error"""
    warn('Error: {}'.format(msg))
    sys.exit(1)

# --------------------------------------------------
def find_input_files(query):
    """Find input files from list of files/dirs"""
    files = []
    for qry in query:
        if os.path.isdir(qry):
            for filename in os.scandir(qry):
                if filename.is_file():
                    files.append(filename.path)
        elif os.path.isfile(qry):
            files.append(qry)
        else:
            die('--query "{}" neither file nor directory'.format(qry))
    return files

# --------------------------------------------------
def line_count(fname):
    """Count the number of lines in a file"""
    n = 0
    for _ in open(fname):
        n += 1

    return n

# --------------------------------------------------
def run_job_file(jobfile, msg='Running job', procs=1):
    """Run a job file if there are jobs"""
    num_jobs = line_count(jobfile)
    warn('{} (# jobs = {})'.format(msg, num_jobs))
#--halt now,fail=1 causes parallel to exit with code 1 if any subjobs exit with code 1
    if num_jobs > 0:
        cmd = 'parallel --halt now,fail=1 -P {} < {}'.format(procs, jobfile)
        warn(cmd)
        completed_process = subprocess.run(cmd, shell=True)
    else:
        warn('No jobs to run!')
        return 0

    os.remove(jobfile)

    return completed_process.returncode

# --------------------------------------------------
def split_files(out_dir, files, max_seqs, file_format, procs):
    """Split input files by max_sequences"""
    split_dir = os.path.join(out_dir, 'split')

    if not os.path.isdir(split_dir):
        os.makedirs(split_dir)

    jobfile = tmp.NamedTemporaryFile(delete=False, mode='wt')
    bin_dir = os.path.dirname(os.path.realpath(__file__))
    tmpl = '{}/fxsplit.py -i {} -f {} -o {} -n {}\n'

    for input_file in files:
        out_file = os.path.join(split_dir, os.path.basename(input_file))
        if not os.path.isdir(out_file):
            jobfile.write(tmpl.format(bin_dir,
                                      input_file,
                                      file_format,
                                      out_file,
                                      max_seqs))


    jobfile.close()

    return_code = run_job_file(jobfile=jobfile.name, msg='Splitting input files', procs=procs)

    if return_code != 0:
        die()

    return list(filter(os.path.isfile,
                       glob.iglob(split_dir + '/**', recursive=True)))

# --------------------------------------------------
def run_centrifuge(file_format, files, exclude_ids, 
        index_name, index_dir, out_dir, threads, procs):
    """Run Centrifuge"""
    reports_dir = os.path.join(out_dir, 'reports')

    if not os.path.isdir(reports_dir):
        os.makedirs(reports_dir)

    jobfile = tmp.NamedTemporaryFile(delete=False, mode='wt')
    exclude_arg = '--exclude-taxids ' + exclude_ids if exclude_ids else ''
    if file_format == 'fasta':
        tmpl = 'CENTRIFUGE_INDEXES={} centrifuge {} -f -p {} -x {} -U {} -S {} --report-file {}\n'
    elif file_format == 'fastq':
        tmpl = 'CENTRIFUGE_INDEXES={} centrifuge {} -p {} -x {} -U {} -S {} --report-file {}\n'
    else:
        die('Need to specifiy read format for centrifuge to work')

    for file in files:
        basename = os.path.basename(file)
        tsv_file = os.path.join(reports_dir, basename + '.tsv')
        sum_file = os.path.join(reports_dir, basename + '.sum')
        if not os.path.isfile(tsv_file):
            jobfile.write(tmpl.format(index_dir,
                                      exclude_arg,
                                      threads,
                                      index_name,
                                      file,
                                      sum_file,
                                      tsv_file))
    jobfile.close()

    return_code = run_job_file(jobfile=jobfile.name, msg='Running Centrifuge', procs=procs)
    
    if return_code != 0:
        die()

    return list(filter(os.path.isfile,
                       glob.iglob(reports_dir + '/**', recursive=True)))

def run_cent_paired(file_format, paired_reads, exclude_ids, 
        index_name, index_dir, out_dir, threads, procs):
    """Run Centrifuge"""
    reports_dir = os.path.join(out_dir, 'reports')

    if not os.path.isdir(reports_dir):
        os.makedirs(reports_dir)

    jobfile = tmp.NamedTemporaryFile(delete=False, mode='wt')
    exclude_arg = '--exclude-taxids ' + exclude_ids if exclude_ids else ''
    if file_format == 'fasta':
        tmpl = 'CENTRIFUGE_INDEXES={} centrifuge {} -f -p {} -x {} -1 {} -2 {} -S {} --report-file {}\n'
    elif file_format == 'fastq':
        tmpl = 'CENTRIFUGE_INDEXES={} centrifuge {} -p {} -x {} -1 {} -2 {} -S {} --report-file {}\n'
    else:
        die('Need to specifiy read format for centrifuge to work')

    for f, r in paired_reads:
        basename = os.path.basename(f)
        tsv_file = os.path.join(reports_dir, basename + '.tsv')
        sum_file = os.path.join(reports_dir, basename + '.sum')
        if not os.path.isfile(tsv_file):
            jobfile.write(tmpl.format(index_dir,
                                      exclude_arg,
                                      threads,
                                      index_name,
                                      f,
                                      r,
                                      sum_file,
                                      tsv_file))
    jobfile.close()

    return_code = run_job_file(jobfile=jobfile.name, msg='Running Centrifuge', procs=procs)

    if return_code != 0:
        die()

    return list(filter(os.path.isfile,
                       glob.iglob(reports_dir + '/**', recursive=True)))

def get_genomes(reports_dir, genome_dir, min_abundance, annotation_type, procs):
    """Get genomes from PATRIC"""

    if not os.path.isdir(genome_dir):
        os.makedirs(genome_dir)
    
    jobfile = tmp.NamedTemporaryFile(delete=False, mode='wt')
    bin_dir = os.path.dirname(os.path.realpath(__file__))
    tmpl='{}/cfuge_to_genome.py --report {} --output {} --min_abundance {} --annotation_type {}\n'

    for report in glob.iglob(reports_dir + '/*.tsv'):
        jobfile.write(tmpl.format(bin_dir,
                                  report,
                                  genome_dir,
                                  min_abundance,
                                  annotation_type))

    jobfile.close()

    return_code = run_job_file(jobfile=jobfile.name, msg='Getting genomes', procs=procs)
    
    if return_code != 0:
        die()


# --------------------------------------------------
def get_excluded_tax(ids):
    """Verify the ids look like numbers"""
    tax_ids = []

    if ids:
        for s in [x.strip() for x in ids.split(',')]:
            if s.isnumeric():
                tax_ids.append(s)
            else:
                warn('tax_id "{}" is not numeric'.format(s))

    return ','.join(tax_ids)

# --------------------------------------------------
def collapse_reports(input_files, reports, out_dir):
    """Collapse the split reports"""
    collapse_dir = os.path.join(out_dir, 'collapsed')
    if not os.path.isdir(collapse_dir):
        os.makedirs(collapse_dir)

    for i, filename in enumerate(input_files):
        basename = os.path.basename(filename)
        print('{:4}: {}'.format(i + 1, basename))
        basename, ext = os.path.splitext(basename)
        splits = {'tsv': [], 'sum': []}

        for report in reports:
            regex = r'{}\.(\d+)\.{}\.(tsv|sum)'.format(basename, ext[1:])
            match = re.match(regex, os.path.basename(report))

            if match:
                num, type_ext = match.groups()
                # storing a tuple with the num first for sorting
                splits[type_ext].append((int(num), report))

        for file_type in splits:
            # sort on the fst of the tuple but take the snd
            files = [a[1] for a in \
                     sorted(splits[file_type], key=lambda a: a[0])]

            if len(files) < 1:
                msg = 'WARNING: No files ending with "{}" for "{}"'
                print(msg.format(file_type, basename))
            else:
                out_path = os.path.join(collapse_dir,
                                        basename + '.' + file_type)

                if os.path.isfile(out_path):
                    print('      "{}" exists, skipping'.format(out_path))
                else:
                    print('      Writing to "{}"'.format(out_path))
                    func = write_tsv if file_type == 'tsv' else write_sum
                    func(files, open(out_path, 'w'))


    return collapse_dir

# --------------------------------------------------
def write_tsv(files, out_fh):
    """collapse tsv files"""
    tax = dict()
    num_flds = ['numReads', 'numUniqueReads']

    for fnum, file in enumerate(files):
        print('    {:4}: {}'.format(fnum + 1, os.path.basename(file)))
        with open(file) as csvfile:
            reader = csv.DictReader(csvfile, delimiter='\t')
            for row in reader:
                tax_id = row['taxID']
                if not tax_id in tax:
                    tax[tax_id] = dict()
                    for fld in row.keys():
                        tax[tax_id][fld] = int(row[fld]) \
                                if fld in num_flds else row[fld]
                else:
                    for fld in num_flds:
                        tax[tax_id][fld] += int(row[fld])

    # Write the headers
    flds = ['name', 'taxID', 'taxRank', 'genomeSize'] + \
            num_flds + ['abundance']
    out_fh.write("\t".join(flds) + '\n')

    total_reads = sum([tax[s]['numReads'] for s in tax])

    for tax_id in sorted(tax.keys(), key=int):
        species = tax[tax_id]

        species['abundance'] = round(species['numReads'] / total_reads, 2)

        out_fh.write('\t'.join([str(species[f]) for f in flds]) + '\n')

# --------------------------------------------------
def write_sum(files, out_fh):
    """collapse sum files"""
    for fnum, file in enumerate(files):
        print('    {:4}: {}'.format(fnum + 1, os.path.basename(file)))
        in_fh = open(file, 'r')
        hdr = in_fh.readline()
        if fnum == 0:
            out_fh.write(hdr)

        for line in in_fh:
            out_fh.write(line)

# --------------------------------------------------
def make_bubble(collapse_dir, out_dir):
    """Make bubble chart"""
    fig_dir = os.path.join(out_dir, 'figures')

    if not os.path.isdir(fig_dir):
        os.makedirs(fig_dir)

    cur_dir = os.path.dirname(os.path.realpath(__file__))
    bubble = os.path.join(cur_dir, 'centrifuge_bubble.r')
    job = '{} --dir {} --outdir {}'.format(bubble,
                                           collapse_dir,
                                           fig_dir)

    subprocess.run(job, shell=True)

    return fig_dir

# --------------------------------------------------
def main():
    """main"""
    args = get_args()
    out_dir = args.out_dir
    index_dir = args.index_dir
    index_name = args.index

    if not index_dir:
        print('--index_dir is required')
        sys.exit(1)

    if not index_name:
        print('--index_name is required')
        sys.exit(1)

    valid_index = set(['p_compressed', 'p_compressed+h+v', 'p+h+v']) #since this is getting genomes from PATRIC db of bacteria, we do not want nt
    if not index_name in valid_index:
        tmpl = '--index "{}" is not valid, please choose from: {}'
        die(tmpl.format(index_name, ', '.join(sorted(valid_index))))

    if not os.path.isdir(index_dir):
        die('--index_dir "{}" is not a directory'.format(index_dir))

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    exclude_ids = get_excluded_tax(args.exclude_taxids)

    if args.query or args.unpaired:
        if args.query and not args.unpaired:
            input_files = find_input_files(args.query)
        elif args.unpaired and not args.query:
            input_files = args.unpaired.split(',')
        else:
            die('Using query input AND unpaired is unsupported\n' +
                    'Please choose one or the other')
        
        num_files = len(input_files)
        warn('Found {} input file{}'.format(num_files,
                                            '' if num_files == 1 else 's'))

        if num_files == 0:
            die('No usable files from --query')

        split_file_names = split_files(out_dir=out_dir,
                                       files=input_files,
                                       file_format=args.format,
                                       max_seqs=args.max_seqs_per_file,
                                       procs=args.procs)

        reports = run_centrifuge(file_format=args.format,
                                 files=split_file_names,
                                 out_dir=out_dir,
                                 exclude_ids=exclude_ids,
                                 index_dir=index_dir,
                                 index_name=index_name,
                                 threads=args.threads,
                                 procs=args.procs)

        collapse_dir = collapse_reports(input_files=input_files,
                                        reports=reports,
                                        out_dir=out_dir)

        fig_dir = make_bubble(collapse_dir=collapse_dir, out_dir=out_dir)

        print('Done, reports in "{}", figures in "{}"'.format(collapse_dir,
                                                              fig_dir))
    elif args.forward and args.reverse:

        f_list = args.forward.split(',')
        r_list = args.reverse.split(',')

        if len(f_list) != len(r_list):
            die('Number of forward reads is not the same as reverse reads!!')

        num_files = len(f_list)

        warn('Found {} paired read{}'.format(num_files,
                                            '' if num_files == 1 else 's'))

        if num_files == 0:
            die('No usable files from -1/-2')

        paired = [] #list of tuples (forward_read,reverse_read)
        for read in f_list:
            paired.append((read, r_list[f_list.index(read)]))

        #NOTE: reports get named after the forward reads so that also has to be the input
        #for "collapsing"
        reports = run_cent_paired(file_format=args.format,
                                 paired_reads=paired,
                                 out_dir=out_dir,
                                 exclude_ids=exclude_ids,
                                 index_dir=index_dir,
                                 index_name=index_name,
                                 threads=args.threads,
                                 procs=args.procs)

        reports_dir = os.path.dirname(reports[0])
        warn('Reports dir: {}'.format(reports_dir))
        fig_dir = make_bubble(collapse_dir=reports_dir, out_dir=out_dir)

#        print('Done, reports in "{}", figures in "{}"'.format(reports_dir,
#                                                              fig_dir))

    else:
        die('Need a query or paired forward and reverse reads\n' +
                'This is what is have: query {}\n'.format(args.query) +
                'forward {}\n'.format(args.forward) +
                'reverse {}\n'.format(args.reverse) +
                'unpaired {}\n'.format(args.unpaired))
    #getting genomes using cfuge_to_genome.py 
    min_abundance = args.min_abundance
    annotation_type = args.annotation_type
    genome_dir = args.genome_dir
    get_genomes(reports_dir=reports_dir,
                            genome_dir=genome_dir,
                            min_abundance=min_abundance,
                            annotation_type=annotation_type,
                            procs=args.procs)

        
# --------------------------------------------------
if __name__ == '__main__':
    main()
