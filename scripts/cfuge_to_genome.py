#!/usr/bin/env python3
'''
Script that takes a centrifuge report and download genomes and annotations
'''
import glob
import time
import sys
import argparse
import os
import errno
from plumbum import local
import pandas as pd

#make programs easily accessible in python using plumbum!
p3_all_genomes = local['p3-all-genomes']
#by default we do not need to download things twice!
wget = local['wget']['--no-clobber']
wc = local['wc']

#digest those arguments
if __name__ == "__main__":
    parser = \
    argparse.ArgumentParser(description="Script that takes centrifuge report(s) and download(s) genomes and annotations")
    
    parser.add_argument("-r", "--report", action="store", \
            help="Directory containing centrifuge reports or single centrifuge report", \
            default='./')
    parser.add_argument("-o", "--output", action="store", \
            help="Output directory for genomes and annotations", \
            default='./')
    parser.add_argument("-m", "--min_abundance", action="store", \
            help="Minimum abundance needed to download a species\' genome", \
            default='.01', type=float)
    parser.add_argument("-a", "--annotation_type", \
            default="refseq", choices=['refseq','patric'], \
            help="Which type of annotation to get: refseq or patric. " \
            "NOTE: Choosing \"refseq\" will discard genomes that only have PATRIC annotaion. " \
            "NOTE2: patric annotations will generally include refseq genes. " \
            "Default is \"refseq\" as it is generally better curated.")
        
    args = parser.parse_args()

def filter_report(report):
    
    holder = []

    for row in report.itertuples(index=True, name='Pandas'):
        if getattr(row, 'abundance') > args.min_abundance:
#            print("{:d}".format(getattr(row, 'taxID')))
            holder.append(getattr(row, 'taxID'))
    
    print("After filtering, these are the NCBI genome id's: {}".format(holder))
    return holder

def download_genomes(filtered_list):
    
    #change working directory
    try:
        os.chdir(args.output)
    except FileNotFoundError as e:
        print("{}, creating {} for you".format(e, args.output))
        os.mkdir(args.output)
        os.chdir(args.output)

    print("Now trying to get the PATRIC id's from p3.theseed.org")
    for taxID in filtered_list:
        for i in range(0,3):
            try:
                print("Attempt {} of 3".format(i+1))
                list_of_patricIDs = p3_all_genomes('-e', 'taxon_id'+','+str(taxID))
    #                '-e', 'genome_status,Complete')
                list_of_patricIDs = list_of_patricIDs.split('\n') 
                list_of_patricIDs = [x for x in list_of_patricIDs if x != '' and x != 'genome.genome_id']
                print("These are the PATRIC id's I got from {}: {}".format(taxID,list_of_patricIDs))
            
            except Exception as e:
                if i == 2:
                    print("patricbrc.org or p3.theseed.org might be down")
                    print("Contact patricbrc.org")
                    print("Exiting...")
                    sys.exit(1)
                else:
                    print("Something went wrong with the PATRIC cli. Error: {}".format(e))
                    print("Can not continue without PATRIC cli, retrying after 5 seconds")
                    time.sleep(5)
                    continue

            break

        for patricID in list_of_patricIDs:
            if patricID != '' and patricID != 'genome.genome_id':
                if args.annotation_type == 'refseq':
                    print('Getting PATRIC genome_id {}'.format(patricID))

                    try:
                        print('Getting the genome nucleotide sequence ".fna"')
                        wget('ftp://ftp.patricbrc.org/genomes/' + patricID + '/' + patricID + '.fna')
                    except Exception as e:
                        print("Something went wrong with downloading the .fna. Error {}".format(e))
                        
                    try:
                        print('Getting all the RefSeq annotation for {}'.format(patricID))
                        wget('ftp://ftp.patricbrc.org/genomes/' + patricID + '/' + \
                            patricID + '.RefSeq*')
                    except Exception as e:
                        print("Something went wrong with downloading the Refseq annotation Error {}".format(e))

                    #check to make sure its not a bogus refseq.gff
                    if os.path.isfile(patricID + '.RefSeq.gff'):
                        line_count = int(wc('-l', patricID + '.RefSeq.gff').strip().split()[0])
                        if line_count < 5:
                            print('The RefSeq.gff is less than 5 lines, it is probably bogus.')
                            [os.remove(x) for x in glob.glob(patricID + '*')]
    #                        os.remove(patricID + '.RefSeq*')
    #                        os.remove(patricID + '.fna' )
                    else:
                        print('There was no RefSeq.gff for genome_id {}, therfore deleting the genome'.format(patricID))
                        [os.remove(x) for x in glob.glob(patricID + '*')]

                elif args.annotation_type == 'patric':
                    print('Getting PATRIC genome_id {}'.format(patricID))

                    try:
                        print('Getting the genome nucleotide sequence ".fna"')
                        wget('ftp://ftp.patricbrc.org/genomes/' + patricID + '/' + patricID + '.fna')
                    except Exception as e:
                        print("Something went wrong with downloading the .fna. Error {}".format(e))
                     
                    try:
                        print('Getting all the PATRIC annotation for {}'.format(patricID))
                        wget('ftp://ftp.patricbrc.org/genomes/' + patricID + '/' + \
                        patricID + '.PATRIC*')
                    except Exception as e:
                        print("Something went wrong with downloading the PATRIC annotation Error {}".format(e))

                else:
                    print("The annotation type {} does not compute".format(args.annotation_type))
                    sys.exit(1)
                    
def get_reports(report_file_or_dir):
    
    if os.path.isfile(args.report):

        report = pd.read_table(args.report,delimiter='\t')
        list_of_taxids = filter_report(report)
        download_genomes(list_of_taxids)

    elif os.path.isdir(args.report):
        tsv_files = []
        for root, _, filenames in os.walk(args.report):
            for filename in filenames:
                if 'tsv' == filename.split('.')[-1]:
                    tsv_files.append(os.path.join(root, filename))
            if len(tsv_files) < 1:
                print('Found no files in --report {}'.format(args.report))
                sys.exit(1)
            else:
                for tsv in tsv_files:
                    report = pd.read_table(tsv,delimiter='\t')
                    list_of_taxids = filter_report(report)
                    download_genomes(list_of_taxids)

    else:
        print('{} does not seem to be a file or a directory!'.format(args.report))
        sys.exit(1)

#All the program besides the functions and setup
print("Start! {:s}".format(time.ctime()))

get_reports(args.report)

print("Done! {:s}".format(time.ctime()))
