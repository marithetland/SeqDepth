#!/usr/bin/env python3

import os, sys, re
import logging, time
import glob
import csv
import shutil
import datetime
from argparse import ArgumentParser
from subprocess import call
from subprocess import call
from subprocess import check_output, CalledProcessError, STDOUT

# Exception classes

def parse_args():
    #Version
    parser = ArgumentParser(description='SeqDepth')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + "v.1.0.0")

    #Argsgroups
    input_args = parser.add_argument_group('Input options (required)')
    optional_args = parser.add_argument_group('Optional flags')
    output_args = parser.add_argument_group('Output options')

    #Input
    input_args.add_argument('-r', '--reads', nargs='+', type=str, required=True, help='Provide the full file path to the (trimmed) read reads (ending in e.g. *fastq.gz, fq.gz)')
    input_args.add_argument('-a', '--assemblies', nargs='+', type=str, required=True, help='Provide the full file path to the assembly-files (ending in *fasta)')
    #Options
    optional_args.add_argument('-k', '--keep', action='store_true', required=False, help='Use this flag if you want to keep the BAM files used in the calculation.')
    #Output
    output_args.add_argument('-o', '--output', type=str, required=False, default='', help='Optional outfile name. Default="seqdepth_results.csv"')

    return parser.parse_args()

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)   

class CommandError(Exception): #yekwahs
    pass

def run_command(command, **kwargs): #yekwahs
    command_str = ''.join(command)
    #logging.info('Running shell command: {}'.format(command_str))
    try:
        exit_status = call(command_str, **kwargs)
    except OSError as e:
        message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
        raise CommandError({"Error:": message})
    if exit_status != 0:
        message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
        raise CommandError({"Error:": message})

def file_exists(seqlist, program, path, extention):
    re_run = []
    for seq in seqlist:
        if os.path.exists(path + seq + extention):
            echo = ("File-extention " + extention + " already exists for " + seq)
        else:
            echo = ("No such file present. Running " + program + " now on " + seq)
            re_run.append(seq)
    if re_run:
        return re_run

def convert(list): 
    ''' Convert a list to a string '''
    s = [str(i) for i in list] 
    res = (",".join(s)) 
    return(res)

def main():
    #Setup 
    print("SeqDepth: Calculate the sequencing depth (overall coverage) and standard deviation of trimmed FASTQ-files. \n")
    args = parse_args()

    logging.basicConfig(
        level=logging.DEBUG,
        filemode='w',
        format='%(asctime)s %(message)s',
        datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('Running SeqDepth v.1.0.0')
    logging.info('Command line: {0}'.format(' '.join(sys.argv)))
    #Sort arguments

    reads=args.reads
    fasta=args.assemblies
    current_dir = os.getcwd() + '/'
    #Create output directories
    if args.output and not os.path.exists(args.output):
        outfile=current_dir+args.output
    elif args.output and os.path.exists(args.output):
        logging.info('Output file "' + args.output + '" already exists. Will append to this file.')
        outfile=current_dir+args.output
    elif not args.output:
        outfile=(current_dir+'seqdepth_results.txt')


#Check input
    sequence_list = []
    for sequence in reads:
        name = sequence.replace(".gz","")
        if name.find('_1.fastq') != -1 and name[:-8] not in sequence_list:
            sequence_list.append(name[:-8])
        if name.find('_2.fastq') != -1 and name[:-8] not in sequence_list:
            sequence_list.append(name[:-8])  
    for sequence in fasta:
        name = sequence.replace("","")
        if name.find('.fasta') != -1 and name[:-6] not in sequence_list:
            sequence_list.append(name[:-6])  

    #Check that all reads have pairs
    missing_pairs = []
    for seqName in sequence_list:
        sequence_1 = False
        sequence_2 = False
        fasta_f = False
        for sequence in reads:
            if sequence.find(seqName+'_1.fastq.gz') != -1:
                sequence_1 = True
            if sequence.find(seqName+'_2.fastq.gz') != -1:
                sequence_2 = True
        for sequence in fasta:
            if sequence.find(seqName+'.fasta') != -1:
                fasta_f = True
        if sequence_1 == False or sequence_2 == False or fasta_f == False:
            missing_pairs.append(seqName)

    if missing_pairs != []:
        logging.info("\nNot all read sequence sets have pairs:")
        for seq in missing_pairs:
            print(seq)
        logging.info("Pipeline Stopped: please fix sequence pairs.")
        sys.exit()
    else:
        logging.info("\nAll read-sets matched. Sequences to be processed:")
        print(sequence_list)

#Run program          
    for item in sequence_list: 
        logging.info(item)
        if not os.path.isfile(item+"_Coverage.Success"):
            try:
                fasta_file=(item+'.fasta')
                read_1=(item+'_1.fastq.gz')
                read_2=(item+'_2.fastq.gz')
                run_command(['echo "eee ',fasta_file,'"'], shell= True)
                run_command(['bwa index ', fasta_file], shell= True)
                run_command(['bwa mem -t 8 ',fasta_file,' ',read_1,' ',read_2,' > input_c.sam  ; \
                    picard SamFormatConverter INPUT=input_c.sam VALIDATION_STRINGENCY=SILENT OUTPUT=input_c.bam ; \
                    picard SortSam INPUT=input_c.bam OUTPUT=input_2_c.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate ; \
                    picard MarkDuplicates INPUT=input_2_c.bam VALIDATION_STRINGENCY=SILENT OUTPUT=final_cont.bam METRICS_FILE=dup_metrics ; \
                    picard BuildBamIndex INPUT=final_cont.bam VALIDATION_STRINGENCY=SILENT OUTPUT=final_cont.bam.bai ' ], shell= True)
                run_command(["echo -n '",item," \t' >> ",outfile," ; \
                    tot_size=$(samtools view -H final_cont.bam | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}') ;\
                    echo $tot_size ; \
                    samtools depth final_cont.bam | awk -v var=$tot_size '{sum+=$3; sumsq+=$3*$3} END {print sum/var \"\t\" sqrt(sumsq/var - (sum/var)**2)}' >> ",outfile," \
                    ;  touch ",item,"_Coverage.Success"  ], shell= True)
                if args.keep:
                    for bamfile in glob.glob('*final_cont.bam'):
                        shutil.move(bamfile, item+'_'+bamfile)
                else:        
                    for bamfile in glob.glob('*final_cont.bam'):
                        os.remove(bamfile)
                for tempfile in glob.glob('dup_metrics'):
                    os.remove(tempfile)
                for tempfile in glob.glob('input*'):
                    os.remove(tempfile)
                for tempfile in glob.glob('*.fasta.*'):
                    os.remove(tempfile)
                for tempfile in glob.glob('final_cont.bam.bai'):
                    os.remove(tempfile)
                if os.path.isfile(item+"_Coverage.Success"):
                    logging.info(item+": Coverage calculation success.")
            except:
                logging.info(item+": Coverage-calculation unsuccessful.")
        else:
            logging.info(item+": Coverage calculation already performed. If this is wrong, make sure you do not have a file for your sample ending in '.Success' in the working directory.")
    

if __name__ == '__main__':
    main()
