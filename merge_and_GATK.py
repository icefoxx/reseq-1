#!/usr/bin/env python
import sys, os
import argparse
import textwrap

args=argparse.ArgumentParser(description='***This program runs the bowtie alignment pipeline on a dataset as specified in listfile.')
args.add_argument('inputfile', help='Specify the file that lists all of the read data to be aligned')
args.add_argument('inputdir', help='Specify the input files that have been deduped, etc. from previous pipeline')
args.add_argument('output',help='Specify the output directory for the alignment files')
theinput=args.parse_args()

inputfile = theinput.inputfile
inputdir = theinput.inputdir.rstrip('/')
outputdir = theinput.output.rstrip('/')

pools_dict={}


with open(inputfile) as fd:
  for count,line in enumerate(fd):
    columns=line.strip().split("\t")
    samplename=columns[0]
    info=columns[1]
    path,filename=os.path.split(info)
    base_filename,ext =os.path.splitext(filename)
    species,junk,pool_id,sampleid,read_num=base_filename.split('_')

    if pools_dict.has_key(samplename)==False:
        pools_dict[samplename]=[pool_id]
    else:
        pools_dict[samplename]+= [pool_id]


for k,v in pools_dict.iteritems():
    inputline=''
    for item in v:
        inputline+=' INPUT={2}/{0}_{1}.deduped.sorted.reordered.bam'.format(k,item,inputdir)
    voorline='qsub -V -N {0}_merge -j y -o {1}/{0}.merge.log -pe omp 2 -b y java -Djava.io.tmpdir=/projectnb/mullenl/tmp -jar /projectnb/mullenl/programs/picard-tools-1.74/MergeSamFiles.jar'.format(k,outputdir)
    outputline=' OUTPUT={1}/{0}.deduped.merged.bam'.format(k,outputdir)
    command = voorline + inputline + outputline + "\n"
    print command
    os.system(command)

