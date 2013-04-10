#!/usr/bin/env python
import sys, os
import argparse
import textwrap

args=argparse.ArgumentParser(description='***This program runs the bowtie alignment pipeline on a dataset as specified in listfile.')
args.add_argument('inputfile', help='Specify the file that lists all of the read data to be aligned')
args.add_argument('output',help='Specify the output directory for the alignment files')
args.add_argument('reference', help='Specify the reference file that the data will be aligned to')
theinput=args.parse_args()

inputfile = theinput.inputfile
reference = os.path.abspath(theinput.reference)
refbase=os.path.basename(reference)
refbase=os.path.splitext(refbase)[0]
outputdir = theinput.output.rstrip('/')
count= 1
bowtie2indexcounter=0

#check for REFERENCE file dictionary
refpath=os.path.splitext(reference)[0]
dictfile=refpath + '.dict'

if os.path.exists(dictfile) is False:
  createdict="""java -jar /projectnb/mullenl/programs/picard-tools-1.74/CreateSequenceDictionary.jar REFERENCE=%s O=%s""" %(reference,dictfile)
  os.system(createdict)

#check for BOWTIE2 INDEX
therefpath=os.path.dirname(reference)
print(therefpath)
root,dirs,files=os.walk(therefpath).next()
for f in files:
  if f.endswith('.bt2'):
    bowtie2indexcounter+=1

if bowtie2indexcounter < 6:
  bowtierefname=os.path.join(therefpath,refbase)
  bowtie2indexcommand="""$BT2_HOME/bowtie2-build %s %s"""%(reference,bowtierefname)
  os.system(bowtie2indexcommand)

with open(inputfile) as fd:
  for count,line in enumerate(fd):
    columns=line.strip().split("\t")
    samplename=columns[0]
    info=columns[1]
    path,filename=os.path.split(info)
    base_filename,ext =os.path.splitext(filename)
    species,junk,pool_id,sampleid,read_num=base_filename.split('_')

    with open(columns[1], 'r') as f:
        first_line = f.readline()
        headerdata=first_line.split(':')
        flowcell=headerdata[2]
        lane=headerdata[3]
        f.close()

    rg_info = """\
    --rg SM:{0} \
    --rg LB:{0} \
    --rg ID:{1}.{2}.{0} \
    --rg PL:ILLUMINA \
    --rg CN:HARVARDFAS \
    --rg PU:{1}.{2} \
    --rg PI:500 \
    --rg-id ID:{1}.{2}.{0} """.format(samplename, flowcell, lane)
    rg_info = textwrap.dedent(rg_info)

    forward_reads=path+"/"+species+"_"+junk+"_"+pool_id+"_"+sampleid+"_"+"R1.adapter-trimmed.fastq"
    reverse_reads=path+"/"+species+"_"+junk+"_"+pool_id+"_"+sampleid+"_"+"R2.adapter-trimmed.fastq"
    unpaired_reads=path+"/"+species+"_"+junk+"_"+pool_id+"_"+sampleid+"_"+"R1.adapter-trimmed.read-singleton.fastq"
    name = sampleid+"_"+pool_id

    qsubscript= '''
#! /bin/bash
#$ -V
#$ -N {0}
#$ -j y
#$ -o {5}/{0}.log
#$ -l h_rt=08:00:00
#$ -pe omp 4
$BT2_HOME/bowtie2 -x {1} -p 3 --very-sensitive-local -I 0 -X 1000 -1 {2} -2 {3} -U {4} {6} | samtools view -bS - > {5}/{0}.bam
java -Xmx8g -Djava.io.tmpdir=/scratch/ -jar /projectnb/mullenl/programs/picard-tools-1.74/ReorderSam.jar INPUT={5}/{0}.bam OUTPUT={5}/{0}.reordered.bam REFERENCE={1}.fasta
java -Xmx8g -Djava.io.tmpdir=/scratch/ -jar /projectnb/mullenl/programs/picard-tools-1.74/SortSam.jar INPUT={5}/{0}.reordered.bam OUTPUT={5}/{0}.sorted.reordered.bam SORT_ORDER=coordinate
java -Xmx8g -Djava.io.tmpdir=/scratch/ -jar /projectnb/mullenl/programs/picard-tools-1.74/MarkDuplicates.jar INPUT={5}/{0}.sorted.reordered.bam OUTPUT={5}/{0}.deduped.sorted.reordered.bam METRICS_FILE={5}/{0}.metrics.txt CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500 ASSUME_SORTED=true 
'''.format(name,refbase,forward_reads,reverse_reads,unpaired_reads,outputdir,rg_info)

   
    print qsubscript
    #if count > 1: break
    
    qsfname = outputdir+"/"+ name + "_qsub" + ".sh"
    qsub = open(qsfname, 'w')
    qsub.write(qsubscript)
    qsub.close()
    qsubscript = ''  #reset the script

  # submit it to SGE with qsub
    cmd = "qsub " + qsfname
    os.system(cmd)

