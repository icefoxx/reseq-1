#!/usr/bin/env python
 
import sys, os
import argparse
 
args=argparse.ArgumentParser(description='***This program runs the bowtie alignment pipeline on a dataset as specified in listfile.')
args.add_argument('inputfile', help='Specify the file that lists all of the read data to be aligned')
args.add_argument('output',help='Specify the output directory for the alignment files')
args.add_argument('reference', help='Specify the reference file that the data will be aligned to')
args.add_argument('illuminaid',help='ID Number for Illumina Run (e.g. CDH7ACXX)')
theinput=args.parse_args()
 
inputfile = theinput.inputfile
reference = os.path.abspath(theinput.reference)
refbase=os.path.basename(reference)
refbase=os.path.splitext(refbase)[0]
outputdir = theinput.output
outputdir = outputdir.rstrip('/')
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
#print(therefpath)
root,dirs,files=os.walk(therefpath).next()
for f in files:
	if f.endswith('.bt2'):
		if f.startswith(refbase):
			bowtie2indexcounter+=1
 
if bowtie2indexcounter < 6:
	bowtierefname=os.path.join(therefpath,refbase)
	bowtie2indexcommand="""$BT2_HOME/bowtie2-build %s %s"""%(reference,bowtierefname)
	os.system(bowtie2indexcommand)
 
with open(inputfile) as fd:
	for count,line in enumerate(fd):
		columns=line.strip().split("\t")
		name=os.path.basename(columns[0])
		name=name.split(".")
		sample=name[0].split("_")
		samplename=sample[1]
		read1=columns[0]
		read2=columns[1]
		unpaired=columns[2]
		#print name, read1,read2,unpaired
		#commandline="""qsub -l h_rt=08:00:00 -V -N %s -b y $BT2_HOME/bowtie2 -x %s --fast-local -I 0 -X 1000 -1 %s -2 %s -U %s -S %s.sam"""%(name,reference,read1,read2,unpaired,name)
		qsubscript= '''
#! /bin/bash
#$ -V
#$ -N %s_%s
#$ -j y
#$ -l h_rt=08:00:00
#$ -pe omp 2
$BT2_HOME/bowtie2 -x %s --fast-local -I 0 -X 1000 -1 %s -2 %s -U %s | samtools view -bS - > %s_%s.bam
java -jar /projectnb/mullenl/programs/picard-tools-1.74/AddOrReplaceReadGroups.jar INPUT=%s_%s.bam OUTPUT=%s_%s.RGadded.bam ID=%s PL=ILLUMINA PU=CDHP7ACXX LB=1 SM=%s
java -jar /projectnb/mullenl/programs/picard-tools-1.74/ReorderSam.jar INPUT=%s_%s.RGadded.bam OUTPUT=%s_%s.reordered.bam REFERENCE=%s
samtools sort %s_%s.reordered.bam %s_%s.sorted
samtools index %s_%s.sorted 
'''%(samplename,refbase,refbase,read1,read2,unpaired,samplename,refbase,samplename,refbase,samplename,refbase,samplename,samplename,samplename,refbase,samplename,refbase,reference,samplename,refbase,samplename,refbase,samplename,refbase)
 
		#print qsubscript
		#if count > 1: break
 
		qsfname = samplename + "_" + refbase + "_qsub" + ".sh"
		qsub = open(qsfname, 'w')
		qsub.write(qsubscript)
		qsub.close()
		qsubscript = ''  #reset the script
 
	# submit it to SGE with qsub
		cmd = "qsub " + qsfname
		os.system(cmd)