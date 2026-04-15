#!/usr/bin/env python
# Title: Align all FASTQ files in DIR, sort and index
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos
# Dependencies: STAR, bedtools, argparse, subprocess, os
import argparse
import subprocess
import os

# Arguments
parser = argparse.ArgumentParser(description=
"****************************************************************************\n"
"STAR alignment with default parameters for all files in specified directory,"
" compressed FASTA files with R1(/R2).fastq.gz extensions are expected\n"
"****************************************************************************")
parser.add_argument("cores", type = int, help="Number of cores used (Default = 2)", default=2)
parser.add_argument("genomeDir", help="Reference genome directory")
parser.add_argument("-v", "--verbose", action="store_true", help="Prints what is going on")
parser.add_argument("-w", "--workdir", action = "store", help = "Input directory",
                    default = os.getcwd())
parser.add_argument("-o", "--outdir", action = "store", help = "Output directory",
                    default = os.getcwd())
args = parser.parse_args()

# Main program
iniDir = os.getcwd()
dirFiles = os.listdir(iniDir)
R1files = [f for f in dirFiles if "R1.fastq.gz" in f]

for file in R1files:
    read1 = file
    read2 = file.replace("R1.fastq.gz", "R2.fastq.gz")
    outPreffix = file.replace("R1.fastq.gz", "")
    if args.verbose:
        print("{} alignment started...".format(outPreffix))
    subprocess.call("STAR --runMode alignReads "
                    "--runThreadN {} "
                    "--genomeDir {} "
                    "--readFilesCommand zcat "
                    "--readFilesIn {} {} "
                    "--outFileNamePrefix {} "
                    "--outSAMtype BAM Unsorted"
                    "".format(args.cores, args.genomeDir, read1, read2,
                              outPreffix), shell = True)
    if args.verbose:
        print("... alignment done!")
if args.verbose:
    print("All alignments done!")

# Make summary log
if args.verbose:
    print("Making alignment summary log...")
subprocess.call("Rscript mapEfficiency_summary.R", shell = True)

# Remove all temporary and unnecessary files
if args.verbose:
    print("Removing temporary files...")
dirFiles = os.listdir(iniDir)
garbage = [f for f in dirFiles if "SJ.out.tab" in f] + \
          [f for f in dirFiles if "Log.progress.out" in f] + \
          [f for f in dirFiles if "Log.out" in f] + \
          [f for f in dirFiles if "Log.final.out" in f]
for file in garbage:
    os.remove(file)
tmpFolders = [f for f in dirFiles if "STARtmp" in f]
for folder in tmpFolders:
    subprocess.call("rm -r {}".format(folder), shell = True)

# Sort and index with samtools
if args.verbose:
    print("Sorting and indexing...")
bamFiles = [f for f in dirFiles if ".bam" in f]
for file in bamFiles:
    subprocess.call("samtools sort -o {} {}".format(
        file.replace(".bam", ".sorted.bam"), file), shell = True)
    subprocess.call("samtools index {}".format(
        file.replace(".bam", ".sorted.bam")), shell = True)
    os.remove(file)

# Finish
if args.verbose:
    print("DONE!")