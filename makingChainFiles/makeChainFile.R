#!/usr/bin/env Rscript
# Title: Make chain file
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############
# Based on Jia-Xing Yue post https://iamphioxus.org/2013/06/25/using-liftover-to-convert-genome-assembly-coordinates/
library("optparse")

# Parsing Arguments ############################################################
option_list = list(
    make_option(c("-o", "--originGenome"), type="character", default = NULL,
                help="Origin genome in FASTA format", metavar="character"),
    make_option(c("-t", "--targetGenome"), type="character", default = NULL,
                help="Target genome in FASTA format", metavar="character"),
    make_option(c("-n", "--fileName"), type="character", default = NULL,
                help="Output filename", metavar="character"),
    make_option(c("-m", "--minIdentity"), type="integer", default = 95,
                help="minIdentity argument for BLAT step [default= %default]", metavar="integer"),
    make_option(c("-d", "--outputDir"), type="character", default = "currentDirectory",
                help="Output directory  [default= %default]", metavar="character"),
    make_option(c("-c", "--nCores"), type="integer", default = 1,
                help="Number of cores to be used [default= %default]", metavar="integer"),
    make_option(c("-v", "--verbose"), type = "logical", default = T,
                help="Tells you what is going on [default= %default]", metavar="logical")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

originGenome <- opt$originGenome
targetGenome <- opt$targetGenome
fileName <- opt$fileName
minIdentity <- opt$minIdentity
outputDir <- opt$outputDir
verbose <- opt$verbose
nCores <- opt$nCores

# Requiring output file name
if(is.null(fileName)){
    stop("Output Filename Required!")
}

#Functions #####################################################################
# source("http://bit.ly/rnaMods")
listFilePatt <- function(pattern, path = "."){
    files <- list.files(path)[grep(pattern = pattern, x = list.files(path))]
    return(files)
}

# Packages #####################################################################
library(parallel)

# MAIN program #################################################################
ORdir <- getwd()

# Append path to file names if target and/or origin genomes are in working directory
if(sum(list.files() %in% originGenome) == 1){
    originGenome <- file.path(ORdir, originGenome)
}
if(sum(list.files() %in% targetGenome) == 1){
    targetGenome <- file.path(ORdir, targetGenome)
}
if(outputDir == "currentDirectory"){
    outputDir <- ORdir
}

# Creating txt  file for report
dummy <- file.create(file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")))

# Move to temporary directory to work in "background"
BGdir <- gsub("^/", "", tempdir())
BGdir <- file.path(outputDir, BGdir)
dir.create(BGdir, recursive = TRUE)
setwd(BGdir)

# Splitting chromosomes from target genome 
write(x = "Splitting chromosomes...", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
system(paste0("/usr/bin/csplit -s -z ", targetGenome, " '/>/' '{*}'"))
system('for i in xx* ; do   n=$(sed "s/>// ; s/ .*// ; 1q" "$i") ;   mv "$i" "$n.fa" ; done')

# Split the G2 chromosomes/scaffolds into 3K chunks and make lift files
dir.create("lift")
dir.create("split")
chroms <- listFilePatt(pattern = ".fa")
for (i in chroms){
    system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/faSplit size ", i, " 3000 ./split/", i,
                  ".split -lift=./lift/", i, ".lft -oneFile"))
}

# run blat
write(x = "Aligning using BLAT...", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
dir.create("psl")

cl <- makeCluster(nCores, type = "FORK")
dummy <- parSapply(cl, chroms, function(i) {
  write(x = paste("Aligning", i),
          file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
  system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/blat ", originGenome,
                  " ./split/", i, 
                  ".split.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=", minIdentity, " -noHead -minScore=100 ./psl/",
                  i, ".psl"))
  })
stopCluster(cl)

# Change coordinates of .psl files to parent coordinate system
dir.create("liftup")
for(i in chroms){
    system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/liftUp -pslQ ./liftup/",
                  i, ".liftup.psl ./lift/", i, ".lft warn ./psl/", i, ".psl"))
}

# Make chain files
write(x = "Making chain files...", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
dir.create("chain_raw")
for(i in chroms){
    system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/axtChain -linearGap=medium -faQ -faT -psl ./liftup/",
                  i, ".liftup.psl ", originGenome, " ", targetGenome, " ./chain_raw/", i, ".chain"))
}

# Merge and sort chain files
write(x = "Merging chain files...", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
system("/apps/RH7U2/general/kentUtils/v377/bin/chainMergeSort ./chain_raw/*.chain | chainSplit chain_split stdin")
system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/faSize ", originGenome, " -detailed > G1.chr_length.txt"))
system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/faSize ", targetGenome, " -detailed > G2.chr_length.txt"))

# Make alignment nets from chain files
dir.create("net")
for(i in listFilePatt("chain", "chain_split/")){
    system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/chainNet ", "chain_split/", i, " ./G1.chr_length.txt ./G2.chr_length.txt ./net/", i, ".net /dev/null"))
}

# Create liftOver chain file
write(x = "Creating liftOver chain file", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
dir.create("over")
for(i in listFilePatt("chain", "chain_split/")){
    system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/netChainSubset ./net/", i, ".net ", "chain_split/", i, " ./over/", i, ".chain"))
}

# Output File
if(!grepl(pattern = ".chain$", fileName)){
    outFileName <- paste0(fileName, ".chain")
}else{outFileName <- fileName}

system(paste0("cat ./over/*.chain > ", file.path(outputDir, outFileName)))
write(x = "...DONE.", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)

# Erasing temporary folder
system(paste0("rm ", BGdir, " -r"))
