#!/usr/bin/env Rscript
# Title: 
# Author: Miguel Angel Garcia-Campos https://github.com/AngelCampos ############
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
                help="Output directory", metavar="character"),
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

# Requiring output file name
if(is.null(fileName)){
    stop("Output Filename Required!")
}

# Testing
#originGenome <- "/home/labs/schwartzlab/miguelg/BIGDATA/genome_references/mouse_strains/CAST_EiJ/completeGenome/CAST_EiJ_genome.fa"
#targetGenome <- "/home/labs/schwartzlab/miguelg/BIGDATA/genome_references/mm9/mm9_genome.fa"
#fileName <- "CAST_mm9"
#minIdentity <- 95
#outputDir <- "currentDirectory"
#verbose <- T

#Functions #####################################################################
# source("http://bit.ly/rnaMods")
listFilePatt <- function(pattern, path = "."){
    files <- list.files(path)[grep(pattern = pattern, x = list.files(path))]
    return(files)
}

# Packages #####################################################################
library(parallel)
#CRAN packages
# CRAN_packs <- c("magrittr")
# dummy <- sapply(CRAN_packs, function(x) installLoad_CRAN(x))

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
# BGdir <- tempdir()
BGdir <- file.path(outputDir, "tmp")
dir.create(BGdir)
setwd(BGdir)

# Splitting chromosomes from target genome 
system(paste0("/usr/bin/csplit -s -z ", targetGenome, " '/>/' '{*}'"))
system('for i in xx* ; do   n=$(sed "s/>// ; s/ .*// ; 1q" "$i") ;   mv "$i" "$n.fa" ; done')
write(x = "Splitting chromosomes...", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)

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

cl <- makeCluster(detectCores() - 1, type = "FORK")
parSapply(cl, chroms, function(i) {
  write(x = paste("Aligning", i),
          file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
  system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/blat ", originGenome,
                  " ./split/", i, 
                  ".split.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=95 -noHead -minScore=100 ./psl/",
                  i, ".psl"))
  })
stopCluster(cl)

# for(i in chroms){
#     write(x = paste("Aligning", i),
#           file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)
#     system(paste0("/apps/RH7U2/general/kentUtils/v377/bin/blat ", originGenome,
#                   " ./split/", i, 
#                   ".split.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=95 -noHead -minScore=100 ./psl/",
#                   i, ".psl"))
#     
# }

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
if(!grepl(pattern = ".over.chain$", fileName)){
    outFileName <- paste0(fileName, ".over.chain")
}else{outFileName <- fileName}

system(paste0("cat ./over/*.chain > ", file.path(outputDir, outFileName)))
write(x = "...DONE.", file = file.path(outputDir, paste0("makeChainReport_", fileName, ".txt")), append = T)

# Erasing temporary folder
# system(paste0("rm ", BGdir, " -r"))
