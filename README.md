# SchwartzLab_ToolBox

A set of programs, scripts, and tutorials that make our life easier.

## Tutorials

A collection of notebooks meant to instruct how to use the main features of
packages useful for genomic research.

* **Bioconductor Genomic Packages**: An introduction and example of usage to 
*GenomicRanges*, *GenomicAlignments*, and
*Gviz* packages.

## ggPlotShortcuts.R

A collection of functions that resemble base R plots syntax but outputs ggplots!

## makeChainFile.R

A command line tool that wraps Kent command line utilities for creating chain 
files, which are needed to convert genomic coordinates from one genomic build to
another using [liftover](https://genome-store.ucsc.edu/).

## STAR_alignAll

A command line tool that wraps STAR and some R code to align all fastq data
present in the working directory (more for illustrative purposes and as a 
template, that you will need to adapt to your needs).
