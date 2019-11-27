# Making Chain Files

makeChainFile.R is an R wrapper for creating chain files for usage
with liftover to convert genomic coordinates from one genomic build 
to another build.

Example of usage: Creating a chain file to convert genomeA coordinates
to genomeB coordinates.

`Rscript makeChainFile.R -o genomeA.fa - t genomeB.fa -n genomeA_to_genomeB`

Help is shown using the argument `--help`

`Rscript makeChainFile.R --help`

```
Usage: makeChainFile.R [options]


Options:
        -o CHARACTER, --originGenome=CHARACTER
                Origin genome in FASTA format

        -t CHARACTER, --targetGenome=CHARACTER
                Target genome in FASTA format

        -n CHARACTER, --fileName=CHARACTER
                Output filename

        -m INTEGER, --minIdentity=INTEGER
                minIdentity argument for BLAT step [default= 95]

        -d CHARACTER, --outputDir=CHARACTER
                Output directory

        -v LOGICAL, --verbose=LOGICAL
                Tells you what is going on [default= TRUE]

        -h, --help
                Show this help message and exit
```

References: 
* https://iamphioxus.org/2013/06/25/using-liftover-to-convert-genome-assembly-coordinates/
* https://hgwdev.gi.ucsc.edu/~kent/src/unzipped/hg/doc/liftOver.txt
* http://genomewiki.ucsc.edu/index.php/LiftOver_Howto
* http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver

