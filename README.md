# utilities
Collection of (hopefully) useful things

# descriptions

## ```screener``` & ```re-screen```

Submit jobs e.g. in a loop on a server without job scheduler.
- usage: ```screener <"command"> <job name>```
- creates directory .screen with files .screen/screen.name.submit and .screen/screen.name.out
- .screen/screen.name.submit contains submission time and command
- .screen/screen.name.out contains all STDERR and STDOUT
- example: screener "echo test" test
- make executable (chmod +x screener) and add do $PATH to call it from anywhere
- ```re-screen <.screen/screener.*.submit>``` re-launches an existing submission file from a screener job

  
## ```simplify_bam.py```

Simplify CIGAR string and remove tags in SAM/BAM file to reduce size and IGV load time; write read length to 9th field
- dependency: pysam (https://github.com/pysam-developers/pysam)
- reads input SAM/BAM
- simplifies input records by removing query alignment details and writing them to indexed output BAM
- removes tags to reduce size 
- can reduce BAM file size to <0.01 %
- usage: ```simplify_bam.py <input sam/bam> <output bam>```

                                    
## ```unscaffold.py```

Splits sequences in a FASTA at stretches of Ns of specified size and outputs a new FASTA
- Treats 'N' and 'n' the same
- does not modify gaps of other than the specified size
- optional size sort of output contigs
- reads all sequences into memory (--> therefore requires ~ 2x the memory size of the uncompressed input FASTA)
- usage: ```unscaffold.py <input path> <output path> <contig-prefix> <gap size> <out fasta line length> <size sort True/False>```
- for details check ```unscaffold.py --help```


