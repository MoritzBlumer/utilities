# utilities
Collection of (hopefully) useful things

# descriptions

## screener

Submit jobs e.g. in a loop on a server without job scheduler.
- usage: screener <"command"> <name>
- creates directory .screen with files .screen/screen.name.submit and .screen/screen.name.out
- .screen/screen.name.submit contains submission time and command
- .screen/screen.name.out contains all STDERR and STDOUT
- example: screener "echo test" test
- make executable (chmod +x screener) and add do $PATH to call it from anywhere

  
## simplify_bam.py

Simplify CIGAR string and remove tags in SAM/BAM file to reduce size and IGV load time (using pysam functions); write read length to 9th field
- reads input SAM/BAM
- simplifies input records by removing query alignment details and writing them to indexed output BAM
- removes tags to reduce size 
- can reduce BAM file size to <0.01 %
- usage: simplify_bam.py <input sam/bam> <output bam>

                                    
## unscaffold.py

Splits sequences in a FASTA at stretches of Ns of specified size and outputs a new FASTA
- Treats 'N' and 'n' the same
- does not modify gaps of other than the specified size
- optional size sort of output contigs
- reads all sequences into memory (--> therefore requires ~ 2x the memory size of the uncompressed input FASTA)


