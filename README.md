# utilities
Collection of (hopefully) useful things

# description

## unscaffold.py

Splits sequences in a FASTA at stretches of Ns of specified size and outputs a new FASTA
- Treats 'N' and 'n' the same
- does not modify gaps of other than the specified size
- optional size sort of output contigs
- reads all sequences into memory (--> therefore requires ~ 2x the memory size of the uncompressed input FASTA)

## screener

Submit jobs e.g. in a loop on a server without job scheduler.
- usage: screener "command" name
- creates directory .screen with files .screen/screen.name.submit and .screen/screen.name.out
- .screen/screen.name.submit contains submission time and command
- .screen/screen.name.out contains all STDERR and STDOUT
- example: screener "echo test" test
