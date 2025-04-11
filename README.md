# utilities
Collection of (hopefully) useful things

# descriptions

## ```slurmer``` & ```re-slurm```

Wrapper to submit SLURM/SBATCH and organize the submission and output files
- usage: ```slurmer <prompt> <jobname> <runtime> <memory> <num_cores_per_task> <num_tasks> <num_nodes> <partition>```
- set defaults (such as run time, memory, partition, etc) and user email address in the SETUP section of the script before usage
- creates directory .slurm with files .slurm/slurm.name.submit and .slurm/slurm.name.out
- .slurm/slurm.name.submit contains submission time and command
- .slurm/slurm.name.out contains all STDERR and STDOUT
- example: slurm "echo test" test
- make executable (chmod +x slurmer) and add do $PATH to call it from anywhere
- ```re-slurm <.slurm/slurm.*.submit>``` re-launches an existing submission file from a slurmer job


## ```screener``` & ```re-screen```

Submit jobs e.g. in a loop on a server without job scheduler
- usage: ```screener <"command"> <job name>```
- creates directory .screen with files .screen/screen.name.submit and .screen/screen.name.out
- .screen/screen.name.submit contains submission time and command
- .screen/screen.name.out contains all STDERR and STDOUT
- example: screener "echo test" test
- make executable (chmod +x screener) and add do $PATH to call it from anywhere
- ```re-screen <.screen/screen.*.submit>``` re-launches an existing submission file from a screener job

  
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


## ```mito_contamination_scan.py```

Script to extract allelic depth from VCF at mitochondrial variant sites to check for sample contamination.
- rationale: Mitochondria are haploid, i.e. there should be no variation except random sequencing error.
- exception: [heteroplasmy](https://en.wikipedia.org/wiki/Heteroplasmy) (rare but possible)
- usage: ```python extract_allelic_depths.py <vcf_path> <primary_ids_path> <mito_name> <mito_len> <base_error_rate> <ad_output_path> <stats_output_path>```


## ```plot_bam_wga.py```

Plot whole genome alignment from a BAM file.
- visualize pairwise alignments of multiple chromosomes or contigs from a BAM file
- collinear alignments in black, inversions highlighted in red, unaligned regions in grey.
- usage: ```python plot_bam_alignments.py <input_path> <output_path> <faidx_path> <assoc_file_path> <regions_file_path>```
- for details check ```plot_bam_alignments.py --help```


## ```subset_beagle.py```

Subset a BEAGLE (or BEAGLE.gz)
- subset to a specified set of samples (-s) and/or every nth line (-l)
- usage: ```python subset_beagle.py <input_beagle_path> <output_beagle_path>```
- for details and options check ```subset_beagle.py --help```


## ```vcf_pca.py```

Perform PCA on genotypes (GT) in a VCF file using scikit-allel and generate interactive/annotated plotly plots.
- dependencies: scikit-allel and plotly
- run ```vcf_pca.py pca``` first to perform PCA and generate output text files, then optionally ```vcf_pca.py plot``` for plotting
- ```vcf_pca.py pca``` reads variants from a VCF file, mean-imputes potential missing variant calls and performs PCA using scikit-allel. It can be run either for a specific region of a chromosome (```-r chr:start-stop```) or for all variants contained in the VCF if ```-r``` is not specified. A subset of samples in the VCF can be specified as a comma-separated list supplied to ```-s``` and the minor allele frequency (MAF) threshold set with ```-m``` (default: 0.01). The ```pca```module generates two output files, one containing the eigenvectors, the other containing the eigenvalues.
- ```vcf_pca.py plot``` reads these output files and creates an interactive HTML plot (and/or regular PDF, see ```-f```). ```-p```can be used to specify the principal components to plot, ```-m``` to supply a metadata file which will be used to annotate the HTML plot (the first column must contain the same sample IDs as used in the VCF file). See ```vcf_pca.py plot -h``` for more option, including reflecting data along PC axes, setting plot colors for groups of samples or using a continuous color scheme based on a numeric metadata column.



