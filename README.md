# utilities
Collection of (hopefully) useful things

# descriptions

## ```count_alleles.py```

Simple script to count missing (.), A, C, G, T alleles per site from a VCF.
- reads VCF (with or without header) from STDIN and prints to STDOUT
- input: SNPs (invariant, bialellic, multiallelic)
- only for diploid samples, phased (|) or unphased (/)
- output: TSV with CHROM, POS . A C G T counts
- set '-h' flag to output a header line


## ```type_winpca_locus.py```

Script to genotype PC values in a single WinPCA output window using two thresholds.
- usage: ```python type_winpca_locus.py <pc_tsv> <window_pos> <upper_threshold> <lower_threshold>```
- reads WinPCA PC output file (```.pc_*.tsv.gz```) (```<pc_tsv>```)
- takes a single window position (```<window_pos>```) plus  ```<upper_threshold>``` and ```<lower_threshold>```
- writes ```sample_id  pc_value  typing``` to ```STDOUT``` where pc_value larger than ```<upper_threshold>``` is typed as ```2```, lower than ```<lower_threshold>``` typed as ```0```, and anything in between typed as ```1``` (```1``` typing is inclusive of both thresholds)
- can be used to type inversions or other divergent haplotypes


## ```mito_contamination_scan.py```

Script to extract allelic depth from VCF at mitochondrial variant sites to check for sample contamination.
- rationale: Mitochondria are haploid, i.e. there should be no variation except random sequencing error.
- exception: [heteroplasmy](https://en.wikipedia.org/wiki/Heteroplasmy) (rare but possible)
- usage: ```python extract_allelic_depths.py <vcf_path> <primary_ids_path> <mito_name> <mito_len> <base_error_rate> <ad_output_path> <stats_output_path>```


## ```plot_assoc.py```
- create interactive or static plots of association data, (e.g. GWAS, r2)
- uses FAIDX index (or any other TSV with chrosome name and size as the first two columns) to infer chromosome boundaries
- ```-c``` to specify column contains chromosome name (default: chrom)
- ```-v``` to specify column contains variant position on the chromosome (default: pos)
- ```-p``` to specify column containing data to plot (default: p)
- ```-l``` to control whether to compute and plot -log10 of the data
- ```-a``` to control alpha value to infer and draw significance threshold (default: 0.01)
- ```-r``` to remove n-th quantile of weakes associations from plot, mainly to reduce file size (removing 90% of data by specifying 0.9 will often produce a nearly identically looking plot with lots of overlapping low associations not shown)
- ```-x``` and ```-y``` to set figure width (default: 1500 p) and height (default: 500 p)
- ```-f``` to select one or a comma-separated list of output formats (default: HTML)


## ```plot_bam_wga.py```

Plot whole genome alignment from BAM.
- if you can, I'd recommend to use ```plot_paf_wga.py``` instead
- visualize pairwise alignments of multiple chromosomes or contigs from a BAM file
- collinear alignments in black, inversions highlighted in red, unaligned regions in grey.
- usage: ```python plot_bam_alignments.py <input_path> <output_path> <faidx_path> <assoc_file_path> <regions_file_path>```
- for details check ```plot_bam_alignments.py --help```


## ```plot_paf_wga.py```

Plot whole genome alignment(s) from PAF.
-  successor for ```plot_bam_wga.py```
- useful to examine assembly completeness by aligning a primary and secondary assembly to each other, or both to a shared reference
-  input pairwise genome alignments in ```.PAF``` format may be created with [minimap](https://github.com/lh3/minimap2) or the even faster [FastGA](https://github.com/thegenemyers/FASTGA)
-  visualizes all reference genome sequences (chromosomes + contigs) and likewise the aligned assembly's sequences on top and connects pairwise alignment blocks
-  optionally, a second aligned assembly can be visualized by specifying a second ```.PAF``` file to ```-s```
-  sequences of the aligned assemblies are automatically ordered using the alignments to the reference genome, the sequence order of the reference genome can be set using ```-p```, which will be shown as a third panel below
- by default, all sequences contained in the PAF files are shown, which means that if there are sequences without any alignments they might be missing
- therefore, ```.FAIDX``` files may be soecified for the ref/primary/secondary assemblies using ```-a```/```-b```/```-c``` to provide all sequences (this adds unaligned sequences to the plot)
- gaps may be optionally plotted for the ref/primary/secondary assemblies if gap locations are provided (```-x```/```-y```/```-z```); they can be easily generated with [seqtk](https://github.com/lh3/seqtk) (```seqtk gap <FASTA>```)
- custom labels for the assemblies may be added with ```-n```
- collinear alignments are shown in grey (alternating light/dark, according to reference chromosomes), inversions in blue, gaps in red (colors can be changed in the CONFIG block)
- by default, only aligments of MQ>=60 are visualized (this may be adjusted in the CONFIG block by changing ```MIN_ALN_SCORE```)
- supported output formats (```-f```) are ```HTML``` (interactive), ```PDF```, ```SVG``` and ```PNG``` (more than one can be specified as a comma-separated list)
- in HTML output, individual chromosomes (and gaps if plotted) can be toggled on and off, and one can zoom in and out and hover displays provide additional information
- more details are in the help message ```plot_paf_wga.py --help```


## ```separate_colors.py```

Plot different color ranges (red, yellow, blue) and shadows as separate panels below the original image.
- usage: ```python separate_colors.py <input_file> <output_file>```
- for details check ```separate_colors.py --help```
- color ranges and transformations can be adjusted in the CONFIG block of the .py file
- dependency: opencv ([https://github.com/pysam-developers/pysam](https://opencv.org))


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


## ```subset_beagle.py```

Subset a BEAGLE (or BEAGLE.gz)
- subset to a specified set of samples (-s) and/or every nth line (-l)
- usage: ```python subset_beagle.py <input_beagle_path> <output_beagle_path>```
- for details and options check ```subset_beagle.py --help```


## ```unscaffold.py```

Splits sequences in a FASTA at stretches of Ns of specified size and outputs a new FASTA
- Treats 'N' and 'n' the same
- does not modify gaps of other than the specified size
- optional size sort of output contigs
- reads all sequences into memory (--> therefore requires ~ 2x the memory size of the uncompressed input FASTA)
- usage: ```unscaffold.py <input path> <output path> <contig-prefix> <gap size> <out fasta line length> <size sort True/False>```
- for details check ```unscaffold.py --help```


## ```vcf_pca.py```

Perform PCA on genotypes (GT) in a VCF file using scikit-allel and generate interactive/annotated plotly plots.
- dependencies: scikit-allel and plotly
- run ```vcf_pca.py pca``` first to perform PCA and generate output text files, then optionally ```vcf_pca.py plot``` for plotting
- ```vcf_pca.py pca``` reads variants from a VCF file, mean-imputes potential missing variant calls and performs PCA using scikit-allel. It can be run either for a specific region of a chromosome (```-r chr:start-stop```) or for all variants contained in the VCF if ```-r``` is not specified. A subset of samples in the VCF can be specified using ```-s``` and the minor allele frequency (MAF) threshold set with ```-m``` (default: 0.01). The ```pca```module generates two output files, one containing the eigenvectors, the other containing the eigenvalues.
- ```vcf_pca.py plot``` reads these output files and creates an interactive HTML plot (and/or regular PDF, see ```-f```). ```-p```can be used to specify the principal components to plot, ```-m``` to supply a metadata file which will be used to annotate the HTML plot (the first column must contain the same sample IDs as used in the VCF file). See ```vcf_pca.py plot -h``` for more option, including reflecting data along PC axes, setting plot colors for groups of samples or using a continuous color scheme based on a numeric metadata column.


## ```trim_circ_seq.py```

Trim trailing sequence from a linearized circular DNA assembly by identifying and removing the region that matches the start of the sequence.
Tries to find a single secondary unique match of the leading k base pairs, assuming that this is where the assembly would begin to overlap if circularized. Starting from k=50, k is reduced by one with each unsuccessful iteration. If a unique secondary match if found it is trimmed from the output as well as any trailing sequence. If more than one unique secondary match is found before a k yields exactly one secondary unique match, there might be a mismatch between in the flanking/overlapping assembly portions. Therefore, the search string is shifted to k+1 and another run is performed. If a unique secondary match is now found for the shifted search string, this region plus the leading offset (which includes the mismatch) and any trailing region is trimmed off. If a mismatch is found, this is reported in the output. If no second unique match can be found for any k >= 5 (allowing for one mismatch) the program exits with an error message.
- usage: ```python trim_circ_seq.py <input_path> <output_path>```
- input/output FASTA may be gzipped
- the search range for k is set to 50..5 and can be updated in the script (seed_len_max and seed_len_min)
