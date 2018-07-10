# genomescans
R/Python modules and scripts for the analysis of aligned genomic sequences, to analyse of phylotype/haplotype structure of a genomic dataset and to perform genome-wide scans for the detection of hotspots of diversity, LD, etc.

This software suite was used in the following paper: Lassalle, F. et al. (2016) Islands of linkage in an ocean of pervasive recombination reveals two-speed evolution of human cytomegalovirus genomes. *Virus Evolution* 2 (1), vew017.[doi:10.1093/ve/vew017](http://dx.doi.org/10.1093/ve/vew017), where mehtods are described in more detail. Please cite this paper if using any of the software below.

## bayesbipartprofile suite
The bayesbipartprofile suite intends to explore the local phylogenetic structure within genomes of recombining species. It reconstructs haplotypes spanning genome regions, looking for any conserved phylogenetetic relationships between variable sets of strains/species/isolates across loci. From a dataset of bayesian samples of gene trees, the [bayesbipartprofile.py] script generates a database of bipartiations and search for similarities between them. (parser last tested and working on ouput from MrBayes 3.2.2)

Then, the [bayesbipartprofile.r] script builds matrices of bipartition support (two metrics are reported: posterior probability and compatibility score) across loci to detect conserved tracks of clonal phylogenetic structure, i.e. haplotypes, and provide text table and graphic output.
This includes a map of phylogenetic haplotypes, and heatmaps of split profile correlation matrices, that provide a direct insight into *long-range linkage disequilibrium (LD)*.

Example:
```sh
# generate database of splits
python bayesbipartprofile.py --ltax /path/to/orderred_taxon_list --bs.thresh.ref.bip=0.35 /path/to/orderred_gene_list /path/to/mrbayes_result_directory /path/to/bayesbipartprofile_output_directory
# make correlation matrices and plots 
bayesbipartprofile.r --gene.list.is.ordered --new /path/to/bayesbipartprofile_output_directory
```
Below is an example of the matrix of bipartition compatibility score correlation between genes in a dataset of 142 HCMV genomes (data from [Lassalle et al. (2016)]):

![HCMV_bipart_compat_cor]


## [genome-wide_localLD_scan.r](https://github.com/flass/genomescans/blob/master/genome-wide_localLD_scan.r)
The [genome-wide_localLD_scan.r] script performs a genome-wide search for *localized LD* based on the distribution of alleles at biallelic polymorphic sites (SNPs) in a multiple sequence alignment.
It mainly rely on `linkageDisequilibrium()` and `rollStats()` functions in [utils-phylo.r] module.

Many options are available, as described below (result of call with `--help` option):

```
Usage: ./genome-wide_localLD_scan.r [-[-genomic.aln|a] <character>] [-[-out.dir|o] <character>] [-[-LD.metric|D] [<character>]] [-[-ref.label|r] [<character>]] [-[-window.size|w] [<integer>]] [-[-step|s] [<integer>]] [-[-nb.snp|m] [<integer>]] [-[-signif.thresh|t] [<double>]] [-[-feature.table|f] [<character>]] [-[-chomp.gene.names|C] [<character>]] [-[-threads|T] [<integer>]] [-[-max.dist.ldr|d] [<integer>]] [-[-max.gap|g] [<integer>]] [-[-min.allele.freq|q] [<integer>]] [-[-crazy.plot|K] [<integer>]] [-[-help|h]]
    -a|--genomic.aln         path to genomic alignment from which biallelic sites will be searched and LD tested
    -o|--out.dir             path to an existing ouput folder; a prefix to give to oupout files can be appended, e.g.: 
                               '/path/to/ouput/folder/file_prefix'
    -D|--LD.metric           metric to report from measurement of LD, one of: 'r2' (correlation coeff.) 
                               or 'Fisher' (Local LD Index [default]: -log(10) p.value of a Man-Witney-Wilcoxon U-test 
                               comparing the local distribution of p-values of Fisher exact tests for pairs of neighbour 
                               biallelic sites within the window vs. in the whole genome)
    -r|--ref.label           (comma-separated) label(s) of the genome sequence(s) to exclude from the analysis; 
                               the first is also assumed to be the reference genome and is used to translate alignment 
                               coordinates into reference genome coordinates
    -w|--window.size         physical size (bp) of the sliding windows in which LD is evaluated [default: 3000]
    -s|--step                step (bp) of the sliding windows in which LD is evaluated [default: 10]
    -m|--nb.snp              number of biallelic SNP within each window that are used for LD computation
                               (windows with less than that are excluded from the report) [default: 20]
    -t|--signif.thresh       threshold of Local LD Index above which the observed LD is significant [default: 5]
    -f|--feature.table       path to GenBank feature table file indicating CDSs and matching the reference sequence coordinates
    -C|--chomp.gene.names    regular expression pattern to shortern gene names; default to '^.+_(.+)_.+$';
                               disable by providing all-matching pattern '(.+)'
    -T|--threads             number of parallel threds used for computation (beware of memory use increase)
                               [default to 1: no parallel computation]
    -d|--max.dist.ldr        maximum distance between pairs of bi-allelic sites for which to compute LD
                               (in mumber of intervening bi-allelic sites ; Beware: this is not uniform !!! 
                               polymorphism density varry across the genome !!!); by default compute the full matrix, 
                               which is rarely that much bigger.
    -g|--max.gap             maximum number of allowed missing sequences to keep a site in the alignement for LD and NucDiv computations
    -q|--min.allele.freq     minimum allele frequency (in count of sequences) in bi-alelic sites to be retained
                               (minalfreq = 1 => all bi-allelic sites [default])
    -K|--crazy.plot          plots the genome-wide pairwise site LD matrix to a PDF file; parameter value gives the number of sites
                               to represent per page in a sub-matrix; number of cells to plot grows quickly, this can be very long
                               to plot, and tedious to read as well [not done by default]
    -h|--help
```

## detect_recomb scripts

These scripts provide wrappers for recombination detection programs, including GeneConv, PHI and HyPhy's SBP/GARD in particular.

Both script allow the execution of many unique jobs. For that you've got to provide a list of tasks, i.e. a file in which each line is a path to a sequence alignment you desire to scan = one task. The Python script run tasks sequentially, and the qsub script (shell wrapper to submit parallel jobs of the Python script to a SGE type of computer cluster) can schedule them for execution in parallel, each job dealing with a chunk of tasks to be executed sequentially).

A requirement to use GARD is to specify the location of the `mpirun` command and of the `hyphy/` root folder (set `mpipath` and `hyphypath` variables directly in the Python script, or set the `$mympi` and `$myhyphy` variables in the qsub script).

[Lassalle et al. (2016)] http://dx.doi.org/10.1093/ve/vew017
[bayesbipartprofile.py] https://github.com/flass/genomescans/blob/master/bayesbipartprofile.py
[bayesbipartprofile.r] https://github.com/flass/genomescans/blob/master/bayesbipartprofile.r
[genome-wide_localLD_scan.r] https://github.com/flass/genomescans/blob/master/genome-wide_localLD_scan.r
[utils-phylo.r] https://github.com/flass/genomescans/blob/master/utils-phylo.r
[HCMV_bipart_compat_cor]: https://github.com/flass/genomscan/blob/master/figures/HCMV_longLD_bipart_compat_score_correlation.png
