# genomescans
R/Python modules and scripts for the analysis of aligned genomic sequences, to analyse of phylotype/haplotype structure of a genomic dataset and to perform genome-wide scans for the detection of hotspots of diversity, LD, etc.

This software suite was used in the following paper: [F. Lassalle et al. (2016) Islands of linkage in an ocean of pervasive recombination reveals two-speed evolution of human cytomegalovirus genomes. Virus Evolution 2 (1), vew017.](http://dx.doi.org/10.1093/ve/vew017), where mehtods are described in detail. Please cite the paper if using any of the software below.

## bayesbipartprofile suite
The bayesbipartprofile suite intends to explore the local phylogenetic structure within genomes of recombining species. It reconstructs haplotypes spanning genome regions, looking for any conserved phylogenetetic relationships between variable sets of strains/species/isolates across loci. From a dataset of bayesian samples of gene trees, the '[bayesbipartprofile.py](https://github.com/flass/genomescans/blob/master/bayesbipartprofile.py)' script generates a database of bipartiations and search for similarities between them. Then, the [bayesbipartprofile.r](https://github.com/flass/genomescans/blob/master/bayesbipartprofile.r)' script builds matrices of bipartition support (using posterioir probability and compatibility metrics) across loci to detect conserved tracks of clonal phylogenetic stracture, i.e. haplotypes, and provide text table and graphic output.
(parser last tested and working on ouput from MrBayes 3.2.2)

## [genome-wide_LD_template.r](https://github.com/flass/genomescans/blob/master/genome-wide_LD_template.r)
The '[genome-wide_LD_template.r](https://github.com/flass/genomescans/blob/master/genome-wide_LD_template.r)' script performs a genome-wide search for linkage disequilibrium (LD) based on the distribution of alleles at biallelic polymorphic sites (SNPs) in a multiple sequence alignment.
It mainly rely on `linkageDisequilibrium()` and `rollStats()` functions in '[utils-phylo.r](https://github.com/flass/genomescans/blob/master/utils-phylo.r)' module.

Many options are available, and have to be specified by changing the values of several variables directly in the script;this includes the paths to the relevant input and output files. Here is a sample of them:

### environment variables
```R
nfsource = 'path/to/utils-phylo.r'
nbcores = 1		# number of cores to be used in parallel ; for LD computation on big datasets (>100kb are big), prefer use nbcores=1 (see below), unless large memory is available
resultdir = 'path/to/data/and/results/directory'
full.aln = read.dna('path/to/genome_alignment.fasta', format='fasta')
reflabel = 'reflabel'	# label in the alignement of the strain to use for reference genome coordinates
```
### parameters of genome-wide LD computation
```R
max.dist.ldr = 3000	# maximum distance between pairs of bi-allelic sites for LD computation (in mumber of intervening bi-allelic sites ; not uniform !!! polymorphism varry in density across the genome !!!) 
# the last comment is good reason not to do it, so alternatively:
max.dist.ldr = NULL	# compute the full matrix; not that much bigger depending on the datset
maxgap = 1			# maximum number of allowed missing sequences to keep a site in the alignement for LD and NucDiv computations
minalfrq = 1		# minimum allele frequency (in count of sequences) in bi-alelic sites to be retained (minalfreq = 1 => all bi-allelic sites)
gapchars = c('-', 'N', 'n')	# various characters to consider as missing data
LDmetric = 'r2' # if one want to report the r-squared statistic of LD; the Chi-squared approximation can be used a posteriori (Chi-squared only depend on r2 value) when the minor allele counts are not too low, or alternatively:
LDmetric = 'Fisher' # if one want to report the p-value of a Fisher exact test for significance of the LD; recommended as the situation above is rarely met in most of the microbial pathogen genomes
```

### parameters for window scans looking for local hotspots of LD or nucleotidic diversity
```R
# parameters for scans with windows with a fixed number of biallelic SNPs, variable physical size
ldsearchpar = list(20, 5, 1e-10)
names(ldsearchpar) = c('windowsize', 'step', 'signifthresh')
# parameters for scans with windows with a fixed physical size, variable number of biallelic SNPs but sub-sampled to a maximum to get the closer to a  homogeneous statistical power along the genome; any wndow with lower number of SNP than the max has a large drop in sensitivity for high LD
# must manage a trade-off between genome coverage (achived by enlarging the windowsize and lowering the maxsize) and resolution using small and SNP-dense winndows (the inverse)
ldsearchparsub = list(700, 10, 20, 1e-5)
names(ldsearchparsub) = c('windowsize', 'step', 'maxsize', 'signifthresh')
# parameters for window scan for nucleotidic diversity; no fancy concept here
nucdivsearchpar = list(100, 1, 0.1)
names(nucdivsearchpar) = c('windowsize', 'step', 'signifthresh')
```
### plotting and reporting options
```R
nfmapcds = '/path/to/map/genes/location/to/ref/genome.RData'
# must be generated in a separate script (or add code to this one) a list called 'lcds.ref.i' with each element refering to a gene, and named accordingly (choose unique names as they are indexes)
# each element is a vector of positions in the genomic alignment corresponding to the segment of the reference sequence where is annotated the gene, i.e. it is a subset of global map object 'map.full2ref'
```
