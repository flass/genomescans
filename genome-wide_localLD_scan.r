#!/usr/bin/Rscript --vanilla

library('ape')
#library('RColorBrewer')
library('ade4')
library('getopt')

raw.args = commandArgs(trailingOnly=F)
thisscript = sub("--file=", "", grep("--file=", raw.args, value=T))
scriptdir = dirname(thisscript)
### define path to R source file of core functions
nfsource = file.path(scriptdir, 'utils-phylo.r')
if (!file.exists(nfsource)){
	# try with 'generic' path
	user = Sys.info()['user']
	homedir = file.path('/home', user)
	nfsource = file.path(homedir, 'software/genomescans/utils-phylo.r')
}
source(nfsource)

### other functions
chompnames = function(fullnames, genepat='^.+_(.+)_.+$'){
	# to take a subart of the full gene name; pattern to be modified depending how you call them in th first place and want them then
	gsub(genepat, fullnames) 
} 

LDI.aliases = c('Fisher', 'LDI', 'local_LD_index', 'Wilcox_Test_Fisher_pval')
fishermetrics = c(LDI.aliases, 'median_Fisher_pval', 'top_Fisher_pval', 'sitepairs_signif_Fisher')
all.measure.vals = c('r2', fishermetrics)

### script options

spec = matrix(c(
  'genomic.aln',     'a', 1, "character", "path to genomic alignment from which biallelic sites will be searched and LD tested",
  'out.dir',         'o', 1, "character", paste("path to an existing ouput folder; a prefix to give to oupout files can be appended, e.g.:",
                                                "'/path/to/ouput/folder/file_prefix'", sep='\n\t\t\t\t'),
  'LD.metric',       'D', 2, "character", paste("metric to report from measurement of LD, one of:",
                                                "'r2', the mean site-to-site polymorphism correlation coeff. in the window;",
                                                "'median_Fisher_pval', the median of the -log(10) p-values of Fisher's exact tests in the window;",
                                                "'top_Fisher_pval', the most significant of the -log(10) p-values of Fisher's exact tests in the window;",
                                                "'sitepairs_signif_Fisher', the number of sites in the window involved in pairs with significant",
                                                "  Fisher's exact tests (p-value scaled by the number of comparisons in the window < 0.05);",
                                                "'Local_LD_Index' (aliases: 'LDI', 'Wilcox_Test_Fisher_pval', 'Fisher') [default],",
                                                "  the -log(10) p-value of a Man-Witney-Wilcoxon U-test comparing the local distribution",
                                                "  of p-values of Fisher exact tests for pairs of neighbour biallelic sites",
                                                "  within the window vs. in the whole genome)", sep='\n\t\t\t\t'),
  'excl.ref.label',  'r', 2, "character", paste("(comma-separated) label(s) of the genome sequence(s) to exclude from the analysis;",
                                                "the first is also assumed to be the reference genome and is used to translate alignment",
                                                "coordinates into reference genome coordinates", sep='\n\t\t\t\t'),
  'window.size',     'w', 2, "integer",   "physical size (bp) of the sliding windows in which LD is evaluated [default: 3000]",
  'step',            's', 2, "integer",   "step (bp) of the sliding windows in which LD is evaluated [default: 10]",
  'nb.snp',          'm', 2, "integer",   "number of biallelic SNP within each window that are used for LD computation\n(windows with less than that are excluded from the report) [default: 20]",
  'signif.thresh',   't', 2, "double",    "threshold of Local LD Index above which the observed LD is significant [default: 5]",
  'feature.table',   'f', 2, "character", "path to GenBank feature table file indicating CDSs and matching the reference sequence coordinates",
  'chomp.gene.names','C', 2, "character", paste("regular expression pattern to shortern gene names; default to '^.+_(.+)_.+$';",
                                                "disable by providing all-matching pattern '(.+)'", sep='\n\t\t\t\t'),
  'threads',         'T', 2, "integer",   paste("number of parallel threds used for computation (beware of memory use increase)",
                                                "[default to 1: no parallel computation]", sep='\n\t\t\t\t'),
  'max.dist.ldr',    'd', 2, "integer",   paste("maximum distance between pairs of bi-allelic sites for which to compute LD",
                                                "(in mumber of intervening bi-allelic sites ; Beware: this is not uniform !!!",
                                                "polymorphism density varry across the genome !!!); by default compute the full matrix,",
                                                "which is rarely that much bigger", sep='\n\t\t\t\t'),
  'max.gap',         'g', 2, "integer",   "maximum number of allowed missing sequences to keep a site in the alignement for LD and NucDiv computations",
  'min.allele.freq', 'q', 2, "integer",   paste("minimum allele frequency (in count of sequences) in bi-alelic sites to be retained",
                                                "(minalfreq = 1 => all bi-allelic sites [default])", sep='\n\t\t\t\t'),
  'nuc.div',         'n', 2, "integer",   "window size (bp) enabless computation of nucleotidic diversity within windows of specified size",
  'crazy.plot',      'K', 2, "integer",   paste("plots the genome-wide pairwise site LD matrix to a PDF file; parameter value gives the number of sites",
                                               "to represent per page in a sub-matrix; number of cells to plot grows quickly, this can be very long",
                                                "to plot, and tedious to read as well [not done by default]", sep='\n\t\t\t\t'),
  'help',            'h', 0, "logical",   ""
), byrow=TRUE, ncol=5)
opt = getopt(spec, opt=commandArgs(trailingOnly=T))


# if help was asked for print a friendly message 
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

if ( is.null(opt$LD.metric         ) ){ opt$LD.metric         = 'Fisher'  }
if ( is.null(opt$window.size       ) ){ opt$window.size       = 3000      }
if ( is.null(opt$step              ) ){ opt$step              = 10        }
if ( is.null(opt$nb.snp            ) ){ opt$nb.snp            = 20        }
if ( is.null(opt$signif.thresh     ) ){ opt$signif.thresh     = 5         }
if ( is.null(opt$threads           ) ){ opt$threads           = 1         }
if ( is.null(opt$max.gap           ) ){ opt$max.gap           = 1         }
if ( is.null(opt$min.allele.freq   ) ){ opt$min.allele.freq   = 1         }
if ( is.null(opt$crazy.plot        ) ){ opt$crazy.plot        = -1        }

signifthresh = opt$signif.thresh
gapchars = c('-', 'N', 'n')	# various characters to consider as missing data

opt$excl.labels = strsplit(opt$excl.ref.label, split=',')[[1]]
opt$ref.label = opt$excl.labels[1]
genome.coord.str = sprintf('%s genome coordinates', opt$ref.label)

print("will use the following parameters:", quote=F)
print(opt)

if (!(opt$LD.metric %in% all.measure.vals)){
	# main reported LD metric
	stop(sprintf("Error: the value specified for the option '--LD.metric' must be among the following values: '%s'", paste(all.measure.vals, collapse="', '")))
}


## define path to result directory and optionally output file prefix

datasettag = ''

firesdir = file.info(opt$out.dir)
if (!is.na(firesdir$isdir)){
	if (firesdir$isdir & !grepl('/$', opt$out.dir)[1]){
		opt$out.dir = paste(opt$out.dir, '/', sep='')
	}
}else{
	firesdir = file.info(dirname(opt$out.dir))
	if (is.na(firesdir$isdir)){
		stop(sprintf("output folder '%s' does not exist", dirname(opt$out.dir)))
	}else{
		datasettag = basename(opt$out.dir)
	}
}
print(sprintf("result folder is: '%s'", opt$out.dir), quote=F)

## parameters
print(paste("will use", opt$threads, "cores"), quote=F)





## look for local hotspots of LD
#~ ldsearchpar = list(20, 5, 1e-10)
#~ # scans with windows with a fixed number of biallelic SNPs, variable physical size
#~ names(ldsearchpar) = c('windowsize', 'step', 'signifthresh')
# scans with windows with a fixed physical size, variable number of biallelic SNPs but sub-sampled to a maximum to get the closer to a  homogeneous statistical power along the genome; any wndow with lower number of SNP than the max has a large drop in sensitivity for high LD
# must manage a trade-off between genome coverage (achived by enlarging the windowsize and lowering the maxsize) and resolution using small and SNP-dense winndows (the inverse)
ldsearchparsub = list(opt$window.size, opt$step, opt$nb.snp, opt$signif.thresh)
names(ldsearchparsub) = c('windowsize', 'step', 'maxsize', 'signifthresh')

# load alignment file
if (!is.null(opt$genomic.aln)){
	print(sprintf("load genomic alignment from: '%s'", opt$genomic.aln), quote=F)
	prefull.aln = read.dna(opt$genomic.aln, format='fasta', as.matrix=T)
	full.aln = prefull.aln[setdiff(rownames(prefull.aln), opt$excl.labels),]
}else{
	print("no alignment file was provided; assume all the data have been pre-computed and can be loaded from RData archive files")
}
# compute correspondency between raw reference sequence and gapped alignement coordinates
nfmap2ref = paste(opt$out.dir, paste('coordinatesFullAln2Ref', 'RData', sep='.'), sep='')
if (file.exists(nfmap2ref)){ load(nfmap2ref)
}else{
	map.full2ref = fullAln2RefCoords(prefull.aln, reflabel=opt$ref.label)
	map.full2refnona = fullAln2RefCoords(prefull.aln, reflabel=opt$ref.label, bijective=F)
	save(map.full2ref, map.full2refnona, file=nfmap2ref)
}

minalfrqset = paste('minalfrq', opt$min.allele.freq, sep='')
siteset = paste('biallelicsites.max', opt$max.gap, 'gaps', sep='')

LDmetric1 = ifelse(opt$LD.metric %in% LDI.aliases, 'Fisher', opt$LD.metric)
LDmetric2 = ifelse(opt$LD.metric %in% fishermetrics, 'Fisher', opt$LD.metric)
nfoutldrad1 = paste(sprintf("LD_%s", LDmetric1), minalfrqset, siteset, ifelse(is.null(opt$max.dist.ldr), 'whole-matrix', paste('maxdist', opt$max.dist.ldr, sep='')), sep='.')
nfoutldrad2 = paste(sprintf("LD_%s", LDmetric2), minalfrqset, siteset, ifelse(is.null(opt$max.dist.ldr), 'whole-matrix', paste('maxdist', opt$max.dist.ldr, sep='')), sep='.')


## selects bi-allelic sites

# define path to file where bi-allelic site indexes are to be saved
# !!! when not satisfied by the course of the script, take care of removing saved data to mot reload them (default) next time !

nfbiali = paste(opt$out.dir, paste("BiallelicIndexes", minalfrqset, siteset, 'RData', sep='.'), sep='')
if (file.exists(nfbiali)){ 
	print(c('load \'bialraregap.i\' from', nfbiali), quote=F)
	load(nfbiali)
}else{
	print('compute \'bialraregap.i\'', quote=F)
	# biallelic sites including gaps (to be filtered out at each pairwise site comparison)
	bialgap.i = getBiAllelicIndexes(full.aln, as.logical=FALSE, minallelefreq=opt$min.allele.freq, nonsingle=FALSE, consider.gaps.as.variant=FALSE, multiproc=opt$threads)
	# biallelic sites with at most 'maxgap' missing strains
	bialraregap.i = bialgap.i[getGaplessIndexes(full.aln[,bialgap.i], as.logical=FALSE, maxgap=opt$max.gap, gapchar=c('-', 'N', 'n'), multiproc=opt$threads)]
	save(bialgap.i, bialraregap.i, file=nfbiali)
}

if ( !is.null(opt$nuc.div) ){
	## compute Nucleotidic diversity
	# parameters window scan for nucleotidic diversity; no fancy concept here
	nucdivsearchpar = list(opt$nuc.div, 1, 0.1)
	names(nucdivsearchpar) = c('windowsize', 'step', 'signifthresh')
	nfnd = paste(opt$out.dir, paste("NucDiv", 'all', 'RData', sep='.'), sep='')
	if (file.exists(nfnd)){ 
		print(c('load \'rollnucdiv\' from', nfnd), quote=F)
		load(nfnd)
	}else{
		print('compute \'rollnucdiv\'', quote=F)
		rollnucdiv = rollStats(full.aln, windowsize=opt$window.size, step=opt$step, fun=NULL, measures=c("nucdiv"), dist.model='raw', rmAbsSeq=TRUE, propgap=0.35, relpropgap=10, multiproc=opt$threads, quiet=FALSE)
		rollnucdiv$reference.position =  map.full2ref[rollnucdiv$foci]
		rollnucdiv = rollnucdiv[!is.na(rollnucdiv$reference.position),]
		rollnucdiv$scaled.nucdiv = rollnucdiv$nucdiv / max(rollnucdiv$nucdiv, na.rm=T)
		write.table(rollnucdiv, file=paste(opt$out.dir, paste("genomic_NucDiv_win100_rmAbsSeq", "allsites", 'tab', sep='.'), sep=''), row.names=F, sep='\t')
		save(rollnucdiv, file=nfnd)
		# make a nuc.div table the size of LD tables to follow
		write.table(rollnucdiv[bialraregap.i,], file=paste(opt$out.dir, paste("genomic_NucDiv_win100_rmAbsSeq", siteset, 'tab', sep='.'), sep=''), row.names=F, sep='\t')
	}
}
#ref.bial.i = map.full2ref[bialraregap.i]

nfldr = paste(opt$out.dir, paste(nfoutldrad2, 'RData', sep='.'), sep='')
if (file.exists(nfldr)){ 
	print(c('load \'bial.ldr2\' from', nfldr), quote=F)
	load(nfldr)
}else{
	nfldrtmp = paste(nfldr, 'tmp', sep='.')
	if (file.exists(nfldrtmp)){ 
		print(c('load \'lbial.ldr2\' from', nfldrtmp), quote=F)
		load(nfldrtmp)
	}else{
		## compute LD
		bial.aln = full.aln[, bialraregap.i]
		print(sprintf('compute \'lbial.ldr2\' with metric \'%s\', using %d cores', LDmetric2, opt$threads), quote=F)
		# !!! very big memory use (ex : for 42 strains * with ~10,000 biallelic SNPs: >20Gb Mem)
		# !!! use of parallism lead to big overhead, probably due to hidden duplication of the (gigantic) data matrix ; code may be optimized, but could not find the way yet. if possible, run with nbcores=1
		# or do it, just use a high memeory machine; 100G free mem to run 16 cores over a 20,000 alignment sites works OK
		lbial.ldr2 = linkageDisequilibrium(as.character(bial.aln), LDmetric2, gapchar=gapchars, multiproc=opt$threads, max.dist=opt$max.dist.ldr, discard.gaps=TRUE, mem.light=TRUE, upper.triangular=FALSE)
		gc()
		# temporary save
		save(lbial.ldr2, file=nfldrtmp)
	}
	print('expand data into a matrix', quote=F)
	bial.ldr2 = as.matrix.LDlist(lbial.ldr2)
	rm(lbial.ldr2) ; gc()
	save(bial.ldr2, file=nfldr)
}

if (opt$LD.metric %in% LDI.aliases){
	measure = "logcompldfisub" # main reported LD metric
}else{ if (opt$LD.metric == 'median_Fisher_pval'){
	measure = "medlogldfisub" # main reported LD metric
}else{ if (opt$LD.metric == 'sitepairs_signif_Fisher'){
	measure = "signifsitepairs" # main reported LD metric
}else{
	measure = "meanldrsub" # main reported LD metric
}}}
# The Fisher test p-values from a window are compared to the distribution of p-values for the genome-wide set of comparisons made in this range (i.e. within the diagonal ribbon of the square matrix of all-versus-all sites LD tests, sampled with a quantile function to get a representative sample of the same size than the window's) using a Mann-Whitney-Wilcoxon U test.
nflocld = paste(opt$out.dir, paste(sprintf("LD_%s", LDmetric1), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep='')
if (file.exists(nflocld)){ 
	load(nflocld)
	if (is.null(ldrollsub$reference.position)){ ldrollsub$reference.position = map.full2refnona[ldrollsub$foci] }
	if (is.null(rollsubsnpdens$reference.position)){ rollsubsnpdens$reference.position = map.full2refnona[rollsubsnpdens$foci] }
}else{
	print("get local LD intensity", quote=F)
	if (opt$LD.metric %in% LDI.aliases){		
		measure = "logcompldfisub" # main reported LD metric
		print(sprintf(" within windows of fixed physical size %dbp to sample bialelic sites within windows at a fixed (maximum) density of %d/window", ldsearchparsub$windowsize, ldsearchparsub$maxsize), quote=F)
		# for the null distribution, uses the total empiric distribution, taking only a sample to have comparable effective size for the test
		# sample the matrix in the diagonal ribbon that will be explored by the sliding window (i.e. short-distance comparisons, should be higher values than on the whole matrix)
		subribbonheight = ldsearchparsub$windowsize/2
		nfsubribbon = paste(opt$out.dir, paste(sprintf("LDmatrix_LowerRibbonIndexes_%dsites_physical-window%d", dim(bial.ldr2)[1], ldsearchparsub$windowsize), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep='')
		if (file.exists(nfsubribbon)){ load(nfsubribbon)
		}else{
			print("  get site range for computing empirical null distribution", quote=F)
			subribbon.i = conditionalMatrixIndexes(dim(bial.ldr2)[1], getIndexesLowerPhysicalRibbon, max.dist=subribbonheight, pos=bialraregap.i, multiproc=1)
			print(sprintf("   %d sites", length(subribbon.i)), quote=F)
			save(subribbon.i, file=nfsubribbon)
		}
		# take a sample the size of the upper triangular part of the sliding window matrix 
		subsamplesize = ldsearchparsub$maxsize*(ldsearchparsub$maxsize - 1)/2
		# use quantiles tu ensure sample is representative of the distribution
		# excludes NA located in the lower triangular matrix so as to get the right count of informative cells
		subsamplewgfi = quantile(bial.ldr2[subribbon.i], p=(0:samplesize)/subsamplesize, na.rm=T)	
		compldfisub = function(alnrange){ 
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
#~ 				print(bialrange)
				wt = wilcox.test(-log10(bial.ldr2[bialrange, bialrange]), -log10(subsamplewgfi), alternative='greater')
				return(wt$p.val)
		}}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=compldfisub, measures=measure, fun.userange=list(compldfisub=TRUE), multiproc=opt$threads)
		ldrollsub$logcompldfisub = -log10(ldrollsub$compldfisub)
		ldrollsub$compldfisub.nonas  = ldrollsub$compldfisub
		ldrollsub[is.na(ldrollsub$compldfisub), 'compldfisub.nonas'] = 1
		ldrollsub$logcompldfisub.nonas = -log10(ldrollsub$compldfisub.nonas)
	}else{ if (opt$LD.metric == 'median_Fisher_pval'){
		measure = "medlogldfisub" # main reported LD metric
		medlogldfisub = function(alnrange){ 
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
				med = median(-log10(bial.ldr2[bialrange, bialrange]), na.rm=T)
				return(med)
		}}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=medlogldfisub, measures=measure, fun.userange=list(medlogldfisub=TRUE), multiproc=opt$threads)	
	}else{ if (opt$LD.metric == 'top_Fisher_pval'){
		measure = "toplogldfisub" # main reported LD metric
		toplogldfisub = function(alnrange){ 
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
				med = max(-log10(bial.ldr2[bialrange, bialrange]), na.rm=T)
				return(med)
		}}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=toplogldfisub, measures=measure, fun.userange=list(toplogldfisub=TRUE), multiproc=opt$threads)	
	}else{ if (opt$LD.metric == 'sitepairs_signif_Fisher'){
		measure = "signifsitepairs" # main reported LD metric
		signifsitepairs = function(alnrange){ 
			if (length(alnrange)<5){ return(NA) 
			}else{
#~ 				sigwinthresh = 0.05 / (ldsearchparsub$windowsize * (ldsearchparsub$windowsize - 1) / 2) # scale p-value by the number of comparisons within the window
				sigwinthresh = 0.05 
				bialrange = which(bialraregap.i %in% alnrange)
				signif.ai = which(bial.ldr2[bialrange, bialrange] < sigwinthresh, arr.ind=T)
				nr.signif.i = unique(as.vector(signif.ai))
				return(length(nr.signif.i))
		}}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=signifsitepairs, measures=measure, fun.userange=list(signifsitepairs=TRUE), multiproc=opt$threads)	
	}else{
		measure = "meanldrsub" # main reported LD metric
		meanldrsub = function(alnrange){
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
 				return(mean(bial.ldr2[bialrange, bialrange], na.rm=TRUE))
			}
		}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=meanldrsub, measures=measure, fun.userange=list(meanldrsub=TRUE), multiproc=opt$threads)	
	}}}}
	ldrollsub$reference.position = map.full2refnona[ldrollsub$foci]
	
	
	
	print("get local biallelic SNP density", quote=F)
	
	reportsnpdens = function(alnrange){ length(alnrange) }
	rollsubsnpdens = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=reportsnpdens, measures=c("reportsnpdens"), fun.userange=list(reportsnpdens=TRUE), multiproc=opt$threads)	
	rollsubsnpdens$reference.position = map.full2refnona[rollsubsnpdens$foci]
	
	save(ldrollsub, rollsubsnpdens, file=nflocld)
	write.table(ldrollsub, paste(opt$out.dir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD-subsampled", minalfrqset, siteset, 'tab', sep='.'), sep=''))
}
print("'ldrollsub' (head/summary):", quote=F)
print(head(ldrollsub))
print(summary(ldrollsub))

hiLDfocisub = ldrollsub$reference.position[which(ldrollsub$compldfisub > ldsearchparsub$signifthresh)]
if ( !is.null(opt$nuc.div) ){
  hypervarfoci = rollnucdiv$reference.position[which(rollnucdiv$nucdiv > nucdivsearchpar$signifthresh)]
}
### optional if you have an annotation with the coordinates of the genes
### to determine which gene harbour the LD hospots

# load CDS/gene coordinates
if (!is.null(opt$feature.table)){	
	# from a GenBank feature table file
	print(c('parse \'lcds.ref.i\' from', opt$feature.table), quote=F)
	features = strsplit(readLines(opt$feature.table), split='\t')
	lcds.ref.i = list()
	# extract gene feature coordinates from file
	currcoords = NULL
	for (n in 1:length(features)){
		line = features[[n]]
		if (length(line) == 3){
			if (line[3]=='gene'){
				currcoords = sort(as.numeric(sub('[><]', '', line[1:2])))
			}
		}
		if (length(line) == 5){
			if (line[4]=='gene'){
				gene = line[5]
				if (!is.null(currcoords)){
					lcds.ref.i[[gene]] = currcoords
					currcoords = NULL
				} #else{ stop('no coordinates recorded for this gene feature') }
			}
		}
	}

	# if no feature table available:
	# must generate in a separate script (or add to this one here) a list called 'lcds.ref.i' with each element refering to a gene, and named accordingly (unique names as they are indexes)
	# each element is a vector of positions in the genomic alignment corresponding to the segment of the reference sequence where is annotated the gene, i.e. it is a subset of 'map.full2ref'

	genenames = chompnames(names(lcds.ref.i))
	# loop below to ensure that gene names are not redundant; useful if you have split genes from after a GARD breakpoint analysis
	gnames = genenames
	gcount = table(gnames)
	for (g in names(gcount)){
		k = 1 
		gi = which(gnames==g)
		gstart = sapply(genenames[gi], function(x){ as.numeric(strsplit(strsplit(x, split='_')[[1]][3], split='\\.')[[1]][1]) })
		for (ign in gi[order(gstart)]){
			gnames[ign] = paste(gnames[ign], k, sep='-')
			k = k + 1
		}
	}
	names(gnames) = names(lcds.ref.i)
	if ( !is.null(opt$nuc.div) ){
		# get locus-wise stats
		lcds.maxnucdiv = t(sapply(lcds.ref.i, function(cds.ref.i){ 
			cds.pos.i = rollnucdiv$reference.position %in% cds.ref.i
			m = max(rollnucdiv$nucdiv[cds.pos.i], na.rm=T)
			p = rollnucdiv$reference.position[cds.pos.i & rollnucdiv$nucdiv==m][1]
			return(c(p,m))
		}))
		hypervarloci = names(lcds.ref.i)[which(lcds.maxnucdiv[,2] > nucdivsearchpar$signifthresh)]
	}
	
	# report the maximum value of the target LD metric for windows covered by each CDS
	lcds.maxmeasure = t(sapply(lcds.ref.i, function(cds.ref.i){ 
		cds.pos.i = ldrollsub$reference.position %in% cds.ref.i
		m = max(ldrollsub[,measure][cds.pos.i], na.rm=T)
		p = ldrollsub$reference.position[cds.pos.i & ldrollsub[,measure]==m][1]
		return(c(p,m))
	}))
	
	hiLDlocisub = rownames(lcds.maxmeasure)[which(lcds.maxmeasure[,2] > ldsearchparsub$signifthresh)]

	hiLDgenessub =  gnames[hiLDlocisub]
	write(hiLDlocisub, file=paste(opt$out.dir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD-subsampled", minalfrqset, siteset, sprintf('hiLDloci.LDI%g', ldsearchparsub$signifthresh), sep='.'), sep=''))
	# ranked genes 
	write(unique(chompnames(hiLDlocisub)[order(lcds.maxmeasure[hiLDlocisub,2], decreasing=T)]), file=paste(opt$out.dir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD-subsampled", minalfrqset, siteset, sprintf('hiLDgenes.LDI%g', ldsearchparsub$signifthresh), sep='.'), sep=''))

}else{
	lcds.ref.i = NULL
	hiLDgenes = character(0)
}

# estimate the effect of ploymorphism distribution across site
# variation of physical size of the window for a fixed number of biallelic sites
plotphysize = function(w, plotfun, ...){
	wrange = ((w/2)+1):(length(bialraregap.i)-(w/2))
	physize = sapply(wrange, function(i){ (bialraregap.i[i+(w/2)-1] - bialraregap.i[i-(w/2)]) })
	plotfun(x=map.full2ref[bialraregap.i[wrange]], y=physize, ...)
	return(physize)
}


pdf(file=paste(opt$out.dir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, 'pdf', sep='.'), sep=''), width=21, height=13,
 title=sprintf('%s - local LD %s - windows %dbp', datasettag, opt$LD.metric, opt$physicalwindowsize))

if (opt$LD.metric %in% fishermetrics){
	sigthresh = ldsearchparsub$signifthresh
	ylabld = sprintf("LD significance in\n%dbp-wide windows [%s]", ldsearchparsub$windowsize, measure)
	histbreaks = 0:40
}else{
	sigthresh = quantile(ldrollsub$meanldrsub, p=.95, na.rm=T)
	ylabld = sprintf("LD strength (r^2) in\n%dbp-wide windows", ldsearchparsub$windowsize)
	histbreaks = (0:40)/40
}
### summary and genomic map plots
par(mar=c(8,8,8,8))
# subsampled data scan
xmax = max(map.full2ref)
if (xmax > 50000){
	ntickintervals = xmax %/% 10000
	toptick = ntickintervals * 10000
}else{ if (xmax > 20000){
	ntickintervals = xmax %/% 5000
	toptick = ntickintervals * 5000
}else{
	ntickintervals = xmax %/% 2000
	toptick = ntickintervals * 2000
}}
plot(map.full2ref[ldrollsub$foci], ldrollsub[,measure], ylab=ylabld, main=paste(datasettag, "LD scan with fixed-size windows", sep='\n'),
 xlab=genome.coord.str, col='white', xaxp=c(0, toptick, ntickintervals))
# plot areas of NA's and low-coverage
windowrect = function(pos){
	c(pos-ldsearchparsub$step/2, pos+ldsearchparsub$step/2)
}
winrects.nodata = lapply(ldrollsub$reference.position[is.na(ldrollsub[,measure])], windowrect)
pb.nodata = sapply(winrects.nodata, plotbound0, coul='grey')
winrects.fewsnp = lapply(rollsubsnpdens$reference.position[which(rollsubsnpdens$reportsnpdens < ldsearchparsub$maxsize)], windowrect)
pb.fewsnp = sapply(winrects.fewsnp, plotbound0, coul='pink')
abline(h=0:10, col=ifelse((0:10)%%5==0, 'grey', 'lightgrey'))
points(map.full2ref[ldrollsub$foci], ldrollsub[,measure], col=ifelse(ldrollsub[,measure] > sigthresh, 'red', rgb(0,0,0,0.3)))
if (!is.null(lcds.ref.i) & length(hiLDgenes)>0){
	text(labels=paste(hiLDgenessub[!is.na(hiLDgenessub)], '\'', sep='\n'), x=lcds.maxmeasure[hiLDlocisub[!is.na(hiLDgenessub)],1], y=lcds.maxmeasure[hiLDlocisub[!is.na(hiLDgenessub)],2]+.5)
}
legend('topright', fill=c('grey', 'pink'), legend=c("no data", sprintf("< %d biallelic SNP / window\n(low test power)", ldsearchparsub$maxsize)), bg='white')
plot(ldrollsub[,measure] ~ rollsubsnpdens$reportsnpdens, ylab=ylabld, xlab="Biallelic SNP density in fixed-size windows")
hist(rollsubsnpdens$reportsnpdens, breaks=0:20, xlab="Biallelic SNP density in fixed-size windows", main=paste(datasettag, "Distribution of SNP densities genome-wide", sep='\n'))
qsnpdens = quantile(rollsubsnpdens$reportsnpdens, p=c(.01, .05), na.rm=T)
abline(v=qsnpdens, col='red')
mtext(names(qsnpdens), at=qsnpdens, side=1, col='red')
hist(ldrollsub[,measure], breaks=histbreaks, xlab=ylabld, main=paste(datasettag, "Distribution of LD genome-wide", sep='\n'))
qlogcompldfisub = quantile(ldrollsub[,measure], p=c(.95, .98), na.rm=T)
abline(v=qlogcompldfisub, col='red')
mtext(names(qlogcompldfisub), at=qlogcompldfisub, side=1, col='red')

dev.off()

# test correlation of SNP density to 
if ( !is.null(opt$nuc.div) ){
	print(cor.test(rollnucdiv$nucdiv[rollnucdiv$foci %in% intersect(rollnucdiv$foci, ldrollsub$foci)], ldrollsub[ldrollsub$foci %in% intersect(rollnucdiv$foci, ldrollsub$foci), measure]))
}
print(cor.test(ldrollsub[,measure], rollsubsnpdens$reportsnpdens))


## find site-to-site significant LD (use Bonferonni correction for genome-wide multiple testing)
# !!! long and memory intensive, though less than the primary LD computations as it re-usues pre-computed data

nfsignifpairs = paste(opt$out.dir, paste(nfoutldrad2, 'significant-sitepairs.RData', sep='.'), sep='')
if (file.exists(nfsignifpairs)){ load(nfsignifpairs)
}else{
  threshpval = 0.05
  m = dim(full.aln)[1] 
  matsize =  length(which(!is.na(bial.ldr2)))
  minorallelefreqs = getMinorAlleleFreq(full.aln, multiproc=opt$threads, consider.gaps.as.variant=FALSE)
  gc()
  if (opt$LD.metric=='r2'){
	## use the Chi-squared approximation; WRONG if any of the rare alleles at biallelic sites are observed less than 5 times, e.g. with singletons
	signifr2 = qchisq(1-(threshpval/matsize), df=1)/m
	siginfpairs = which(bial.ldr2 > signifr2, arr.ind=T)
	print(sprintf("found %d pairs of biallelic sites (out of %d) in significant linkage (p < %f after Bonferonni correction for genome-wide multiple testing)", dim(siginfpairs)[1], matsize, threshpval))
	gc()
	starttime = Sys.time()
	posr2pvals = as.data.frame(t(simplify2array(mclapply(1:dim(siginfpairs)[1], function(i){
		biali = siginfpairs[i,]
		pos = bialraregap.i[biali]
		refpos = map.full2ref[pos]
		minfreqs = minorallelefreqs[pos]
		r2 = bial.ldr2[biali[1], biali[2]][[1]]
		pval = 1-pchisq(m*r2, df=1)
		qval = pval*matsize
		printProgressUpperMatrix(i, dim(siginfpairs)[1], step=1000, initclock=starttime) 
		return(c(refpos, pos, biali, r2, minfreqs, pval, qval))
	}, mc.cores=opt$threads, mc.preschedule=T))))
	colnames(posr2pvals) = c('ref.pos.1', 'ref.pos.2', 'aln.pos.1', 'aln.pos.2', 'bial.pos.1', 'bial.pos.2', 'r2', 'min.allele.freq.1', 'min.allele.freq.2', 'p.val', 'q.val')
  }else{
	## already get a p-value from using the Fisher exact test (recomended)
	siginfpairs = which(bial.ldr2*matsize < threshpval, arr.ind=T)
	print(sprintf("found %d pairs of biallelic sites (out of %d) in significant linkage (p < %f after Bonferonni correction for genome-wide multiple testing)", dim(siginfpairs)[1], matsize, threshpval))
	gc()
	starttime = Sys.time()
	posr2pvals = as.data.frame(t(simplify2array(mclapply(1:dim(siginfpairs)[1], function(i){
		biali = siginfpairs[i,]
		pos = bialraregap.i[biali]
		minfreqs = minorallelefreqs[pos]
		refpos = map.full2ref[pos]
		pval = bial.ldr2[biali[1], biali[2]][[1]]
		qval = pval*matsize
		printProgressUpperMatrix(i, dim(siginfpairs)[1], step=1000, initclock=starttime) 
		return(c(refpos, pos, biali, minfreqs, pval, qval))
	}, mc.cores=opt$threads, mc.preschedule=T))))
	colnames(posr2pvals) = c('ref.pos.1', 'ref.pos.2', 'aln.pos.1', 'aln.pos.2', 'bial.pos.1', 'bial.pos.2', 'min.allele.freq.1', 'min.allele.freq.2', 'p.val', 'q.val')
  }
  gc()
  posr2pvals$site.dist = abs(posr2pvals$ref.pos.2 - posr2pvals$ref.pos.1)
  write.table(posr2pvals, file=gsub('\\.RData', '.tab', nfsignifpairs))
  # for all comparisons
  if (opt$LD.metric=='r2'){
    allqvals = log10((1-pchisq(m*bial.ldr2, df=1))*matsize)
  }else{
    allqvals = log10(as.numeric(bial.ldr2)*matsize)
  }
  gc()
  rangeqvals = seq(from=floor(min(allqvals, na.rm=T)), to=floor(max(allqvals, na.rm=T)), by=0.5)
  gc()
  freqqvals = sapply(rangeqvals, function(k){ length(which(allqvals>=k & allqvals<k+0.5)) })
  gc()
  save(matsize, posr2pvals, siginfpairs, allqvals, rangeqvals, freqqvals, file=nfsignifpairs)
}
if (dim(posr2pvals)[1] > 0){
	# only for the significant comparisons
	intersitedistclases = seq(from=0, to=floor(max(posr2pvals$site.dist/1000, na.rm=T)), by=1)
	qvalbyintersitedistclases = lapply(intersitedistclases, function(k){ log10(posr2pvals$q.val[posr2pvals$site.dist >= k*1000 & posr2pvals$site.dist < (k+1)*1000]) })
	names(qvalbyintersitedistclases) = sapply(intersitedistclases, function(k){ paste(k, k+1, sep='-') })
	#~ sitedistbyqvalclass = sapply(seq(from=floor(min(allqvals, na.rm=T)), to=floor(max(allqvals, na.rm=T)), by=0.5), function(k){ length(which(allqvals>=k & allqvals<k+0.5)) })

	pdf(file=paste(opt$out.dir, nfoutldrad2, 'significant-pairs.pdf', sep='.'), sep=''), height=12, width=12)
	barplot(log10(freqqvals), names.arg=paste('[', rangeqvals, '; ', rangeqvals+0.5, ']', sep=''), ylim=c(0,10), las=2, xlab='log10(corrected p-values)', ylab='log10(frequency)',  cex.lab=1.5,main='Distribution of p-values of Fisher\'s exact test for LD significance')
	hist(posr2pvals$site.dist, breaks=(0:24)*10000, cex.lab=1.5, xlab='inter-site distance')
	plot(x=((1:100)/100), y=quantile(posr2pvals$site.dist, p=(1:100)/100), xlab='quantiles', ylab='inter-site distance',  cex.lab=1.5,main=sprintf('Cumulative distribution function of inter-site distances\nfor significantly linked site (p < %g after Bonferonni cor.)', threshpval))
	boxplot(qvalbyintersitedistclases, xlab='inter-site distance (kb)', ylab='log10(corrected p-values)', cex.lab=1.5, main='Distribution of significant p-values with site distance')
	#~ hist(log10(posr2pvals$q.val), breaks=101)
	#~ plot(log10(posr2pvals$q.val) ~ posr2pvals$site.dist)
	dev.off()
}else{
	print("found nosugnificant site associations", quote=F)
}

if (crazy.plot > 0){
	## plots the LD matrix (the one computed with potentially included rare gaps) to PDF
	# USE ONLY ON SMALL MATRICES !!!!!!
	area.mat.i = map.full2ref[bialraregap.i]
	pdf(paste(opt$out.dir, paste("genome-wide_LD", minalfrqset, siteset, 'pdf', sep='.'), sep=''), width=20, height=10)
	plotLDHeatMaps(simplify2array(bial.ldr2), area.mat.i, nlabs=opt$crazy.plot, main=paste(datasettag, "correlation coefficient r^2 between sites", sep='\n'))
	dev.off()
}

