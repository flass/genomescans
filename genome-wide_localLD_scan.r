#!/usr/bin/Rscript --vanilla

library('ape')
library('RColorBrewer')
library('gplots')
library('ade4')
library('getopt')

### define path to R source file of core functions
user = Sys.info()['user']
homedir = paste('/home', user, sep='/')
nfsource = paste(homedir, 'scripts/misc/utils-phylo.r', sep='/')
#~ nfsource = 'path/to/utils-phylo.r'
source(nfsource)

### other functions
chompnames = function(fullnames, genepat='^.+_(.+)_.+$'){
	# to take a subart of the full gene name; pattern to be modified depending how you call them in th first place and want them then
	gsub(genepat, fullnames) 
} 

### script options

spec = matrix(c(
  'genomic.aln',   'a', 1, "character", "path to genomic alignment from which biallelic sites will be searched and LD tested",
  'out.dir',       'o', 1, "character", "path to ouput folder; can be appended with the prefix to give to oupout files, e.g.: '/path/to/ouput/folder/file_prefix'",
  'LD.metric',     'D', 2, "character", "metric to report from measurement of LD, one of: 'r2' (correlation coeff.), 'Fisher' (Local LD Index [default]: -log(10) p.value of a Man-Witney-Wilcoxon U-test comparing the local distribution of p-values of Fisher exact tests for pairs of neighbour biallelic sites within the window vs. in the whole genome)",
  'excl.ref.label','r', 2, "character", "(comma-separated) name(s) of the (reference) genome to exclude from the analysis",
  'window.size',   'w', 2, "integer",   "physical size (bp) of the sliding windows in which LD is evaluated [default: 3000]",
  'step',          's', 2, "integer",   "step (bp) of the sliding windows in which LD is evaluated [default: 10]",
  'nb.snp',        'm', 2, "integer",   "number of biallelic SNP within each window that are used for LD computation\n(windows with less than that are excluded from the report) [default: 20]",
  'signif.thresh', 't', 2, "double",    "threshold of Local LD Index above which the observed LD is significant [default: 5]",
  'feature.table', 'f', 2, "character", "path of ouput PDF file for plots",
  'chomp.gene.names', 'C', 2, "character", "regular expression pattern to shortern gene names; default to '^.+_(.+)_.+$'; disable by providing all-matching pattern '(.+)'",
  'threads',       'T', 2, "integer",    "number of parallel threds used for computation (beware of memory use increase) [default to 1: no parallel compuation]",
  'max.dist.ldr',  'd', 2, "integer",    "maximum distance between pairs of bi-allelic sites for which to compute LD (in mumber of intervening bi-allelic sites ; Beware: this is not uniform !!! polymorphism density varry across the genome !!!); by default compute the full matrix, which is rarely that much bigger",
  'max.gap',       'g', 2, "integer",    "maximum number of allowed missing sequences to keep a site in the alignement for LD and NucDiv computations",
  'min.allele.freq','q', 2, "integer",    "minimum allele frequency (in count of sequences) in bi-alelic sites to be retained (minalfreq = 1 => all bi-allelic sites [default])",
  'crazy.plot',    'K', 2, "integer",    "plots the genome-wide pairwise site LD matrix to a PDF file; parameter value gives the number of sites to represent per page in a sub-matrix; number of cells to plot grows quickly, this can be very tedious to plot, and to read as well [not doone by default]",
  'help',          'h', 0, "logical",    ""
), byrow=TRUE, ncol=5);
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
if ( is.null(opt$max.dist.ldr      ) ){ opt$max.dist.ldr      = 3000      }
if ( is.null(opt$max.gap           ) ){ opt$max.gap           = 1         }
if ( is.null(opt$min.allele.freq   ) ){ opt$min.allele.freq   = 1         }
if ( is.null(opt$crazy.plot        ) ){ opt$crazy.plot        = -1        }



signifthresh = 10^(-opt$signif.thresh)
gapchars = c('-', 'N', 'n')	# various characters to consider as missing data

print("will use the following parameters:", quote=F)
print(opt)

exclreflabels = strsplit(opt$excl.ref.label, split',')[[1]]

## define path to result directory

datasettag = ''

firesdir = file.info(opt$resultdir)
if (!is.na(firesdir$isdir)){
	if (firesdir$isdir & !grepl('/$', opt$resultdir)[1]){
		opt$resultdir = paste(opt$resultdir, '/', sep='')
	}
}else{
	firesdir = file.info(dirname(opt$resultdir))
	if (is.na(firesdir$isdir)){
		stop(sprintf("output folder '%s' does not exist", dirname(opt$resultdir)))
	}else{
		datasettag = basename(opt$resultdir)
	}
}
print(sprintf("result folder is: '%s'", opt$resultdir), quote=F)

nfoutldrad = paste(sprintf("LD_%s", opt$LD.metric), minalfrqset, siteset, ifelse(is.null(opt$max.dist.ldr), 'whole-matrix', paste('maxdist', opt$max.dist.ldr, sep=''))

## parameters
print(paste("will use", opt$threads, "cores"), quote=F)





## look for local hotspots of LD
#~ ldsearchpar = list(20, 5, 1e-10)
#~ # scans with windows with a fixed number of biallelic SNPs, variable physical size
#~ names(ldsearchpar) = c('windowsize', 'step', 'signifthresh')
# scans with windows with a fixed physical size, variable number of biallelic SNPs but sub-sampled to a maximum to get the closer to a  homogeneous statistical power along the genome; any wndow with lower number of SNP than the max has a large drop in sensitivity for high LD
# must manage a trade-off between genome coverage (achived by enlarging the windowsize and lowering the maxsize) and resolution using small and SNP-dense winndows (the inverse)
ldsearchparsub = list(opt$physicalwindowsize, opt$step, opt$nb.snp, signifthresh)
names(ldsearchparsub) = c('windowsize', 'step', 'maxsize', 'signifthresh')
# paprmaeters window scan for nucleotidic diversity; no fancy concept here
nucdivsearchpar = list(100, 1, 0.1)
names(nucdivsearchpar) = c('windowsize', 'step', 'signifthresh')


#~ opt$reflabel = 'Reference'	# label in the alignement of the strain to use for reference genome coordinates

# define path to alignment file
print(c('load genomic alignment from', opt$nfgenomeali), quote=F)
datasetname = strsplit(basename(opt$nfgenomeali), '\\.')[[1]][1]
prefull.aln = read.dna(opt$nfgenomeali, format='fasta', as.matrix=T)
reflabels.i = which(rownames(prefull.aln) in exclreflabels)
full.aln = prefull.aln[-reflabels.i,]

# compute correspondency between raw reference sequence and gapped alignement coordinates
nfmap2ref = paste(opt$resultdir, paste('coordinatesFullAln2Ref', 'RData', sep='.'), sep='')
if (file.exists(nfmap2ref)){ load(nfmap2ref)
}else{
	map.full2ref = fullAln2RefCoords(prefull.aln, opt$reflabel=opt$reflabel)
	map.full2refnona = fullAln2RefCoords(prefull.aln, opt$reflabel=opt$reflabel, bijective=F)
	save(map.full2ref, map.full2refnona, file=nfmap2ref)
}


minorallelefreqs = getMinorAlleleFreq(full.aln, multiproc=opt$threads, consider.gaps.as.variant=FALSE)
minalfrqset = paste('minalfrq', opt$min.allele.freq, sep='')
siteset = paste('biallelicsites.max', opt$max.gap, 'gaps', sep='')


## selects bi-allelic sites

# define path to file where bi-allelic site indexes are to be saved
# !!! when not satisfied by the course of the script, take care of removing saved data to mot reload them (default) next time !

nfbiali = paste(opt$resultdir, paste("BiallelicIndexes", minalfrqset, siteset, 'RData', sep='.'), sep='')
if (file.exists(nfbiali)){ 
	print(c('load \'bialraregap.i\' from', nfbiali), quote=F)
	load(nfbiali)
}else{
	print('compute \'bialraregap.i\'', quote=F)
	# biallelic sites including gaps (to be filtered out at each pairwise site comparison)
	bialgap.i = getBiAllelicIndexes(full.aln, as.logical=FALSE, minallelefreq=NULL, nonsingle=FALSE, consider.gaps.as.variant=FALSE, multiproc=opt$threads)
	# biallelic sites with at most 'maxgap' missing strains
	bialraregap.i = bialgap.i[getGaplessIndexes(full.aln[,bialgap.i], as.logical=FALSE, maxgap=opt$max.gap, gapchar=c('-', 'N', 'n'), multiproc=opt$threads)]
	save(bialgap.i, bialraregap.i, file=nfbiali)
}


## compute Nucleotidic diversity
nfnd = paste(opt$resultdir, paste("NucDiv", 'all', 'RData', sep='.'), sep='')
if (file.exists(nfnd)){ 
	print(c('load \'rollnucdiv\' from', nfnd), quote=F)
	load(nfnd)
}else{
	print('compute \'rollnucdiv\'', quote=F)
	rollnucdiv = rollStats(full.aln, windowsize=100, step=1, fun=NULL, measures=c("nucdiv"), dist.model='raw', rmAbsSeq=TRUE, propgap=0.35, relpropgap=10, multiproc=opt$threads, quiet=FALSE)
	rollnucdiv$reference.position =  map.full2ref[rollnucdiv$foci]
	rollnucdiv = rollnucdiv[!is.na(rollnucdiv$reference.position),]
	rollnucdiv$scaled.nucdiv = rollnucdiv$nucdiv / max(rollnucdiv$nucdiv, na.rm=T)
	write.table(rollnucdiv, file=paste(opt$resultdir, paste("genomic_NucDiv_win100_rmAbsSeq", "allsites", 'tab', sep='.'), sep=''), row.names=F, sep='\t')
	save(rollnucdiv, file=nfnd)
	# make a nuc.div table the size of LD tables to follow
	write.table(rollnucdiv[bialraregap.i,], file=paste(opt$resultdir, paste("genomic_NucDiv_win100_rmAbsSeq", siteset, 'tab', sep='.'), sep=''), row.names=F, sep='\t')
}

## compute LD
bial.aln = full.aln[, bialraregap.i]
ref.bial.i = map.full2ref[bialraregap.i]

nfldr = paste(opt$resultdir, nfoutldrad, 'RData', sep='.'), sep='')
if (file.exists(nfldr)){ 
	print(c('load \'bial.ldr2\' from', nfldr), quote=F)
	load(nfldr)
}else{
	nfldrtmp = paste(nfldr, 'tmp', sep='.')
	if (file.exists(nfldrtmp)){ 
		print(c('load \'lbial.ldr2\' from', nfldrtmp), quote=F)
		load(nfldrtmp)
	}else{
		print(sprintf('compute \'lbial.ldr2\' with metric \'%s\', using %d cores', opt$LD.metric, opt$threads), quote=F)
		# !!! very big memory use (ex : for 42 strains * with ~10,000 biallelic SNPs: >20Gb Mem)
		# !!! use of parallism lead to big overhead, probably due to hidden duplication of the (gigantic) data matrix ; code may be optimized, but could not find the way yet. if possible, run with nbcores=1
		# or do it, just use a high memeory machine; 100G free mem to run 16 cores over a 20,000 alignment sites works OK
		lbial.ldr2 = linkageDisequilibrium(as.character(bial.aln), opt$LD.metric, gapchar=gapchars, multiproc=opt$threads, max.dist=opt$max.dist.ldr, discard.gaps=TRUE, mem.light=TRUE, upper.triangular=FALSE)
		gc()
		# temporary save
		save(lbial.ldr2, file=nfldrtmp)
	}
	print('expand data into a matrix', quote=F)
	bial.ldr2 = as.matrix.LDlist(lbial.ldr2)
	rm(lbial.ldr2) ; gc()
	save(bial.ldr2, file=nfldr)
}

# The Fisher test p-values from a window are compared to the distribution of p-values for the genome-wide set of comparisons made in this range (i.e. within the diagonal ribbon of the square matrix of all-versus-all sites LD tests, sampled with a quantile function to get a representative sample of the same size than the window's) using a Mann-Whitney-Wilcoxon U test.
nflocld = paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep='')
if (file.exists(nflocld)){ load(nflocld)
}else{
	print("get local LD intensity", quote=F)
	if (opt$LD.metric == 'Fisher'){
#~ 		print(sprintf(" within windows of variable physical size, containing always %d biallelic sites", ldsearchpar$windowsize), quote=F)
#~ 		# for the null distribution, uses the total empiric distribution, taking only a sample to have comparable effective size for the test
#~ 		# sample the matrix in the diagonal ribbon that will be explored by the sliding window (i.e. short-distance comparisons, should be higher values than on the whole matrix)
#~ 		ribbonheight = ldsearchpar$windowsize/2
#~ 		nfribbon = paste(opt$resultdir, paste(sprintf("LDmatrix_LowerRibbonIndexes_%dsites_window%d", dim(bial.ldr2)[1], ldsearchpar$windowsize), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep='')
#~ 		if (file.exists(nfribbon)){ load(nfribbon)
#~ 		}else{
#~ 			print("  get site range for computing empirical null distribution", quote=F)
#~ 			ribbon.i = conditionalMatrixIndexes(dim(bial.ldr2)[1], getIndexesLowerRibbon, max.dist=ribbonheight, multiproc=1)
#~ 			print(sprintf("   %d sites", length(ribbon.i)), quote=F)
#~ 			save(ribbon.i, file=nfribbon)
#~ 		}
#~ 		# take a sample the size of the upper triangular part of the sliding window matrix 
#~ 		samplesize = ldsearchpar$windowsize*(ldsearchpar$windowsize - 1)/2
#~ 		# use quantiles tu ensure sample is representative of the distribution
#~ 		# excludes NA located in the lower triangular matrix so as to get the right count of informative cells
#~ 		samplewgfi = quantile(bial.ldr2[ribbon.i], p=(0:samplesize)/samplesize, na.rm=T)	
#~ 		measure = "compldfi"
#~ 		compldfi = function(alnrange){ if (length(alnrange)<5){ return(NA) }else{ wilcox.test(-log10(bial.ldr2[alnrange, alnrange]), -log10(samplewgfi), alternative='greater')$p.val }}		
#~ 		ldroll = rollStats(bial.aln, windowsize=ldsearchpar$windowsize, step=ldsearchpar$step, fun=compldfi, measures=c("compldfi"), fun.userange=list(compldfi=TRUE), multiproc=1)
#~ 		ldroll$logcompldfi = -log10(ldroll$compldfi)
		
		print(sprintf(" within windows of fixed physical size %dbp to sample bialelic sites within windows at a fixed (maximum) density of %d/window", ldsearchparsub$windowsize, ldsearchparsub$maxsize), quote=F)
		# for the null distribution, uses the total empiric distribution, taking only a sample to have comparable effective size for the test
		# sample the matrix in the diagonal ribbon that will be explored by the sliding window (i.e. short-distance comparisons, should be higher values than on the whole matrix)
		subribbonheight = ldsearchparsub$windowsize/2
		nfsubribbon = paste(opt$resultdir, paste(sprintf("LDmatrix_LowerRibbonIndexes_%dsites_physical-window%d", dim(bial.ldr2)[1], ldsearchparsub$windowsize), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep='')
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
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=compldfisub, measures=c("compldfisub"), fun.userange=list(compldfisub=TRUE), multiproc=opt$threads)
		ldrollsub$logcompldfisub = -log10(ldrollsub$compldfisub)
		ldrollsub$compldfisub.nonas  = ldrollsub$compldfisub
		ldrollsub[is.na(ldrollsub$compldfisub), 'compldfisub.nonas'] = 1
		ldrollsub$logcompldfisub.nonas = -log10(ldrollsub$compldfisub.nonas)
	}else{
		measure = "meanldr"
#~ 		meanldr = function(alnrange){ mean(bial.ldr2[alnrange, alnrange], na.rm=TRUE) }
#~ 		ldroll = rollStats(bial.aln, windowsize=ldsearchpar$windowsize, step=ldsearchpar$step, fun=meanldr, measures=c("meanldr"), fun.userange=list(meanldr=TRUE), multiproc=1)	
		meanldrsub = function(alnrange){
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
 				return(mean(bial.ldr2[bialrange, bialrange], na.rm=TRUE))
			}
		}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=meanldrsub, measures=c("meanldrsub"), fun.userange=list(meanldrsub=TRUE), multiproc=opt$threads)	
	}
#~ 	ldroll$reference.position = map.full2ref[bialraregap.i[ldroll$foci]]
	ldrollsub$reference.position = map.full2refnona[ldrollsub$foci]
	
	print("get local biallelic SNP density", quote=F)
	
	reportsnpdens = function(alnrange){ length(alnrange) }
	rollsubsnpdens = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=reportsnpdens, measures=c("reportsnpdens"), fun.userange=list(reportsnpdens=TRUE), multiproc=opt$threads)
	
#~ 	save(ldroll, ldrollsub, rollsubsnpdens, file=nflocld)
	save(ldrollsub, rollsubsnpdens, file=nflocld)
#~ 	write.table(ldroll, paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, 'tab', sep='.'), sep=''))
	write.table(ldrollsub, paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD-subsampled", minalfrqset, siteset, 'tab', sep='.'), sep=''))
}


hiLDfoci    = ldroll$reference.position[which(ldroll$compldfi < ldsearchpar$signifthresh)]
hiLDfocisub = ldrollsub$reference.position[which(ldrollsub$compldfisub < ldsearchparsub$signifthresh)]
hypervarfoci = rollnucdiv$reference.position[which(rollnucdiv$nucdiv > nucdivsearchpar$signifthresh)]

### optional if you have an annotation with the coordinates of the genes
### to determine which gene harbour the LD hospots

# load CDS/gene coordinates
if (file.exists(opt$feature.table)){	
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
	
	# get locus-wise stats
	lcds.maxnucdiv = t(sapply(lcds.ref.i, function(cds.ref.i){ 
		cds.pos.i = rollnucdiv$reference.position %in% cds.ref.i
		m = max(rollnucdiv$nucdiv[cds.pos.i], na.rm=T)
		p = rollnucdiv$reference.position[cds.pos.i & rollnucdiv$nucdiv==m][1]
		return(c(p,m))
	}))
	
	hypervarloci = names(lcds.ref.i)[which(lcds.maxnucdiv[,2] > nucdivsearchpar$signifthresh)]
	
	if (opt$LD.metric=='Fisher'){
		lcds.maxlogcompldfi = t(sapply(lcds.ref.i, function(cds.ref.i){ 
			cds.pos.i = ldroll$reference.position %in% cds.ref.i
			m = max(ldroll$logcompldfi[cds.pos.i], na.rm=T)
			p = ldroll$reference.position[cds.pos.i & ldroll$logcompldfi==m][1]
			return(c(p,m))
		}))
		lcds.maxlogcompldfisub = t(sapply(lcds.ref.i, function(cds.ref.i){ 
			cds.pos.i = ldrollsub$reference.position %in% cds.ref.i
			m = max(ldrollsub$logcompldfisub[cds.pos.i], na.rm=T)
			p = ldrollsub$reference.position[cds.pos.i & ldrollsub$logcompldfisub==m][1]
			return(c(p,m))
		}))
		
		hiLDloci = names(lcds.ref.i)[which(lcds.maxlogcompldfi[,2] > -log10(ldsearchpar$signifthresh))]
		hiLDlocisub = rownames(lcds.maxlogcompldfisub)[which(lcds.maxlogcompldfisub[,2] > -log10(ldsearchparsub$signifthresh))]

		hiLDgenes    = gnames[hiLDloci]
		hiLDgenessub =  gnames[hiLDlocisub]
		write(hiLDloci, file=paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, sprintf('hiLDloci.pval%g', ldsearchpar$signifthresh), sep='.'), sep=''))
		write(hiLDlocisub, file=paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD-subsampled", minalfrqset, siteset, sprintf('hiLDloci.pval%g', ldsearchparsub$signifthresh), sep='.'), sep=''))
		# ranked genes 
		write(unique(chompnames(hiLDloci)[order(lcds.maxlogcompldfi[hiLDloci,2], decreasing=T)]), file=paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, sprintf('hiLDgenes.pval%g', ldsearchpar$signifthresh), sep='.'), sep=''))
		write(unique(chompnames(hiLDlocisub)[order(lcds.maxlogcompldfisub[hiLDlocisub,2], decreasing=T)]), file=paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD-subsampled", minalfrqset, siteset, sprintf('hiLDgenes.pval%g', ldsearchparsub$signifthresh), sep='.'), sep=''))
	}
}else{ lcds.ref.i = NULL }

# estimate the effect of ploymorphism distribution across site
# variation of physical size of the window for a fixed number of biallelic sites
plotphysize = function(w, plotfun, ...){
	wrange = ((w/2)+1):(length(bialraregap.i)-(w/2))
	physize = sapply(wrange, function(i){ (bialraregap.i[i+(w/2)-1] - bialraregap.i[i-(w/2)]) })
	plotfun(x=map.full2ref[bialraregap.i[wrange]], y=physize, ...)
	return(physize)
}


#~ pdf(file=paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, 'pdf', sep='.'), sep=''), width=30, height=20)
pdf(file=paste(opt$resultdir, paste(sprintf("LD_%s", opt$LD.metric), "LocalLD", minalfrqset, siteset, 'pdf', sep='.'), sep=''), width=15, height=10,
 title=sprintf('%s - local LD %s - windows %dbp', datasetname, opt$LD.metric, opt$physicalwindowsize))

if (opt$LD.metric == 'Fisher'){
	### summary and genomic map plots
	par(mar=c(8,8,8,8))
	# full data scan
	plot(map.full2ref[bialraregap.i[ldroll$foci]], ldroll$logcompldfi, ylab=sprintf('LD significance in\n%d biallelic site windows [-log10(p)]', ldsearchpar$windowsize), main=paste(datasettag, "LD scan with variable-size windows", sep='\n'), xlab=sprintf('%s genome coordinates', opt$reflabel), col='white')
	abline(h=1:10, col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
	points(map.full2ref[bialraregap.i[ldroll$foci]], ldroll$logcompldfi, col=ifelse(ldroll$compldfi < ldsearchpar$signifthresh, 'red', 'black'))
	if (!is.null(lcds.ref.i) & length(hiLDgenes)>0){ text(labels=paste(hiLDgenes[!is.na(hiLDgenes)], '\'', sep='\n'), x=lcds.maxlogcompldfi[hiLDloci[!is.na(hiLDgenes)],1], y=lcds.maxlogcompldfi[hiLDloci[!is.na(hiLDgenes)],2]+.5) }
	physize = plotphysize(20, plot, type='l', xlab=sprintf('%s genome coordinates', opt$reflabel), ylab="Physical size of 20-SNP windows")
	snpdens = ldsearchpar$windowsize/(physize[ldroll$foci])
	plot(x=snpdens, y=ldroll$logcompldfi, ylab=sprintf("LD significance in\n%d biallelic site windows [-log10(p)]", ldsearchpar$windowsize), xlab="Biallelic SNP density in variable-size windows")
	abline(lm(ldroll$logcompldfi ~ snpdens), col='red')
	ct = cor.test(x=snpdens, y=ldroll$logcompldfi)
	text(x=0.8, y=max(ldroll$logcompldfi), labels=sprintf("r = %g\np = %g", ct$est, ct$p.val))
	# subsampled data scan
	plot(map.full2ref[ldrollsub$foci], ldrollsub$logcompldfisub, ylab=sprintf('LD significance in\n%dbp-wide windows [-log10(p)]', ldsearchparsub$windowsize), main=paste(datasettag, "LD scan with fixed-size windows", sep='\n'), xlab=sprintf('%s genome coordinates', opt$reflabel), col='white')
	# plot areas of NA's and low-coverage
	pb = sapply(lapply(ldrollsub$reference.position[is.na(ldrollsub$compldfisub)], function(pos){ c(pos-ldsearchparsub$windowsize/2, pos-1+ldsearchparsub$windowsize/2) }), plotbound, col='grey')
	pb = sapply(lapply(rollsubsnpdens$foci[rollsubsnpdens$reportsnpdens < ldsearchparsub$maxsize], function(pos){ c(pos-ldsearchparsub$windowsize/2, pos-1+ldsearchparsub$windowsize/2) }), plotbound, col='pink')
	abline(h=1:10, col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
	points(map.full2ref[ldrollsub$foci], ldrollsub$logcompldfisub, col=ifelse(ldrollsub$compldfisub < ldsearchparsub$signifthresh, 'red', 'black'))
	if (!is.null(lcds.ref.i) & length(hiLDgenes)>0){ text(labels=paste(hiLDgenessub[!is.na(hiLDgenessub)], '\'', sep='\n'), x=lcds.maxlogcompldfisub[hiLDlocisub[!is.na(hiLDgenessub)],1], y=lcds.maxlogcompldfisub[hiLDlocisub[!is.na(hiLDgenessub)],2]+.5) }
	legend('topright', fill=c('grey', 'pink'), legend=c('no data', sprintf("< %d biallelic SNP / window\n(low test power)", ldsearchparsub$maxsize)), bg='white')
	plot(ldrollsub$logcompldfisub ~ rollsubsnpdens$reportsnpdens, ylab=sprintf("LD significance in\n%dbp-wide windows [-log10(p)]", ldsearchparsub$windowsize), xlab="Biallelic SNP density in fixed-size windows")
	hist(rollsubsnpdens$reportsnpdens, breaks=0:20, xlab="Biallelic SNP density in fixed-size windows", main=paste(datasettag, "Distribution of SNP densities genome-wide", sep='\n'))
	qsnpdens = quantile(rollsubsnpdens$reportsnpdens, p=c(.01, .05), na.rm=T)
	abline(v=qsnpdens, col='red')
	mtext(names(qsnpdens), at=qsnpdens, side=1, col='red')
	hist(ldrollsub$logcompldfisub, breaks=0:40, xlab=sprintf("LD significance in\n%dbp-wide windows [-log10(p)]", ldsearchparsub$windowsize), main=paste(datasettag, "Distribution of LD genome-wide", sep='\n'))
	qlogcompldfisub = quantile(ldrollsub$logcompldfisub, p=c(.95, .98), na.rm=T)
	abline(v=qlogcompldfisub, col='red')
	mtext(names(qlogcompldfisub), at=qlogcompldfisub, side=1, col='red')
}else{
	### summary and genomic map plots
	par(mar=c(8,8,8,8))
	# full data scan
	plot(map.full2ref[bialraregap.i[ldroll$foci]], ldroll$meanldr, ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", ldsearchparsub$windowsize), main=paste(datasettag, "LD scan with variable-size windows", sep='\n'), xlab=sprintf('%s genome coordinates', opt$reflabel), col='white')
	abline(h=1:10*(1/10), col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
	qmeanldr = quantile(ldrollsub$meanldr, p=c(.95, .98), na.rm=T)
	points(map.full2ref[ldrollsub$foci], ldrollsub$meanldr, col=ifelse(ldrollsub$meanldr > qmeanldr[1], 'red', 'black'))
	physize = plotphysize(20, plot, type='l', xlab=sprintf('%s genome coordinates', opt$reflabel), ylab="Physical size of 20-SNP windows")
	snpdens = ldsearchpar$windowsize/(physize[ldroll$foci])
	plot(x=snpdens, y=ldroll$meanldr, ylab=sprintf("LD strength (r^2) in\n%d biallelic site windows", ldsearchpar$windowsize), xlab="Biallelic SNP density in variable-size windows")
	abline(lm(ldroll$meanldr ~ snpdens), col='red')
	ct = cor.test(x=snpdens, y=ldroll$meanldr)
	text(x=0.8, y=max(ldroll$meanldr), labels=sprintf("r = %g\np = %g", ct$est, ct$p.val))
	# subsampled data scan
	plot(map.full2ref[ldrollsub$foci], ldrollsub$meanldrsub, ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", ldsearchparsub$windowsize, ldsearchparsub$windowsize), main=paste(datasettag, "LD scan with fixed-size windows", sep='\n'), xlab=sprintf('%s genome coordinates', opt$reflabel), col='white')
	# plot areas of NA's and low-coverage
	pb = sapply(lapply(ldrollsub$reference.position[is.na(ldrollsub$compldfisub)], function(pos){ c(pos-ldsearchparsub$windowsize/2, pos-1+ldsearchparsub$windowsize/2) }), plotbound, col='grey')
	pb = sapply(lapply(rollsubsnpdens$foci[rollsubsnpdens$reportsnpdens < ldsearchparsub$maxsize], function(pos){ c(pos-ldsearchparsub$windowsize/2, pos-1+ldsearchparsub$windowsize/2) }), plotbound, col='pink')
	abline(h=1:10*(1/10), col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
	qmeanldrsub = quantile(ldrollsub$meanldrsub, p=c(.95, .98), na.rm=T)
	points(map.full2ref[ldrollsub$foci], ldrollsub$meanldrsub, col=ifelse(ldrollsub$meanldrsub > qmeanldrsub[1], 'red', 'black'))
	legend('topright', fill=c('grey', 'pink'), legend=c('no data', sprintf("< %d biallelic SNP / window\n(low test power)", ldsearchparsub$maxsize)), bg='white')
	plot(ldrollsub$meanldrsub ~ rollsubsnpdens$reportsnpdens, ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", ldsearchparsub$windowsize), xlab="Biallelic SNP density in fixed-size windows")
	hist(rollsubsnpdens$reportsnpdens, breaks=0:20, xlab="Biallelic SNP density in fixed-size windows", main=paste(datasettag, "Distribution of SNP densities genome-wide", sep='\n'))
	qsnpdens = quantile(rollsubsnpdens$reportsnpdens, p=c(.01, .05), na.rm=T)
	abline(v=qsnpdens, col='red')
	mtext(names(qsnpdens), at=qsnpdens, side=1, col='red')
	hist(ldrollsub$meanldrsub, breaks=0:40*(1/40), xlab=sprintf("LD strength (r^2) in\n%dbp-wide windows", ldsearchparsub$windowsize), main=paste(datasettag, "Distribution of LD genome-wide", sep='\n'))
	abline(v=qmeanldrsub, col='red')
	mtext(names(qmeanldrsub), at=qmeanldrsub, side=1, col='red')

}
dev.off()

# test correlation of SNP density to 
if (opt$LD.metric=='r2'){
	print(cor.test(rollnucdiv$nucdiv[rollnucdiv$foci %in% intersect(rollnucdiv$foci, ldrollsub$foci)], ldrollsub$meanldrsub[ldrollsub$foci %in% intersect(rollnucdiv$foci, ldrollsub$foci)]))
}else{
	print(cor.test(ldroll$logcompldfi[1:length(snpdens)], snpdens))
	print(cor.test(ldrollsub$logcompldfisub[1:length(snpdens)], snpdens))
	print(cor.test(ldrollsub$logcompldfisub, rollsubsnpdens$reportsnpdens))
}

## find site-to-site significant LD (use Bonferonni correction for genome-wide multiple testing)
# !!! long and memory intensive, though less than the primary LD computations as it re-usues pre-computed data

threshpval = 0.05
N = length(ref.bial.i)
m = dim(full.aln)[1] 

nfsignifpairs = paste(opt$resultdir, nfoutldrad, 'significant-sitepairspairs.RData', sep='.'), sep='')
if (file.exists(nfsignifpairs)){ load(nfsignifpairs)
}else{
	matsize =  length(which(!is.na(bial.ldr2)))	
	gc()
	if (opt$LD.metric=='r2'){
	## use the Chi-squared approximation; WRONG if any of the rare alleles at biallelic sites are observed less than 5 times, e.g. with singletons
	signifr2 = qchisq(1-(threshpval/matsize), df=1)/m
	siginfpairs = which(bial.ldr2 > signifr2, arr.ind=T)
	print(sprintf("found %d pairs of biallelic sites (out of %d) in significant linkage (p < %f after Bonferonni correction for genome-wide multiple testing)", dim(siginfpairs)[1], matsize, threshpval))
	gc()
	starttime = Sys.time()
	posr2pvals = as.data.frame(t(simplify2array(mclapply(1:dim(siginfpairs)[1], function(i){
#~ 	posr2pvals = as.data.frame(t(sapply(1:dim(siginfpairs)[1], function(i){
		biali = siginfpairs[i,]
		pos = bialraregap.i[biali]
		refpos = map.full2ref[pos]
		minfreqs = minorallelefreqs[pos]
	#~ 	r2 = lbial.ldr2[[biali[1]]][biali[2]]
		r2 = bial.ldr2[biali[1], biali[2]][[1]]
	#~ 	print(class(r2))
		pval = 1-pchisq(m*r2, df=1)
		qval = pval*matsize
		printProgressUpperMatrix(i, dim(siginfpairs)[1], step=1000, initclock=starttime) 
		return(c(refpos, pos, biali, r2, minfreqs, pval, qval))
	}, mc.cores=opt$threads, mc.preschedule=T))))
#~ 	})))
	colnames(posr2pvals) = c('ref.pos.1', 'ref.pos.2', 'aln.pos.1', 'aln.pos.2', 'bial.pos.1', 'bial.pos.2', 'r2', 'min.allele.freq.1', 'min.allele.freq.2', 'p.val', 'q.val')
	}else{
	## already get a p-value from using th Fisher exact test (recomended)
	siginfpairs = which(bial.ldr2*matsize < threshpval, arr.ind=T)
	print(sprintf("found %d pairs of biallelic sites (out of %d) in significant linkage (p < %f after Bonferonni correction for genome-wide multiple testing)", dim(siginfpairs)[1], matsize, threshpval))
	gc()
	starttime = Sys.time()
	posr2pvals = as.data.frame(t(simplify2array(mclapply(1:dim(siginfpairs)[1], function(i){
#~ 	posr2pvals = as.data.frame(t(sapply(1:dim(siginfpairs)[1], function(i){
		biali = siginfpairs[i,]
		pos = bialraregap.i[biali]
		minfreqs = minorallelefreqs[pos]
		refpos = map.full2ref[pos]
		pval = bial.ldr2[biali[1], biali[2]][[1]]
		qval = pval*matsize
		printProgressUpperMatrix(i, dim(siginfpairs)[1], step=1000, initclock=starttime) 
		return(c(refpos, pos, biali, minfreqs, pval, qval))
	}, mc.cores=opt$threads, mc.preschedule=T))))
#~ 	})))
	colnames(posr2pvals) = c('ref.pos.1', 'ref.pos.2', 'aln.pos.1', 'aln.pos.2', 'bial.pos.1', 'bial.pos.2', 'min.allele.freq.1', 'min.allele.freq.2', 'p.val', 'q.val')
	}
	gc()
	posr2pvals$site.dist = abs(posr2pvals$ref.pos.2 - posr2pvals$ref.pos.1)
	write.table(posr2pvals, file=paste(opt$resultdir, nfoutldrad, 'significant-pairs.tab', sep='.'), sep=''))
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

	pdf(file=paste(opt$resultdir, nfoutldrad, 'significant-pairs.pdf', sep='.'), sep=''), height=12, width=12)
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
	pdf(paste(opt$resultdir, paste("genome-wide_LD", minalfrqset, siteset, 'pdf', sep='.'), sep=''), width=20, height=10)
	plotLDHeatMaps(simplify2array(bial.ldr2), area.mat.i, nlabs=opt$crazy.plot, main=paste(datasettag, "correlation coefficient r^2 between sites", sep='\n'))
	dev.off()
}

