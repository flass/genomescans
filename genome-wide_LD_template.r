#!/usr/bin/Rscript --no-save --no-restore

library('ape')
library('RColorBrewer')
library('gplots')
library('ade4')

## define path to R source file
user = Sys.info()['user']
homedir = paste('/home', user, sep='/')
#~ nfsource = paste(homedir, 'scripts/misc/utils-phylo.r', sep='/')
nfsource = 'path/to/utils-phylo.r'
source(nfsource)

## define path to result directory
resultdir = 'path/to/data/and/results/directory'


## parameters
nbcores = 1		# number of cores to be used in parallel ; for LD computation on big datasets (>100kb are big), prefer use nbcores=1 (see below), unless large memory is available
print(paste("using", nbcores, "cores"), quote=F)


#~ max.dist.ldr=3000	# maximum distance between pairs of bi-allelic sites for LD computation (in mumber of intervening bi-allelic sites ; not uniform !!! polymorphism varry in density across the genome !!!) # good reason not to do it
max.dist.ldr=NULL	# compute the full matrix; not that much bigger
maxgap = 1			# maximum number of allowed missing sequences to keep a site in the alignement for LD and NucDiv computations
minalfrq = 1		# minimum allele frequency (in count of sequences) in bi-alelic sites to be retained (minalfreq = 1 => all bi-allelic sites)
gapchars = c('-', 'N', 'n')	# various characters to consider as missing data
#~ LDmetric = 'r2'
LDmetric = 'Fisher'

reflabel = 'reflabel'	# label in the alignement of the strain to use for reference genome coordinates

# define path to alignment file
full.aln = read.dna('path/to/genome_alignment.fasta', format='fasta')
minorallelefreqs = getMinorAlleleFreq(full.aln, multiproc=nbcores, consider.gaps.as.variant=FALSE)

# compute correspondency between raw reference sequence and gapped alignement coordinates
nfmap2ref = paste(resultdir, paste('coordinatesFullAln2Ref', 'RData', sep='.'), sep='/')
if (file.exists(nfmap2ref)){ load(nfmap2ref)
}else{
	map.full2ref = fullAln2RefCoords(full.aln, reflabel=reflabel)
	map.full2refnona = fullAln2RefCoords(full.aln, reflabel=reflabel, bijective=F)
	save(map.full2ref, map.full2refnona, file=nfmap2ref)
}


minalfrqset = paste('minalfrq', minalfrq, sep='')
siteset = paste('biallelicsites.max', maxgap, 'gaps', sep='')


## selects bi-allelic sites

# define path to file where bi-allelic site indexes are to be saved
# !!! when not satisfied by the course of the script, take care of removing saved data to mot reload them (default) next time !

nfbiali = paste(resultdir, paste("BiallelicIndexes", minalfrqset, siteset, 'RData', sep='.'), sep='/')
if (file.exists(nfbiali)){ 
	print(c('load \'bialraregap.i\' from', nfbiali))
	load(nfbiali)
}else{
	print('compute \'bialraregap.i\'')
	#~ bial.i = as.vector(as.numeric(sapply(readLines(paste(resultdir, rms, paste(minalfrqset, siteset, 'snp', sep='.'), sep='/')), function(x){ strsplit(x, split='\t')[[1]][4] })))
	# biallelic sites including gaps (to be filtered out at each pairwise site comparison)
	bialgap.i = getBiAllelicIndexes(full.aln, as.logical=FALSE, minallelefreq=NULL, nonsingle=FALSE, consider.gaps.as.variant=FALSE, multiproc=nbcores)
	# biallelic sites with at most 'maxgap' missing strains
	bialraregap.i = bialgap.i[getGaplessIndexes(full.aln[,bialgap.i], as.logical=FALSE, maxgap=maxgap, gapchar=c('-', 'N', 'n'), multiproc=nbcores)]
	# ALTERNATIVE  biallelic sites with at most 'maxgap' missing strains and 'maxgap' strains with N character
#~ 	bialraregap.i = bialgap.i[getGaplessIndexes(full.aln[,bialgap.i], as.logical=TRUE, maxgap=maxgappercol, gapchar=c('-', 'N', 'n'), multiproc=nbcores)) & getGaplessIndexes(full.aln[,bialgap.i], as.logical=TRUE, gapchar=c('-', 'N', 'n'), maxgap=maxgappercol, multiproc=nbcores))]
	save(bialgap.i, bialraregap.i, file=nfbiali)
}


## compute Nucleotidic diversity
nfnd = paste(resultdir, paste("NucDiv", 'all', 'RData', sep='.'), sep='/')
if (file.exists(nfnd)){ 
	print(c('load \'rollnucdiv\' from', nfnd))
	load(nfnd)
}else{
	print('compute \'rollnucdiv\'')
	rollnucdiv = rollStats(full.aln, windowsize=100, step=1, fun=NULL, measures=c("nucdiv"), dist.model='raw', rmAbsSeq=TRUE, propgap=0.35, relpropgap=10, multiproc=nbcores, quiet=FALSE)
	rollnucdiv$reference.position =  map.full2ref[rollnucdiv$foci]
	rollnucdiv = rollnucdiv[!is.na(rollnucdiv$reference.position),]
	rollnucdiv$scaled.nucdiv = rollnucdiv$nucdiv / max(rollnucdiv$nucdiv, na.rm=T)
	write.table(rollnucdiv, file=paste(resultdir, paste("genomic_NucDiv_win100_rmAbsSeq", "allsites", 'tab', sep='.'), sep='/'), row.names=F, sep='\t')
	save(rollnucdiv, file=nfnd)
	# make a nuc.div table the size of LD tables to follow
	write.table(rollnucdiv[bialraregap.i,], file=paste(resultdir, paste("genomic_NucDiv_win100_rmAbsSeq", siteset, 'tab', sep='.'), sep='/'), row.names=F, sep='\t')
}

## compute LD
bial.aln = full.aln[, bialraregap.i]
ref.bial.i = map.full2ref[bialraregap.i]
nfldr = paste(resultdir, paste(sprintf("LD_%s", LDmetric), minalfrqset, siteset, ifelse(is.null(max.dist.ldr), 'whole-matrix', paste('maxdist', max.dist.ldr, sep='')), 'RData', sep='.'), sep='/')
if (file.exists(nfldr)){ 
	print(c('load \'bial.ldr2\' from', nfldr))
	load(nfldr)
}else{
	nfldrtmp = paste(nfldr, 'tmp', sep='.')
	if (file.exists(nfldrtmp)){ 
		print(c('load \'lbial.ldr2\' from', nfldrtmp))
		load(nfldrtmp)
	}else{
		print(sprintf('compute \'lbial.ldr2\' with metric \'%s\', using %d cores', LDmetric, nbcores))
		# !!! very big memory use (ex : for 42 strains * with ~10,000 biallelic SNPs: >20Gb Mem)
		# !!! use of parallism lead to big overhead, probably due to hidden duplication of the (gigantic) data matrix ; code may be optimized, but could not find the way yet. if possible, run with nbcores=1
		# or do it, just use a high memeory machine; 100G free mem to run 16 cores over a 20,000 alignment sites works OK
		lbial.ldr2 = linkageDisequilibrium(as.character(bial.aln), LDmetric, gapchar=gapchars, multiproc=nbcores, max.dist=max.dist.ldr, discard.gaps=TRUE, mem.light=TRUE, upper.triangular=FALSE)
		gc()
		# temporary save
		save(lbial.ldr2, file=nfldrtmp)
	}
	print('expand data into a matrix')
#~ 	bial.ldr2 = as.matrix.LDlist(lbial.ldr2, length(bialraregap.i))
	bial.ldr2 = as.matrix.LDlist(lbial.ldr2)
	rm(lbial.ldr2) ; gc()
	save(bial.ldr2, file=nfldr)
}


## look for local hotspots of LD
ldsearchpar = list(20, 5, 1e-10)
# scans with windows with a fixed number of biallelic SNPs, variable physical size
names(ldsearchpar) = c('windowsize', 'step', 'signifthresh')
# scans with windows with a fixed physical size, variable number of biallelic SNPs but sub-sampled to a maximum to get the closer to a  homogeneous statistical power along the genome; any wndow with lower number of SNP than the max has a large drop in sensitivity for high LD
# must manage a trade-off between genome coverage (achived by enlarging the windowsize and lowering the maxsize) and resolution using small and SNP-dense winndows (the inverse)
ldsearchparsub = list(700, 10, 20, 1e-5)
names(ldsearchparsub) = c('windowsize', 'step', 'maxsize', 'signifthresh')
# paprmaeters window scan for nucleotidic diversity; no fancy concept here
nucdivsearchpar = list(100, 1, 0.1)
names(nucdivsearchpar) = c('windowsize', 'step', 'signifthresh')

# The Fisher test p-values from a window are compared to the distribution of p-values for the genome-wide set of comparisons made in this range (i.e. within the diagonal ribbon of the square matrix of all-versus-all sites LD tests, sampled with a quantile function to get a representative sample of the same size than the window's) using a Mann-Whitney-Wilcoxon U test.
#~ samplewgfi = as.vector(sample(bial.ldr2, ldsearchpar$windowsize^2))	# includes NA so take a sample the size of the whole matrix
nflocld = paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep="/")
if (file.exists(nflocld)){ load(nflocld)
}else{
	print("get local LD intensity", quote=F)
	if (LDmetric == 'Fisher'){
		# for the null distribution, uses the total empiric distribution, taking only a sample to have comparable effective size for the test
		# sample the matrix in the diagonal ribbon that will be explored by the sliding window (should be higher values than on the whole matrix)
		ribbonheight = ldsearchpar$windowsize/2
		nfribbon = paste(resultdir, paste(sprintf("LDmatrix_LowerRibbonIndexes_%dsites_window%d", dim(bial.ldr2)[1], ldsearchpar$windowsize), "LocalLD", minalfrqset, siteset, 'RData', sep='.'), sep="/")
		if (file.exists(nfribbon)){ load(nfribbon)
		}else{
			ribbon.i = conditionalMatrixIndexes(dim(bial.ldr2)[1], getIndexesLowerRibbon, max.dist=ribbonheight, multiproc=1)
			save(ribbon.i, file=nfribbon)
		}
		# excludes NA located in the lower triangular matrix so has take a sample the size of the upper traingular part of the sliding window matrix 
		samplesize = ldsearchpar$windowsize*(ldsearchpar$windowsize - 1)/2
		# use quantiles tu ensure sample is representative of the distribution
		samplewgfi = quantile(bial.ldr2[ribbon.i], p=(0:samplesize)/samplesize, na.rm=T)	
		measure = "compldfi"
		compldfi = function(alnrange){ if (length(alnrange)<5){ return(NA) }else{ wilcox.test(-log10(bial.ldr2[alnrange, alnrange]), -log10(samplewgfi), alternative='greater')$p.val }}		
		ldroll = rollStats(bial.aln, windowsize=ldsearchpar$windowsize, step=ldsearchpar$step, fun=compldfi, measures=c("compldfi"), fun.userange=list(compldfi=TRUE), multiproc=1)
		ldroll$logcompldfi = -log10(ldroll$compldfi)
		# use windows of fixed physical size to sample bialelic sites at a fixed (maximum) density within it
		compldfisub = function(alnrange){ 
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
#~ 				print(bialrange)
				wt = wilcox.test(-log10(bial.ldr2[bialrange, bialrange]), -log10(samplewgfi), alternative='greater')
				return(wt$p.val)
		}}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=compldfisub, measures=c("compldfisub"), fun.userange=list(compldfisub=TRUE), multiproc=1)
		ldrollsub$logcompldfisub = -log10(ldrollsub$compldfisub)
		ldrollsub$compldfisub.nonas  = ldrollsub$compldfisub
		ldrollsub[is.na(ldrollsub$compldfisub), 'compldfisub.nonas'] = 1
		ldrollsub$logcompldfisub.nonas = -log10(ldrollsub$compldfisub.nonas)
	}else{
		measure = "meanldr"
		meanldr = function(alnrange){ mean(bial.ldr2[alnrange, alnrange], na.rm=TRUE) }
		ldroll = rollStats(bial.aln, windowsize=ldsearchpar$windowsize, step=ldsearchpar$step, fun=meanldr, measures=c("meanldr"), fun.userange=list(meanldr=TRUE), multiproc=1)	
		meanldrsub = function(alnrange){
			if (length(alnrange)<5){ return(NA) 
			}else{
				bialrange = which(bialraregap.i %in% alnrange)
 				return(mean(bial.ldr2[bialrange, bialrange], na.rm=TRUE))
			}
		}
		ldrollsub = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=meanldrsub, measures=c("meanldrsub"), fun.userange=list(meanldrsub=TRUE), multiproc=1)	
	}
	ldroll$reference.position = map.full2ref[bialraregap.i[ldroll$foci]]
	ldrollsub$reference.position = map.full2refnona[ldrollsub$foci]
	
	reportsnpdens = function(alnrange){ length(alnrange) }
	rollsubsnpdens = rollStats(full.aln, subsample=list(sites=bialraregap.i, maxsize=ldsearchparsub$maxsize), windowsize=ldsearchparsub$windowsize, step=ldsearchparsub$step, fun=reportsnpdens, measures=c("reportsnpdens"), fun.userange=list(reportsnpdens=TRUE), multiproc=1)
	
	save(ldroll, ldrollsub, rollsubsnpdens, file=nflocld)
	write.table(ldroll, paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD", minalfrqset, siteset, 'tab', sep='.'), sep="/"))
	write.table(ldrollsub, paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD-subsampled", minalfrqset, siteset, 'tab', sep='.'), sep="/"))
}


hiLDfoci    = ldroll$reference.position[which(ldroll$compldfi < ldsearchpar$signifthresh)]
hiLDfocisub = ldrollsub$reference.position[which(ldrollsub$compldfisub < ldsearchparsub$signifthresh)]
hypervarfoci = rollnucdiv$reference.position[which(rollnucdiv$nucdiv > nucdivsearchpar$signifthresh)]

### optional if you have an annotation with the coordinates of the genes
### to determine which gene harbour the LD hospots

# must load CDS coordinates
#~ nfmapcds = paste('somepath', "mapCDSalignementsToRef.RData", sep="/")
nfmapcds = '/path/to/map/genes/location/to/ref/genome.RData'

### MUST generate in a separate script (or add to this one here) a list called 'lcds.ref.i' with each elelment refering to a gene, and named accordingly (unique names as they are indexes)
### each element is a vector of positions in the genomic alignment corresponding to the segment of the reference sequence where is annotated the gene, i.e. it is a subset of 'map.full2ref'

if (file.exists(nfmapcds)){
	load(nfmapcds)
	print(c('load \'lcds.ref.i\' from', nfmapcds))
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
	lcds.maxnucdiv = t(sapply(lcds.ref.i, function(cds.ref.i){ 
		cds.pos.i = rollnucdiv$reference.position %in% cds.ref.i
		m = max(rollnucdiv$nucdiv[cds.pos.i], na.rm=T)
		p = rollnucdiv$reference.position[cds.pos.i & rollnucdiv$nucdiv==m][1]
		return(c(p,m))
	}))
	
	chompnames = function(fullnames){
		# to take a subart of the full gene name; pattern to be modified depending how you call them in th first place and want them then
		gsub('^.+_(.+)_.+$', '\\1', fullnames) 
	} 
	genenames = chompnames(names(lcds.ref.i))
#~ 	genenames = names(lcds.ref.i) # or just keep as is, but won't be pretty on the plot if long
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
	
	hiLDloci = rownames(lcds.maxlogcompldfi)[which(lcds.maxlogcompldfi[,2] > -log10(ldsearchpar$signifthresh))]
	hiLDlocisub = rownames(lcds.maxlogcompldfisub)[which(lcds.maxlogcompldfisub[,2] > -log10(ldsearchparsub$signifthresh))]
	hypervarloci = rownames(lcds.maxlogcompldfi)[which(lcds.maxnucdiv[,2] > nucdivsearchpar$signifthresh)]

	hiLDgenes    = gnames[hiLDloci]
	hiLDgenessub =  gnames[hiLDlocisub]
	write(hiLDloci, file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD", minalfrqset, siteset, sprintf('hiLDloci.pval%g', ldsearchpar$signifthresh), sep='.'), sep="/"))
	write(hiLDlocisub, file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD-subsampled", minalfrqset, siteset, sprintf('hiLDloci.pval%g', ldsearchparsub$signifthresh), sep='.'), sep="/"))
	# ranked genes 
	write(unique(chompnames(hiLDloci)[order(lcds.maxlogcompldfi[hiLDloci,2], decreasing=T)]), file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD", minalfrqset, siteset, sprintf('hiLDgenes.pval%g', ldsearchpar$signifthresh), sep='.'), sep="/"))
	write(unique(chompnames(hiLDlocisub)[order(lcds.maxlogcompldfisub[hiLDlocisub,2], decreasing=T)]), file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD-subsampled", minalfrqset, siteset, sprintf('hiLDgenes.pval%g', ldsearchparsub$signifthresh), sep='.'), sep="/"))
	
}else{ lcds.ref.i = NULL }

# estimate the effect of ploymorphism distribution across site
# variation of physical size of the window for a fixed number of biallelic sites
plotphysize = function(w, plotfun, ...){
	wrange = ((w/2)+1):(length(bialraregap.i)-(w/2))
	physize = sapply(wrange, function(i){ (bialraregap.i[i+(w/2)-1] - bialraregap.i[i-(w/2)]) })
	plotfun(x=map.full2ref[bialraregap.i[wrange]], y=physize, ...)
	return(physize)
}


#~ pdf(file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD", minalfrqset, siteset, 'pdf', sep='.'), sep="/"), width=30, height=20)
pdf(file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), "LocalLD", minalfrqset, siteset, 'pdf', sep='.'), sep="/"), width=15, height=10)

### summary and genomic map plots
par(mar=c(8,8,8,8))
# full data scan
plot(map.full2ref[bialraregap.i[ldroll$foci]], ldroll$logcompldfi, ylab=sprintf('LD significance in\n%d biallelic site windows [-log10(p)]', ldsearchpar$windowsize), main="LD scan with variable-size windows", xlab=sprintf('%s genome coordinates', reflabel), col='white')
abline(h=1:10, col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
points(map.full2ref[bialraregap.i[ldroll$foci]], ldroll$logcompldfi, col=ifelse(ldroll$compldfi < ldsearchpar$signifthresh, 'red', 'black'))
if (!is.null(lcds.ref.i)){ text(labels=paste(hiLDgenes[!is.na(hiLDgenes)], '\'', sep='\n'), x=lcds.maxlogcompldfi[hiLDloci[!is.na(hiLDgenes)],1], y=lcds.maxlogcompldfi[hiLDloci[!is.na(hiLDgenes)],2]+.5) }
physize = plotphysize(20, plot, type='l', xlab=sprintf('%s genome coordinates', reflabel), ylab="Physical size of 20-SNP windows")
snpdens = ldsearchpar$windowsize/(physize[ldroll$foci])
plot(x=snpdens, y=ldroll$logcompldfi, ylab=sprintf("LD significance in\n%d biallelic site windows [-log10(p)]", ldsearchpar$windowsize), xlab="Biallelic SNP density in variable-size windows")
abline(lm(ldroll$logcompldfi ~ snpdens), col='red')
ct = cor.test(x=snpdens, y=ldroll$logcompldfi)
text(x=0.8, y=max(ldroll$logcompldfi), labels=sprintf("r = %g\np = %g", ct$est, ct$p.val))
# subsampled data scan
plot(map.full2ref[ldrollsub$foci], ldrollsub$logcompldfisub, ylab=sprintf('LD significance in\n%dbp-wide windows [-log10(p)]', ldsearchparsub$windowsize), main="LD scan with fixed-size windows", xlab=sprintf('%s genome coordinates', reflabel), col='white')
# plot areas of NA's and low-coverage
pb = sapply(lapply(ldrollsub$reference.position[is.na(ldrollsub$compldfisub)], function(pos){ c(pos-ldsearchparsub$windowsize/2, pos-1+ldsearchparsub$windowsize/2) }), plotbound, col='grey')
pb = sapply(lapply(rollsubsnpdens$foci[rollsubsnpdens$reportsnpdens < ldsearchparsub$maxsize], function(pos){ c(pos-ldsearchparsub$windowsize/2, pos-1+ldsearchparsub$windowsize/2) }), plotbound, col='pink')
abline(h=1:10, col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
points(map.full2ref[ldrollsub$foci], ldrollsub$logcompldfisub, col=ifelse(ldrollsub$compldfisub < ldsearchparsub$signifthresh, 'red', 'black'))
if (!is.null(lcds.ref.i)){ text(labels=paste(hiLDgenessub[!is.na(hiLDgenessub)], '\'', sep='\n'), x=lcds.maxlogcompldfisub[hiLDlocisub[!is.na(hiLDgenessub)],1], y=lcds.maxlogcompldfisub[hiLDlocisub[!is.na(hiLDgenessub)],2]+.5) }
legend('topright', fill=c('grey', 'pink'), legend=c('no data', sprintf("< %d biallelic SNP / window\n(low test power)", ldsearchparsub$maxsize)), bg='white')
plot(ldrollsub$logcompldfisub ~ rollsubsnpdens$reportsnpdens, ylab=sprintf("LD significance in\n%dbp-wide windows [-log10(p)]", ldsearchparsub$windowsize), xlab="Biallelic SNP density in fixed-size windows")
hist(rollsubsnpdens$reportsnpdens, breaks=0:20, xlab="Biallelic SNP density in fixed-size windows", main="Distribution of SNP densities genome-wide")
qsnpdens = quantile(rollsubsnpdens$reportsnpdens, p=c(.01, .05), na.rm=T)
abline(v=qsnpdens, col='red')
mtext(names(qsnpdens), at=qsnpdens, side=1, col='red')
hist(ldrollsub$logcompldfisub, breaks=0:40, xlab=sprintf("LD significance in\n%dbp-wide windows [-log10(p)]", ldsearchparsub$windowsize), main="Distribution of LD genome-wide")
qlogcompldfisub = quantile(ldrollsub$logcompldfisub, p=c(.95, .98), na.rm=T)
abline(v=qlogcompldfisub, col='red')
mtext(names(qlogcompldfisub), at=qlogcompldfisub, side=1, col='red')
dev.off()

# test correlation of SNP density to 
if (LDmetric=='r2'){
	print(cor.test(rollnucdiv$nucdiv[rollnucdiv$foci %in% intersect(rollnucdiv$foci, ldrollsub$foci)], ldrollsub$meanldrsub[ldrollsub$foci %in% intersect(rollnucdiv$foci, ldrollsub$foci)]))
}else{
	print(cor.test(ldroll$logcompldfi[1:length(snpdens)], snpdens))
	print(cor.test(ldrollsub$logcompldfisub[1:length(snpdens)], snpdens))
	print(cor.test(ldrollsub$logcompldfisub, rollsubsnpdens$reportsnpdens))
}

## find site-to-site significant LD (use Bonferonni correction for genome-wide multiple testing)
# !!! long and memory intensive, though les than the promary LD computations (can be run on 20G where the previous neded 80G)

threshpval = 0.05
N = length(ref.bial.i)
m = dim(full.aln)[1] 

nfsignifpairs = paste(resultdir, paste(sprintf("LD_%s", LDmetric), minalfrqset, siteset, ifelse(is.null(max.dist.ldr), 'whole-matrix', paste('maxdist', max.dist.ldr, sep='')), 'significant-sitepairspairs.RData', sep='.'), sep='/')
if (file.exists(nfsignifpairs)){ load(nfsignifpairs)
}else{
	#~ matsize = (N*(N-1))/2
	# accounts for discarded comparisons due to non-overlap of nformative (non-gap) rows in some pairwise site comparisons
	matsize =  length(which(!is.na(bial.ldr2)))	
	gc()
	if (LDmetric=='r2'){
	## use the Chi-squared approximation; WRONG if any of the rare alleles at biallelic sites are observed less than 5 times, e.g. with singletons
	signifr2 = qchisq(1-(threshpval/matsize), df=1)/m
	#~ siginfpairs = lapply(1:length(lbial.ldr2), function(i){ which(lbial.ldr2[[i]] > signifr2) })
	siginfpairs = which(bial.ldr2 > signifr2, arr.ind=T)
	print(sprintf("found %d pairs of biallelic sites (out of %d) in significant linkage (p < %f after Bonferonni correction for genome-wide multiple testing)", dim(siginfpairs)[1], matsize, threshpval))
	gc()
	#~ posr2pvals = as.data.frame(t(simplify2array(mclapply(1:dim(siginfpairs)[1], function(i){
	starttime = Sys.time()
	posr2pvals = as.data.frame(t(sapply(1:dim(siginfpairs)[1], function(i){
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
	#~ }, mc.cores=nbcores, mc.preschedule=T))))
	})))
	colnames(posr2pvals) = c('ref.pos.1', 'ref.pos.2', 'aln.pos.1', 'aln.pos.2', 'bial.pos.1', 'bial.pos.2', 'r2', 'min.allele.freq.1', 'min.allele.freq.2', 'p.val', 'q.val')
	}else{
	## already get a p-value from using th Fisher exact test (recomended)
	siginfpairs = which(bial.ldr2*matsize < threshpval, arr.ind=T)
	print(sprintf("found %d pairs of biallelic sites (out of %d) in significant linkage (p < %f after Bonferonni correction for genome-wide multiple testing)", dim(siginfpairs)[1], matsize, threshpval))
	gc()
	#~ posr2pvals = as.data.frame(t(simplify2array(mclapply(1:dim(siginfpairs)[1], function(i){
	starttime = Sys.time()
	posr2pvals = as.data.frame(t(sapply(1:dim(siginfpairs)[1], function(i){
		biali = siginfpairs[i,]
		pos = bialraregap.i[biali]
		minfreqs = minorallelefreqs[pos]
		refpos = map.full2ref[pos]
		pval = bial.ldr2[biali[1], biali[2]][[1]]
		qval = pval*matsize
		printProgressUpperMatrix(i, dim(siginfpairs)[1], step=1000, initclock=starttime) 
		return(c(refpos, pos, biali, minfreqs, pval, qval))
	#~ }, mc.cores=nbcores, mc.preschedule=T))))
	})))
	colnames(posr2pvals) = c('ref.pos.1', 'ref.pos.2', 'aln.pos.1', 'aln.pos.2', 'bial.pos.1', 'bial.pos.2', 'min.allele.freq.1', 'min.allele.freq.2', 'p.val', 'q.val')
	}
	gc()
	posr2pvals$site.dist = abs(posr2pvals$ref.pos.2 - posr2pvals$ref.pos.1)
	write.table(posr2pvals, file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), minalfrqset, siteset, ifelse(is.null(max.dist.ldr), 'whole-matrix', paste('maxdist', max.dist.ldr, sep='')), 'significant-pairs.tab', sep='.'), sep='/'))
	# for all comparisons
	if (LDmetric=='r2'){
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
# for only the significant comparisons
intersitedistclases = seq(from=0, to=floor(max(posr2pvals$site.dist/1000, na.rm=T)), by=1)
qvalbyintersitedistclases = lapply(intersitedistclases, function(k){ log10(posr2pvals$q.val[posr2pvals$site.dist >= k*1000 & posr2pvals$site.dist < (k+1)*1000]) })
names(qvalbyintersitedistclases) = sapply(intersitedistclases, function(k){ paste(k, k+1, sep='-') })
#~ sitedistbyqvalclass = sapply(seq(from=floor(min(allqvals, na.rm=T)), to=floor(max(allqvals, na.rm=T)), by=0.5), function(k){ length(which(allqvals>=k & allqvals<k+0.5)) })

pdf(file=paste(resultdir, paste(sprintf("LD_%s", LDmetric), minalfrqset, siteset, ifelse(is.null(max.dist.ldr), 'whole-matrix', paste('maxdist', max.dist.ldr, sep='')), 'significant-pairs.pdf', sep='.'), sep='/'), height=12, width=12)
barplot(log10(freqqvals), names.arg=paste('[', rangeqvals, '; ', rangeqvals+0.5, ']', sep=''), ylim=c(0,10), las=2, xlab='log10(corrected p-values)', ylab='log10(frequency)',  cex.lab=1.5,main='Distribution of p-values of Fisher\'s exact test for LD significance')
hist(posr2pvals$site.dist, breaks=(0:24)*10000, cex.lab=1.5, xlab='inter-site distance')
plot(x=((1:100)/100), y=quantile(posr2pvals$site.dist, p=(1:100)/100), xlab='quantiles', ylab='inter-site distance',  cex.lab=1.5,main=sprintf('Cumulative distribution function of inter-site distances\nfor significantly linked site (p < %g after Bonferonni cor.)', threshpval))
boxplot(qvalbyintersitedistclases, xlab='inter-site distance (kb)', ylab='log10(corrected p-values)', cex.lab=1.5, main='Distribution of significant p-values with site distance')
#~ hist(log10(posr2pvals$q.val), breaks=101)
#~ plot(log10(posr2pvals$q.val) ~ posr2pvals$site.dist)
dev.off()



## plots the LD matrix (the one computed with potentially included rare gaps) to PDF
# USE ONLY ON SMALL MATRICES !!!!!!
area.mat.i = map.full2ref[bialraregap.i]
nlabs = 200 # number of evenly spaced site coordinates to plot
pdf(paste(resultdir, paste("genome-wide_LD", minalfrqset, siteset, 'pdf', sep='.'), sep='/'), width=20, height=10)
plotLDHeatMaps(simplify2array(bial.ldr2), area.mat.i, nlabs=100, main="correlation coefficient r^2 between sites")
dev.off()

### deprecated code, most functions from external module 'util-phylo.r' should still work though

## export biallelic sites to Plink format readable by Haploview, Plink or EigenStrat (software for analysis of structure of population, mainly using a PCA, Price et al. Nat Genet 38, 904â€“909 (2006), needs a modified format)
# ONLY WORKS ON GAPLESS SITES
#sitesetno = 'biallelicsites.nogaps'
#nfbialnoi = paste(resultdir, paste("BiallelicIndexes", minalfrqset, sitesetno, 'RData', sep='.'), sep='/')
#if (file.exists(nfbialnoi)){ 
#	print(c('load \'bialnogap.i\' from', nfbialnoi))
#	load(nfbiali)
#}else{
#	print('compute \'bialnogap.i\'')
#	bialnogap.i = bialgap.i[getGaplessIndexes(full.aln[,bialgap.i], as.logical=FALSE, maxgap=0, gapchar=c('-', 'N', 'n'), multiproc=nbcores)]
#	save(bialnogap.i, file=nfbialnoi)
#}
#geno.aln = convertTo01genotypes(full.aln[,bialnogap.i], get.biallelic=F, get.gapless=F, multiproc=nbcores)
#exportToPED(full.aln, paste(resultdir, paste(minalfrqset, sitesetno, 'ped', sep='.'), sep='/'), sep=' ', get.biallelic=F, get.gapless=F, STRUCTUREfmt=T, multiproc=6)
#write.table(t(geno.aln), file=paste(resultdir, paste(minalfrqset, sitesetno, 'eigenstratgeno', sep='.'), sep='/'), sep='', row.names=F, col.names=F, quote=F)
#exportSNPpos(bialnogap.i, paste(resultdir, paste(minalfrqset, sitesetno, 'snp', sep='.'), sep='/'), EIGEINSTRATfmt=T)
#exportSNPpos(bialnogap.i, paste(resultdir, paste(minalfrqset, sitesetno, 'pedsnp', sep='.'), sep='/'), EIGEINSTRATfmt=F)
# if you have a table of strain origins called 'strain.ori'
#~ write.table(cbind(rownames(full.aln), rep('U', dim(full.aln)[1]), as.character(strain.ori[rownames(full.aln), 'country'])),
#~  file=paste(resultdir, paste(minalfrqset, sitesetno, 'ind', sep='.'), sep='/'), row.names=F, col.names=F, quote=F)

# site-wise average LD map at varying distances
#coul = c('black', 'blue', 'purple', 'red', 'orange')
#kbs = c(0, 1, 4, 9)
#kbe = c(1, 2, 5, 10)
#
#ref.bial.i = map.full2ref[bialraregap.i]
#
## adapt granularity of data considered for plotting
##~ grain = ifelse(max(ref.bial.i, na.rm=T)>100000, 10000, ifelse(max(ref.bial.i, na.rm=T)>20000, 5000, 2000))
#grain = 1000 # (set to 100 or 1000 for real genome wide plotting)
#pdf(paste(resultdir, paste("genomic_LD_varying_dist", minalfrqset, siteset, 'pdf', sep='.'), sep='/'), width=20, height=10)
#N = length(ref.bial.i)
#for (k in 1:length(kbs)){
#	starttime = Sys.time()
#	gc()
#	print(paste(kbs[k], 'to', kbe[k], 'kb away'))
#	kbmeanld = simplify2array(mclapply(1:length(ref.bial.i), function(i){ 
## 	kbmeanld = sapply(1:N, function(i){ 
#		pair.range = which(!sapply(1:N, getIndexesOutDiagonalRibbon, i=i, max.dist=max.dist.ldr)) 
#		d = abs(ref.bial.i[i] - ref.bial.i[pair.range])	# physical distance map the size of the ith element of LD list
#		stopifnot(length(bial.ldr2[[i]])==length(d))
#		kbdisti = which(d>=kbs[k]*1000 & d<kbe[k]*1000)
#		printProgressUpperMatrix(i, N, step=100, initclock=starttime)
#		mean(bial.ldr2[[i]][kbdisti], na.rm=T) 
## 	})
#	}, mc.cores=nbcores, mc.preschedule=T))	# parallel implementation ; switch to simple 'sapply' loop if problematic
#	
#	# plot intensity of LD  ; maybe heavy to compute
#	if (k==1){
#		plot(ref.bial.i, kbmeanld, col=coul[k], xaxt='n')
#		axis(side=1, at=seq(from=0, to=max(ref.bial.i)%/%grain, by=grain))
#	}else{
#		points(ref.bial.i, kbmeanld, col=coul[k])
#	}
#	tabkbmeanld = cbind(ref.bial.i, kbmeanld) ; colnames(tabkbmeanld) = c('position', 'mean.ld')
#	write.table(tabkbmeanld, file=paste(resultdir, paste("genomic_LD_dist", minalfrqset, siteset, paste(kbs[k], '-', kbe[k], 'kb', sep=''), 'tab', sep='.'), sep='/'),
#	 row.names=F, sep='\t')
#}
#legend('topright', legend=paste(kbs, 'to', kbe, 'kb away'), fill=coul)
#dev.off()
#

