#!/usr/bin/Rscript --no-save --no-restore
library(ade4)
library(RColorBrewer)
library(gplots)

options(width = 120)

minsizclades = 3

homedir = Sys.getenv()['HOME']

source(paste(homedir, 'scripts/misc/utils-phylo.r', sep='/'))

printinfomatrix = function(m){
	mname = deparse(substitute(m))
	print(sprintf("class(%s): %s ; dim(%s): %s", mname, class(m), mname, paste(dim(m), collapse=', ')), quote=F)
	print(sprintf("colnames(%s):", mname), quote=F)
	print(colnames(m))
}

carg = commandArgs(trailingOnly=TRUE)
ordtag = '--gene.list.is.ordered'
ordgenelist = ordtag %in% carg
if (ordgenelist){
	print("# gene list is provided in the genome order; will not reorder it", quote=F)
	carg = carg[-which(carg==ordtag)]
}
newtag = '--new'
recompute = newtag %in% carg
if (recompute){
	print("# ignore formerly saved data archive and recompute all matrices", quote=F)
	carg = carg[-which(carg==newtag)]
}
print(c("# arguments:", carg), quote=F)
dircladeprofile = carg[1]
print(paste('read data in', dircladeprofile))
if (length(carg)>1){
	nfspecialgenes = carg[2]
}else{ nfspecialgenes = NULL }
if (length(carg)>2){
	nfstrainorigin = carg[3]
}else{ nfstrainorigin = NULL }
if (length(carg)>3){
	focal.strains = strsplit(carg[4], split=',')[[1]]
}else{ focal.strains = NULL }

# bipart description
dbbiparts = read.table(paste(dircladeprofile, 'bipart_db.tab', sep='/'), h=T, colClasses=c('character', rep('numeric', 4)))
# bipart posterior probabilities along genes
PPbiparts = as.matrix(read.table(paste(dircladeprofile, 'bipart_PostProbs.tab', sep='/'), h=T, row.names=1))
# bipart occurence and degenerate matching along genes
patbiparts = as.matrix(read.table(paste(dircladeprofile, 'bipart_clusters.tab', sep='/'), row.names=1, h=T))
# gene tree details
gtdet = read.table(paste(dircladeprofile, 'bipart_intrees.tab', sep='/'), h=T)
# ordered list of taxa
ltax = scan(paste(dircladeprofile, 'taxlabels', sep='/'), what='character')
# generic bipartition/gene labels
nrbiparts = as.numeric(rownames(patbiparts))

nfglist = paste(dircladeprofile, 'screened_gene_list.txt', sep='/')
if (!recompute & file.exists(nfglist)){
	genenames = gsub('-', '.', readLines(nfglist))
}else{
	genenames = gsub('\\.', '-', colnames(patbiparts))
}
if (!ordgenelist){
	genecoords = as.data.frame(t(sapply(genenames, function(x){ as.numeric(strsplit(strsplit(x, split='_')[[1]][2], split='-')[[1]]) })))
	genecoordorder = order(genecoords[,1])
	genenames = genenames[genecoordorder]
	PPbiparts = PPbiparts[,genecoordorder]
}
print("genenames:", quote=F)
print(genenames, quote=F)

# make short gene names

## may have to tweak this to parse correctly gene names, depending on their format

gnames = sapply(genenames, function(x){ strsplit(x, split='_')[[1]][1] })
#~ gnames = sapply(genenames, function(x){ strsplit(x, split='_')[[1]][3] })
#~ print(gnames)
gcount = table(gnames)
#~ for (g in names(gcount)[gcount>1]){
for (g in names(gcount)){
	k = 1 
	gi = which(gnames==g)
	gstart = sapply(genenames[gi], function(x){ as.numeric(strsplit(strsplit(x, split='_')[[1]][3], split='\\.')[[1]][1]) })
	for (ign in gi[order(gstart)]){
		gnames[ign] = paste(gnames[ign], k, sep='-')
		k = k + 1
	}
}
#~ gnames = genenames
rownames(dbbiparts) = paste('c', dbbiparts$bipart, sep='')
colnames(patbiparts) = colnames(PPbiparts) = as.vector(gnames[genenames])
printinfomatrix(PPbiparts)
patbiparts.bin = !apply(patbiparts, 2, is.na)
cnrbiparts = rownames(patbiparts.bin) = rownames(patbiparts) = rownames(PPbiparts) = paste('c', nrbiparts, sep='')

# distance between bipartitions

bips = t(sapply(dbbiparts$profile[nrbiparts], function(x){ as.numeric(strsplit(x, split='')[[1]]) }))
rownames(bips) = rownames(patbiparts)

# combine PP with presence/absence of biparts in consensus tree
bipconsabspres = sapply(as.vector(genenames), function(g){as.numeric(nrbiparts %in% gtdet[gtdet$gene_label==g,'bipart'])})
colnames(bipconsabspres) = as.vector(gnames[genenames])
consabspresPPbiparts = bipconsabspres + PPbiparts
rownames(consabspresPPbiparts) = rownames(PPbiparts)
colnames(consabspresPPbiparts) = colnames(PPbiparts)
printinfomatrix(consabspresPPbiparts)

nfmatcompbip = paste(dircladeprofile, 'bipartition_compatibility_matrix.RData', sep='/')
if (!recompute & file.exists(nfmatcompbip)){ load(nfmatcompbip)
}else{
	# matrix of compatible biparts
	bipinPPmat = which(paste('c', dbbiparts$bipart, sep='') %in% rownames(PPbiparts))
	bipinPPmat = sort(bipinPPmat)
	names(bipinPPmat) = paste('c', bipinPPmat, sep='')
	matcompbip = compatibleBiparts(dbbiparts$profile[bipinPPmat], ltax, nbcores=6)
	dimnames(matcompbip) = list(names(bipinPPmat), names(bipinPPmat))
	save(matcompbip, file=nfmatcompbip)
}
# matrix of probability of rejection of one bipartiton in one gene
# = the maximal posterior probability in a gene of all bipartitions uncompatible with the focal bipartition
nfcompPPbip = paste(dircladeprofile, 'bipart_Compat.tab', sep='/')
if (!recompute & file.exists(nfcompPPbip)){ compPPbiparts = read.table(nfcompPPbip)
}else{
	rejPPbiparts = data.matrix(t(sapply(rownames(PPbiparts), function(b){
		uncompbip = colnames(matcompbip)[!matcompbip[b,]]
		sapply(colnames(PPbiparts), function(g){
			max(c(0, PPbiparts[uncompbip,g]), na.rm=T)
		})
	})))
	compPPbiparts = 1 - rejPPbiparts
	colnames(compPPbiparts) = as.vector(gnames[genenames])
	write.table(compPPbiparts, file=nfcompPPbip, row.names=T, sep='\t')
}
printinfomatrix(compPPbiparts)
# combine compatibilty measure with presence/absence of biparts in consensus tree
consabsprescompPPbiparts = bipconsabspres + compPPbiparts
rownames(consabsprescompPPbiparts) = rownames(compPPbiparts)
colnames(consabsprescompPPbiparts) = colnames(compPPbiparts)
printinfomatrix(consabsprescompPPbiparts)


plotHeatMapBipartGenomeDensity = function(mat, bip.filter=NULL, gene.filter=NULL, gene.names=gnames, main=NULL, mark.consensus.bipart=NULL,
 minsumPPovergenes=1, minPPamonggenes=NULL, heatvariable='posterior probability', densfun=function(x){sum(x, na.rm=T)},
 lab.bip=NULL, dendro.dist.bips=FALSE, heatmapengine='heatmap', col.palette=NULL, col.na=par('bg'), side.cols=NULL){
	print(main)
	
	if (is.null(bip.filter)){ bf = rep(TRUE, dim(mat)[1]) }else{ bf = bip.filter }
	if (is.null(gene.filter)){ gf = rep(TRUE, dim(mat)[2]) }else{ gf = gene.filter }
	
	if (is.null(mark.consensus.bipart)){
		bipcols = c(brewer.pal(9, 'Blues'))
		m = mat
	}else{
		bipcols = c(brewer.pal(9, 'Blues'), brewer.pal(9, 'Reds'))
		m = mark.consensus.bipart
	}
	# override color settings ; caution when mixing with 'mark.consensus.bipart'
	if (!is.null(col.palette)){ bipcols = col.palette }
	
	bipdensity = apply(mat[,gf], 1, densfun)
	
	if (!is.null(minsumPPovergenes)){ denseprofile = bipdensity >= minsumPPovergenes
	}else{ denseprofile = rep(TRUE, dim(mat)[1]) }
	if (!is.null(minPPamonggenes)){ minprofile = apply(mat[,gf], 1, function(x){ any(x>minPPamonggenes, na.rm=T) })
	}else{ minprofile = rep(TRUE, dim(mat)[1]) }
	B = denseprofile & minprofile & bf
	
	denscols = c("#FFFFFF", brewer.pal(9, 'Reds'))
	mindens = min(min(bipdensity[B], na.rm=T), 0, na.rm=T)
	
	ngenecol = 5
	geneinfo = apply(mat[B,gf], 2, densfun)
	infocols = c("#FFFFFF", brewer.pal(ngenecol, 'Reds'))
	if (!is.null(side.cols)){
		bipdenscols = sapply(bipdensity[B], function(s){ s = s - mindens ; denscols[ifelse(s>=10, 10, floor(s)+1)] })
		geneinfocols = sapply(geneinfo, function(s){ infocols[ifelse(s<=1, 1, ifelse(s%/%4 >= ngenecol, ngenecol+1, (s%/%4)+2))] })
	}else{
		bipdenscols = geneinfocols = NULL
	}
	if (dendro.dist.bips){
		distbips = dist(bips[B,], method='manhattan')
		dendrobip = as.dendrogram(hclust(distbips, 'average'))
	}else{ dendrobip = NA }
	if (is.null(lab.bip)){ lR = rownames(m)[B]
	}else{ lR = lab.bip[B] }
	
	if (!is.null(minsumPPovergenes)){
		ylabel = paste("bipartitions with cumulative", heatvariable, ">=", minsumPPovergenes)
	}else{ if (!is.null(minPPamonggenes)){
		ylabel = paste("bipartitions with", heatvariable, ">=", minPPamonggenes, "in at least one gene")
	}else{
		ylabel = "bipartitions"
	}}
	if (is.null(main)){ main = heatvariable }
	
	heatmapfun = get(heatmapengine, inherits=T)
	print(dim(m[B,gf]))
	heatmapfun(data.matrix(m[B,gf]), Rowv=dendrobip, Colv=NA, revC=(heatmapengine=='heatmap' & !dendro.dist.bips), scale='none', col=bipcols, labCol=gene.names[gf], na.color=col.na, labRow=lR,
	 main=main, xlab="Loci", ylab=ylabel, cex.main=2, cex.lab=2,
	 trace='none', dendrogram='none')
#~ 	 RowSideColors=bipdenscols, ColSideColors=geneinfocols,
}

plotHMBGDbyGenes = function(lPPB, lcabPPB, cherries, cherrylabs, gene.filter=NULL, focus.strains=NULL){
for (nlPPb in names(lPPb)){
		mainlab = paste(nlPPb, 'of bipartitions across genes')
		# 2-strain clades ("cherries")
		plotHeatMapBipartGenomeDensity(lPPb[[nlPPb]], gene.filter=gene.filter, bip.filter=cherries, gene.names=gnames[genenames], main=paste(mainlab, "(strain pairs)"),
		 mark.consensus.bipart=lcabPPb[[nlPPb]], minsumPPovergenes=1)
		# 2-strain clades ("cherries") with strain names
		plotHeatMapBipartGenomeDensity(lPPb[[nlPPb]], gene.filter=gene.filter, bip.filter=cherries, gene.names=gnames[genenames], main=paste(mainlab, "(strain pairs)"),
		 mark.consensus.bipart=lcabPPb[[nlPPb]], minsumPPovergenes=1, lab.bip=sapply(cherrylabs, paste, collapse=', '))
		if (!is.null(focus.strains)){
			for (fs in focus.strains){
				fsi = sapply(cherrylabs, function(x){ fs %in% x })
				plotHeatMapBipartGenomeDensity(lPPb[[nlPPb]], gene.filter=gene.filter, bip.filter=(cherries & fsi), gene.names=gnames[genenames], main=paste(mainlab, " (strain pairs with ", fs,")", sep=''),
				 mark.consensus.bipart=lcabPPb[[nlPPb]], minsumPPovergenes=1, lab.bip=sapply(cherrylabs, paste, collapse='\n'))
			}
		}
		 
		# "larger" clades
		plotHeatMapBipartGenomeDensity(lPPb[[nlPPb]], gene.filter=gene.filter, bip.filter=largeclades, gene.names=gnames[genenames], main=paste(mainlab, "(clades of >= 3 strains)"), mark.consensus.bipart=lcabPPb[[nlPPb]], minsumPPovergenes=1)
	}
}
# PP support of bipartition accross genes

# clades with >= 3 sample/strains
largeclades = dbbiparts$size_smallpart[as.numeric(nrbiparts)] >= minsizclades
# 2-strain clades ("cherries")
cherries = dbbiparts$size_smallpart[as.numeric(nrbiparts)] == 2
# row labelling showing the small clade's sample names for biparts #cherry [cherries]
cherrylabs = sapply(dbbiparts$profile[as.numeric(nrbiparts)], translateBipartToTaxSet, ltax=ltax, value='smallclade', as.vector=TRUE)


dbbips = dbbiparts[as.numeric(nrbiparts),]
dbbips$labels  = as.vector(sapply(cherrylabs, paste, collapse=', '))
write.table(dbbips[,c('bipart', 'labels')], file=paste(dircladeprofile, 'bipart_labels.tab', sep='/'), sep='\t', quote=T, row.names=F)


pdf(paste(dircladeprofile, 'bipart-profile.pdf', sep='/'), width=20, height=20)

	## PP support
lPPb = list(PPbiparts, compPPbiparts) ; names(lPPb) = c('Support', 'Compatibility')
	## compatibility = 1 - PP support against the bipartition
lcabPPb = list(consabspresPPbiparts, consabsprescompPPbiparts) ; names(lcabPPb) = c('Support', 'Compatibility')
## all genes
plotHMBGDbyGenes(lPPB, lcabPPB, cherries, cherrylabs, focus.strains=focal.strains)
## sub-sample of genes
if (!is.null(nfspecialgenes)){
	specialgenes = sort(readLines(nfspecialgenes))
	spegenes = gnames[specialgenes]
	spegenefi = colnames(patbiparts) %in% spegenes
	plotHMBGDbyGenes(lPPB, lcabPPB, cherries, cherrylabs, gene.filter=spegenefi)
}
dev.off()

gapsizes = 0:2
ppthresh = 0.35
locationtag = 'country'

pdf(paste(dircladeprofile, 'distrib_cherry_track_sizes.pdf', sep='/'), width=10, height=7)
lPPb = list(PPbiparts, compPPbiparts) ; names(lPPb) = c('supported', 'compatible')
congruenttracks = lapply(names(lPPb), function(nlPPb){
	lapply(gapsizes, function(gapsize){

		# draw distribution of sizes of tracks of genes with common pairing of strains
		suptrackcherries = getSupportedBlocks(lPPb[[nlPPb]], subset.bip=which(cherries), gapsize=gapsize, ppthresh=ppthresh)
		sizesuptrack = suptrackcherries$track.size
		h = hist(sizesuptrack, xlab='number of genes', breaks=0:max(sizesuptrack),
		 main=paste('distribution of sizes of tracks of genes\nwith common ', nlPPb, ' pairing of strains\n(min. support threshold = ', ppthresh, '; max. gap size = ', gapsize, ')', sep=''))
		text(x=0.9*max(sizesuptrack), y=0.8*max(h$counts), adj=1, labels=paste('mean =', as.character(signif(mean(sizesuptrack), digits=3))))
		hist(sizesuptrack[sizesuptrack>=2], xlab='number of genes', breaks=0:max(sizesuptrack),
		 main=paste('distribution of sizes of tracks of genes\nwith common ', nlPPb, ' pairing of strains\n(min. support threshold = ', ppthresh, '; max. gap size = ', gapsize, ')', sep=''))

		# list all tracks of genes with common clades
		suptrack = getSupportedBlocks(lPPb[[nlPPb]], gapsize=gapsize, ppthresh=ppthresh, gene.names=gnames[genenames], min.report=5)
		suptrack = merge(suptrack, dbbiparts[dbbiparts$size_smallpart>1,], by.x='bipart', by.y='row.names')
		suptrack$smallclade = sapply(as.character(suptrack$profile), translateBipartToTaxSet, ltax, value='smallclade', as.vector=FALSE)
		write.table(suptrack[order(suptrack$track.size, decreasing=T),], file=paste(dircladeprofile, paste('gene_blocks_with_common_', nlPPb, '_splits_', paste('gapsize', gapsize, sep='='), '.tab', sep=''), sep='/'), sep='\t', row.names=F)
		 
		# try to relate the blocks of strain paring to their sampling provenance
		if (!is.null(nfstrainorigin)){
			strainorigin = read.table(nfstrainorigin, sep='\t', header=T, colClasses='character')
			suptrackcherries = merge(suptrackcherries, dbbiparts[dbbiparts$size_smallpart==2,], by.x='bipart', by.y='row.names')
			suptrackcherries$smallclade = sapply(as.character(suptrackcherries$profile), translateBipartToTaxSet, ltax, value='smallclade', as.vector=FALSE)
			countorigins = sapply(suptrackcherries$smallclade, function(s){
				length(unique(strainorigin[strainorigin[,'seqacc'] %in% strsplit(s, split=', ')[[1]], locationtag]))
			})
			print(table(countorigins))
			suptrackcherries$sameorigin = ifelse(countorigins==1, paste('same', locationtag), paste('different', locationtag))
			sizesuptrack = suptrackcherries$track.size
			bp = boxplot(sizesuptrack ~ suptrackcherries$sameorigin, main=paste('Blocks of genes with common', nlPPb, 'pairing of strain', paste('; gapsize', gapsize, sep='=')), ylab='Size of tracks')
			ttcontsizes = t.test(sizesuptrack ~ suptrackcherries$sameorigin)
			text(x=mean(par('usr')[1:2]), y=max(sizesuptrack)*0.8, cex=1,
			 labels=paste(paste(names(ttcontsizes$estimate), as.character(signif(ttcontsizes$estimate, digits=3)), sep=': ', collapse=','), '; t-test p = ', as.character(signif(ttcontsizes$p.val, digits=3)), ')', sep=''))
			
			# chi-squared test
			nmax = min(sapply(paste(c('same', 'different'), locationtag), function(x){
				sst = sizesuptrack[suptrackcherries$sameorigin==x]
#~ 				min(which(rev(cumsum(rev(table(factor(sst, levels=1:max(sizesuptrack), ordered=T)))))<5))
				min(which(table(factor(sst, levels=1:max(sizesuptrack), ordered=T))<10))-1
			}))
			supnmax = paste('>', nmax, sep='') ; print(supnmax)
			tcontsizes = t(sapply(paste(c('same', 'different'), locationtag), function(x){
				sst = sizesuptrack[suptrackcherries$sameorigin==x]
				table(factor(ifelse(sst<=nmax, as.character(sst), supnmax), levels=c(as.character(1:nmax), supnmax), ordered=T))
			}))
			print(tcontsizes)
			proptcontsizes = tcontsizes / rowSums(tcontsizes)
			at = barplot(proptcontsizes, beside=T, main=paste('Blocks of genes with common', nlPPb, 'pairing of strain', paste('; gapsize', gapsize, sep='=')), xlab='Size of tracks', ylab='Proportion', col=c('white', 'grey'), ylim=c(0, 1))
			legend('topright', fill=c('white', 'grey'), legend=paste(c('same', 'different'), locationtag))
			ct = chisq.test(tcontsizes)
			text(x=mean(at), y=max(proptcontsizes)/2, cex=1,
			 labels=paste('same vs. different ', locationtag, ': Chi-squared (df = ', ct$param, ') = ', as.character(signif(ct$statistic, digits=3)), '; test p = ', as.character(signif(ct$p.val, digits=3)), ')', sep=''))
			text(x=at, y=proptcontsizes+(max(proptcontsizes)*0.05), labels=as.character(signif(ct$residuals, digits=3)), cex=0.5)
			
		}
		return(suptrack)
	})
})
names(congruenttracks) = names(lPPb)
dev.off()

pdf(paste(dircladeprofile, 'top_long_tracks.pdf', sep='/'), width=25, height=21)
nperpage = 7
layout(matrix(1:nperpage, nperpage, 1))
cgt = congruenttracks[['compatible']][[2]]
n = 1
for (bip in unique(cgt[order(cgt$track.size, decreasing=T),'bipart'])[1:21]){
	par(mar=c(5,4,2,2))
	plot(c(0, PPbiparts[bip,], 0), ylim=c(-1,1), type='l', xaxt='n', yaxt='n', main=cgt$smallclade[cgt$bipart==bip][1], xlab='', ylab='Compatibility        Support    ')
	axis(side=2, labels=abs(seq(-1,1,by=0.5)), at=seq(-1,1,by=0.5), las=2)
	polygon(c(1, 1:(dim(PPbiparts)[2]+2), dim(PPbiparts)[2]+2), c(1, 0, PPbiparts[bip,], 0, 1), density=NA, col='lightblue')
	polygon(1:(dim(PPbiparts)[2]+2), c(-1, rejPPbiparts[bip,]-1, -1), density=NA, col='navyblue')
	abline(h=0)
	if (n%%nperpage==0){ axis(side=1, labels=c('', gnames, ''), at=1:(dim(PPbiparts)[2]+2), las=2) }
	n = n + 1
}
dev.off()
