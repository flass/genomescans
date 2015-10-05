#!/usr/bin/Rscript --no-save --no-restore
library('ade4')
library('ape')
library('pegas')
library('gplots')
nbcores = 6

homedir = Sys.getenv()['HOME']
source(paste(homedir, 'scripts/misc/utils-phylo.r', sep='/'))

dark.heat.colors = function(n, k=11){
	ndark = n - k
	darkcol = sapply(1:(ndark), function(i){ rgb((1/ndark)*(i-1), 0, 0,1) })
	print(c(n,k))
	return(c(darkcol, heat.colors(k)))
}

plotCorTestMatrix = function(rvmat, pval='p.value', estr='rv', maxdim.cellnote=30, triangular=T, ...){
	N = dim(rvmat)[1]	
	nrvmat = rvmat
	if (triangular){
		for (i in 1:N){
			for (j in 1:N){
				if (j>=i) nrvmat[i,j]=NA
		}}
	}
	rvpval = log10(sapply(1:N, function(i){ sapply(1:N, function(j){ ifelse((i!=j & !is.na(nrvmat[i,j][[1]][1])), nrvmat[i,j][[1]][[pval]], NA) }) }) * (N*(N-1)/2))
	dimnames(rvpval) = dimnames(nrvmat) 
	colbrk = c(floor(min(rvpval, na.rm=T)), seq(-20, 0, 1), ceiling(max(rvpval, na.rm=T)))
	if (N<=maxdim.cellnote){
		rvest = sapply(1:N, function(i){ sapply(1:N, function(j){ ifelse(!is.na(nrvmat[i,j][[1]][1]), round(nrvmat[i,j][[1]][[estr]], 2), NA) }) })
		dimnames(rvest) = dimnames(nrvmat)
		heatmap.2(rvpval, scale='none', Rowv=NA, Colv=NA, trace='none', cellnote=rvest, notecol='black', na.col='white',
		 breaks=colbrk, cexRow=1.5 , cexCol=1.5, margins = c(12, 12), col='dark.heat.colors', ...)
		# alternate plot with circles the area is proportional to r^2 <=> radius = r/2
		xs = unlist(lapply(1:dim(rvest)[1], function(i){ u = 1:dim(rvest)[1] ; return(u[u>i]-0.5) }))
		ys = unlist(lapply(1:dim(rvest)[1], function(i){ u = 1:dim(rvest)[1] ; return(rep(-i+0.5, length(which(u>i)))) }))
		zs = unlist(lapply(1:dim(rvest)[1], function(i){ u = 1:dim(rvest)[1] ; return(rvest[i,u>i]/4) }))
#~ 		plot(x=xs, y=ys, pch='+')
		symbols(x=xs, y=ys, circles=zs, inches=F, col='turquoise')
	}else{
		heatmap.2(rvpval, scale='none', Rowv=NA, Colv=NA, trace='none', na.col='white', breaks=colbrk, cexRow=1.5, cexCol=1.5, margins = c(12, 12), col='dark.heat.colors', ...) #
	}
	colvec = dark.heat.colors(length(colbrk)-1)
	plot(x=colbrk, y=(0:(length(colbrk)-1))/(length(colbrk)-1), col='white')
	lapply(1:(length(colbrk)-1), function(i, ybt=par()$usr[3:4], ...){
		pair = colbrk[c(i, i+1)]
		rect(xleft=pair[1], xright=pair[2], ybottom=ybt[1], ytop=ybt[2], density=NA, col=colvec[i])
	})
	# use correlation coefficient as a similarity for plotting distances through MDS/PCoA
	print(rvmat)
	rvdist = data.matrix(sapply(1:N, function(i){ sapply(1:N, function(j){ 1-as.numeric(rvmat[i,j][[1]][1]) }) }))
	dimnames(rvdist) = dimnames(rvmat)
	rv.pcoa = pcoa(rvdist)
	print(rv.pcoa$Relative_eig)
	biplot(rv.pcoa)
}

coul = function(n){ c('#FFFFFFFF', rainbow(n-1)) }

plotHiLDSeqTypes = function(genes, dbbiparts, gtdet, ltax, dir.out=NULL, lcds.aln=NULL, minrbl=0.1, minbs=0.75, fixnbcut=7, relative.lengths=TRUE, term.branch.standard=FALSE, size.weights=FALSE, excl.gene=NULL, subsample=NULL, n.subsample=100){	# dump the supported clades in specific genes
	genenames = names(genes)
	if (file.info(dbbiparts)$isdir){
		dirgeno = paste(dbbiparts, '/list_genotypes_rvt', minrbl, '_msw', minbs, '_fnc', fixnbcut, sep='')
		print(dirgeno)
		lseqtypes = lapply(genenames, function(g){
			strsplit(readLines(paste(dirgeno, paste(g, 'geno_labels', sep='.'), sep='/')), split='\t')
		})
		names(lseqtypes) = genenames
		matprofseqtypes = sapply(genenames, function(g){
			print(g)
			if (!is.null(lcds.aln[[g]])){
				# check the presence of all strains in alignment (as it may be eluded in bipart coding)
				excl.strains = setdiff(ltax, rownames(lcds.aln[[g]]))
			}else{ excl.strains = NULL }
			if (!is.null(lcds.aln)){
				codest = codeSequenceTypes(lseqtypes[[g]], ltax, excl.strains=excl.strains, excl.as=0)
				if (length(excl.strains)>0){ print(c(g, 'excl:', excl.strains)) ; print(class(codest)) ; print(codest)}
				return(codest)
			}else{ return(rep(NA, length(ltax))) }
		})
	}else{
		matprofseqtypes = sapply(genenames, function(g){
			ggtdet = gtdet[gtdet$gene_label==g,]
			if (!is.null(lcds.aln)){
				# check the presence of all strains in alignment (as it may be eluded in bipart coding)
				excl.strains = setdiff(ltax, rownames(lcds.aln[[g]]))
			}else{ excl.strains = NULL }
			if (relative.lengths){ branchcol =  ggtdet$rel_branch_length
			}else{ branchcol =  ggtdet$branch_length }
			if (size.weights){ branchweicol = ggtdet$smallerclade_size
			}else{ branchweicol = rep(1, dim(ggtdet)[1]) }
			if (term.branch.standard){ 
				termlens = branchcol[ggtdet$smallerclade_size==1]
				brlenthresh = median(termlens)*minrbl
			}else{ brlenthresh = minrbl }
			ggt = ggtdet[branchcol*branchweicol>brlenthresh & ggtdet$branch_support>minbs, ]
			print(c(g, 'brlenthresh', brlenthresh, 'nlongbranch', dim(ggt)[1]))
			if (!is.null(lcds.aln)){
				if (dim(ggt)[1]>0){
					gbip = dbbiparts$profile[ggt$bipart]
					gsupcla = getSupportedClades(gbip, ltax)
					codest = codeSequenceTypes(gsupcla, ltax, excl.strains=excl.strains, excl.as=0)
					return(codest)
				}else{ return(ifelse(ltax %in% excl.strains, 0, 1)) }
			}else{ return(rep(NA, length(ltax))) }
		})
	}
	genocounts = apply(matprofseqtypes, 2, function(x){nlevels(as.factor(x))})
	barplot(genocounts, names.arg=genes[colnames(matprofseqtypes)], cex.lab=0.5, las=2, ylab='number of distinct genotypes')
	coords = t(simplify2array(strsplit(genenames, split='_|-')))
	colnames(coords) = c('locus_tag', 'b', 'e', 'gene', 'pb', 'pe')
	gene.begin = as.numeric(coords[,'b']) + (as.numeric(coords[,'pb'])-1)
	gene.end = apply(cbind((as.numeric(coords[,'b']) + (as.numeric(coords[,'pe'])-1)), as.numeric(coords[,'e'])), 1, min)
	coordgenocounts = cbind(coords[,c('locus_tag', 'gene')], gene.begin, gene.end, genocounts)
	coordgenocountsext = coordgenocounts[1,]
	for (i in 2:dim(coordgenocounts)[1]){
		intercoord = c(as.numeric(coordgenocounts[i-1,'gene.end'])+1, as.numeric(coordgenocounts[i,'gene.begin'])-1)
		if (intercoord[1] < intercoord[2]){
			interline = c('intergene', 'intergene', intercoord, 0)
			coordgenocountsext = rbind(coordgenocountsext, interline, coordgenocounts[i,])
		}else{
			coordgenocountsext = rbind(coordgenocountsext, coordgenocounts[i,])
		}
	}
	write.table(coordgenocountsext, file=paste(dir.out, '/count_genotypes_rvt', minrbl, '_msw', minbs, '.tab', sep=''), sep='\t', quote=F, row.names=F)
	
	dimnames(matprofseqtypes) = list(ltax, as.vector(genes))
	if (is.logical(dir.out) && !dir.out){ dir.out = NULL }
	if (!is.null(lcds.aln) & !is.null(dir.out)){ write.table(matprofseqtypes, file=paste(dir.out, 'Seq-Type_defining_biparts.tab', sep='/'), sep='\t') }
	# matrix of alleles per gene
	nst = max(matprofseqtypes, na.rm=TRUE)
	heatmap.2(matprofseqtypes, distfun=distcat, Colv=NULL,
	 col='coul', breaks=c(seq(-1, nst)+0.5),
	 trace='none', notecol='black', cexRow=1.5 , cexCol=1.5, margins = c(12, 12),
	 main=paste('Profiles of sequence types defined with\nrelative branch length >', minrbl, 'and branch_support >', minbs)) # , cellnote=matprofseqtypes
	
	
	# drop genes which profile is homogeneous
	mpst = matprofseqtypes[,apply(matprofseqtypes, 2, function(x){ nlevels(as.factor(x)) })>1]
	
	print(sprintf("dropped %d / %d genes which profile is homogeneous", dim(matprofseqtypes)[2] - dim(mpst)[2], dim(matprofseqtypes)[2]))
	
	mpst.coa = dudi.coa(mpst, scannf=F, nf=2)
	print(mpst.coa$eig)
	scatter(mpst.coa)
	
	nfsav = paste(dir.out, 'Seq-Type_correlations.RData', sep='/')
	if (file.exists(nfsav)){ load(nfsav)
	}else{
		# matrix correllation (RV cefficient) of alleles occurence among gene (pair-wise gene comparison)
		trv = testCorSeqTypes(mpst)
		print(trv)
		if (!is.null(lcds.aln) & !is.null(dir.out)){ save(trv, matprofseqtypes, file=nfsav) }
	}
	plotCorTestMatrix(trv, main=paste('Pairwise RV coefficient (correlation) tests\nbetween profiles of sequence types defined with\nrelative branch length >', minrbl, 'and branch_support >', minbs))
	plotCorTestMatrix(trv[1:30,1:30], main=paste('Pairwise RV coefficient (correlation) tests\nbetween profiles of sequence types defined with\nrelative branch length >', minrbl, 'and branch_support >', minbs))
	plotCorTestMatrix(trv[88:104,88:104], main=paste('Pairwise RV coefficient (correlation) tests\nbetween profiles of sequence types defined with\nrelative branch length >', minrbl, 'and branch_support >', minbs))
	return(matprofseqtypes)
}
	




carg = commandArgs(trailingOnly=TRUE)
dircladeprofile = carg[1]
cdsalndir = carg[2]
nfgenenames = carg[3]
dirout = carg[4]
# bipart description
if (length(carg)>4){ 
	print(paste('read genotype classifications in', dircladeprofile))
	dbbiparts = carg[5]
	gtdet = NULL
}else{
	stop("deprecated use!")
	print(paste('read bipartition tables in', dircladeprofile))
	dbbiparts = read.table(paste(dircladeprofile, 'bipart_db.tab', sep='/'), h=T, colClasses=c('character', rep('numeric', 4)))
	gtdet = read.table(paste(dircladeprofile, 'bipart_intrees.tab', sep='/'), h=T)
}
radout = strsplit(nfgenenames, split='/')[[1]] ; radout = strsplit(radout[length(radout)], split='\\.')[[1]][1]
diroutseqtypes = paste(dirout, paste(radout, 'SeqType_alns', sep='-'), sep='/')
# taxon list
ltax = scan(paste(dircladeprofile, 'taxlabels', sep='/'), what='character')
# gene list
genenames = readLines(nfgenenames)
# get order by coordinates
genecoords = as.data.frame(t(sapply(genenames, function(x){ as.numeric(strsplit(strsplit(x, split='_')[[1]][2], split='-')[[1]]) })))
genecoordorder = order(genecoords[,1])
genenames = genenames[genecoordorder]
# make short gene names
gnames = sapply(genenames, function(x){ paste(strsplit(x, split='_')[[1]][3:4], collapse='_') })
genes = gnames[genenames]
names(genes) = genenames
names(genenames) = genes


# get the gene alignments
lnfaln = list.files(cdsalndir)
lnfaln = lnfaln[sapply(lnfaln, function(x){ strsplit(x, split='\\.')[[1]][1] }) %in% genenames]
lalnpath = paste(cdsalndir, lnfaln, sep='/')
names(lalnpath) = lnfaln
lcds.aln = lapply(lnfaln, function(nfaln){
	cds.aln = read.dna(lalnpath[nfaln], format='fasta')
	keep = sapply(rownames(cds.aln), function(rn){ sum(base.freq(cds.aln[rn,], all=T, freq=T)[c('a', 'c', 'g', 't')])>30 })
	cds.aln[keep,]
})
names(lcds.aln) = sapply(strsplit(lnfaln, split='\\.'), `[`, 1)


pdf(paste(dirout, paste(radout, 'SeqType-profile.pdf', sep='-'), sep='/'), width=20, height=20)
llmatprofseqtypes = lapply(c(2), function(minrbl){
	lapply(c(5), function(minbs){
		print(c(minrbl, minbs))
		mainfig = FALSE
#~ 		mainfig = minrbl==3 & minbs==0.9
#~ 		mainfig = minrbl==0.02 & minbs==0.9
#~ 		mainfig = minrbl==0.15 & minbs==0.9
		matprofseqtypes = plotHiLDSeqTypes(genes, dbbiparts, gtdet, ltax, lcds.aln=lcds.aln, minrbl=minrbl, relative.lengths=F, term.branch.standard=T, dir.out=diroutseqtypes, minbs=minbs)
		if (mainfig){
			dir.create(diroutseqtypes, showWarning=F)
			lcds.stconsaln = exportSeqTypesAlignments(matprofseqtypes, lcds.aln, dir.out=diroutseqtypes) #, gname.long.alias=genenames
			lcds.stcons = lcds.stconsaln[1,]
			lcds.staln = lcds.stconsaln[2,]
			pdf(paste(dirout, paste(radout, 'DistConsensus_SeqTypes.pdf', sep='-'), sep='/'), width=20, height=20)
			models = c('nucleotide divergence', 'distances under a F81 substitution model') ; names(models) = c('raw', 'F81')
			for (mod in names(models)){
				lcs.dist.cons = lapply(lcds.stcons, dist.dna, model=mod)
				bxp = boxplot(lcs.dist.cons, names=colnames(matprofseqtypes),
				 main=paste('Distribution of', models[mod], '\nbetween consensus of sequence types'), mar=c(20, 4, 6, 2)+0.1, las=2, cex.main=2)
				points(x=do.call(c, lapply(1:length(lcs.dist.cons), function(i){ rep(i, length(lcs.dist.cons[[i]])) })), y=do.call(c, lcs.dist.cons))
	#~ 			mtext(at=0:length(lcs.dist.cons), text=c('n =',  sapply(lcds.stcons, dim)[2,]), line=18)
				text(x=0:length(lcs.dist.cons), y=0, labels=c('# ST\n# sites',  sapply(lcds.stcons, function(x){paste(dim(x), collapse='\n')})))
				points(x=1:length(lcs.dist.cons), y=sapply(lcs.dist.cons, mean, na.rm=T), pch=3, col='red')
				legend('topright', pch=3, col='red', legend='average', cex=2)
			}
			dev.off()
			pdf(paste(dirout, paste(radout, 'FoldedSFS_SeqTypes.pdf', sep='-'), sep='/'), width=20, height=20)
			for (ge in genenames){
				lstaln = lcds.staln[[ge]]
				several = sapply(lcds.staln[[ge]], dim)[1,] >= 4
				nbst = length(which(several))
				layout(matrix(c(1:nbst, rep(0, 9-nbst)),3,3, byrow=T))
				for (st in names(lstaln)[several]){
					print(st)
					staln = lstaln[[st]]
					print(dim(staln))
					ali = 1:dim(staln)[2]
					sfs123 = t(sapply(1:3, function(pos){
						ipos = ali[ali%%3==(pos-1)]
						site.spectrum(staln[,ipos], folded=T)
					}))
					nby2 = dim(sfs123)[2]
					bp = barplot(sfs123, beside=T, col=rep(c('orange', 'red', 'gold'), nby2), names.arg=1:nby2, cex.names=2, cex.axis=2, cex.main=2, main=paste(ge, st))
					legend('topright', fill=c('orange', 'red', 'gold'), legend=paste('codon pos', 1:3), cex=2)
					nucdiv = nuc.div(staln, pairwise.deletion=T)
					print(nucdiv)
					tajDt = tajima.test(staln)
					tajDt = lapply(tajDt, signif, digits=3)
					print(tajDt)
					text(x=mean(bp), y=max(sfs123)*0.9, labels=paste('Pi =', signif(nucdiv, digits=3)), cex=2)
					text(x=mean(bp), y=max(sfs123)*0.8, labels=paste('Tajima\'s D =', tajDt$D, '; p(Norm) =', tajDt$Pval.norm, '; p(Beta) =', tajDt$Pval.beta), cex=2)
					warnings()
				}
			}
			dev.off()
		}
	})
})

dev.off()
