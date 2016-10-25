#!/usr/bin/R
#~ library('seqinr')
library('ape')
#~ library('pegas')
library('parallel')
#~ library('RColorBrewer')
library('phangorn')

### Misc
printProgressUpperMatrix = function(i, N, step=1000, initclock=FALSE){
	if ( i%%round(N/step)==0 | N<((step/2)+1)){
		s = paste("\r", round(100*(i/N)*(2*N-i-1)/(N-1), digits=1), "%\t")
		if (initclock) s = paste(s, format(Sys.time() - initclock, format="%H:%M:%S"), "\t")
		cat(s)
	}
}


### Distance matrix-related
read.phylipDist = function(nf, asmatrix=FALSE){
	MLdistlines = sapply(readLines(nf), strsplit, split=' +')
	if (length(MLdistlines)<1){return(NULL)}
	matsize = as.numeric(MLdistlines[[1]][2])
	strains = sapply(2:(matsize+1), function(i){MLdistlines[[i]][1]})
	distmat = t(sapply(2:(matsize+1), function(i){
		if (i<=2){ d = numeric() }else{ d = as.numeric(MLdistlines[[i]][2:(i-1)]) }
		c(d, rep(NA, matsize-i+2)) 
	}))
	dimnames(distmat) = list(strains, strains)
	if (asmatrix){ return(distmat) 
	}else{ return(as.dist(distmat)) }
}

# !!! always with i < j <= n
collowmat = function(i, n){ n*(i-1) - i*(i-1)/2 }
getdistmatindex = function(i, j, n){ collowmat(i, n) + j-i }
getdistmatcoords = function(k,n){
	i = 1
	while ( (collowmat(i+1, n) < k) & (i < n-1)){ i = i + 1 }
	j = k + i - collowmat(i, n)
	return(c(i, j))
}

getdistmatindexesfromlabel = function(lab, distmat=NULL, n=NULL, labs=NULL){
	if (!is.null(distmat)){
		labs = attr(distmat, 'Labels')
		n = attr(distmat, 'Size')
	}else{ stopifnot(!is.null(n), !is.null(labs)) }
	k = which(labs==lab)
	if (k>1){	idi = sapply(1:(k-1), function(i){ getdistmatindex(i, k, n) }) 
	}else{ idi = numeric() }
	if (k<n){	idj = sapply((k+1):n, function(j){ getdistmatindex(k, j, n) }) 
	}else{ idj = numeric() }
	return(c(idi, idj))
}

getlabelsfromdistmatindexes = function(k, distmat=NULL, n=NULL, labs=NULL){
	if (!is.null(distmat)){
		labs = attr(distmat, 'Labels')
		n = attr(distmat, 'Size')
	}else{ stopifnot(!is.null(n), !is.null(labs)) }
	ij = sapply(k, getdistmatcoords, n=n)
	return(t(apply(ij, 2, function(x){ labs[x] })))
}

singleLinkageClustersFromMatrix = function(mat, measure='dist', dist.thresh=0.05){
	# gets a distance matrix
	if (class(mat)!='dist')
	if (type=='sim'){ mat = 1 - mat 	
	}else{ if (type!='sim') stop("needs matrix of class 'dist', containing similarity or distance measures (measure='sim' or measure='dist')") }
	clusters = list()
	labs = attr(mat, 'Labels')
	N = attr(mat, 'Size')
	for (i in 1:(N-1)){
		for (j in (i+1):N){
			z = getdistmatindex(i, j, N)
			# print(paste(i, 'x', j, '=>', z, ', dist', mat[z]), quote=F)
			stopifnot(!is.null(mat[z]))
			if ((mat[z] < dist.thresh) & !is.na(mat[z])){
				clust = labs[c(i, j)]
				# print(paste('match', i, j, paste(clust, collapse= ', ')), quote=F)
				K = length(clusters)
				hits = numeric()
				if (K>0){
					for (k in 1:K){
						for (t in 1:2){
							if (clust[t] %in% clusters[[k]]){
								# print(paste(clust[t], '    hit:', k, ':', paste(clusters[[k]], collapse= ', ')), quote=F)
								clusters[[k]] = union(clusters[[k]], clust)
								# print(paste('    new:', paste(clusters[[k]], collapse= ', ')), quote=F)
								hits = c(hits, k)
							}
						}
					}
				}
				hits = unique(hits)
				if (length(hits)==0){ clusters = c(clusters, list(clust)) }
				if (length(hits)>1){ 
					# print(paste('    hits:', paste(hits, collapse= ', ')), quote=F)
					# merge clusters
					ref = hits[1]
					incl = hits[2:length(hits)]
					for (p in incl){
						# print(paste('    merge:', ref, '+', p, ':', paste(clusters[[ref]], collapse= ', '), paste(clusters[[p]], collapse= ', ')), quote=F)
						clusters[[ref]] = union(clusters[[ref]], clusters[[p]])
						# print(paste('    new:', paste(clusters[[ref]], collapse= ', ')), quote=F)
					}
					# remove included clusters
					# print(paste('remove redundant clusters', paste(incl, collapse= ', ')), quote=F)
					clusters = clusters[-incl]
				}
			}
			# print('clusters', quote=F)
			# print(clusters, quote=F)
		}
	}
	clusters = lapply(unique(clusters), sort)
	unclustered = setdiff(labs, c(clusters, recursive=TRUE))
	# exclude those that have mostly missing data (> 50% absent comparisons in dist matrix)
	nazes = character()
	for (lab in unclustered){
		is = getdistmatindexesfromlabel(lab, n=N, labs=labs)
		nas = which(is.na(mat[is]))
		if (length(nas)/(N-1) > 0.5){ nazes = c(nazes, lab) }
	}
	unclustered = setdiff(unclustered,  nazes)
	# add singletons
	clusters = c(clusters, unclustered)
	tagged = sapply(labs, function(lab){ 
		whereis = sapply(clusters, function(x){ lab %in% x })
		ifelse(any(whereis), which(whereis), NA)
	})
	return(tagged)
}

getIndexesOutUpperTriangular = function(i,j, max.dist=NULL){ return(j<=i) }
getIndexesOutUpperRibbon = function(i,j, max.dist){ if (j<=i){ return(TRUE) }else{ return(j > max.dist+i) } }
getIndexesOutDiagonalRibbon = function(i,j, max.dist){ return((j < i-max.dist) | (j > i+max.dist)) }
getIndexesUpperRibbon = function(i,j, max.dist){ return(j>i & j<=max.dist+i) }
getIndexesLowerRibbon = function(i,j, max.dist){ return(j<i & j>=i-max.dist) }

# retrieve indexes of a NxN square matrix using a conditional function ; return them as a Nx2 matrix (directly usable as an index array)
conditionalMatrixIndexes = function(matdim, fun, multiproc=1, ...){
	starttime = Sys.time()
	condi = function(i){
			condok = which(sapply(1:matdim, fun, i=i, ...))
			printProgressUpperMatrix(i, matdim, step=100, initclock=starttime)
			return(cbind(rep(i, length(condok)), condok))
	}
	if (multiproc > 1){
		do.call(rbind, mclapply(1:matdim, condi, mc.cores=multiproc, mc.preschedule=FALSE))
	}else{
		do.call(rbind, lapply(1:matdim, condi))
	}
}

linkageDisequilibrium = function(aln, metric="r", discard.gaps=TRUE, multiproc=1, verbose=FALSE, quiet=FALSE, max.dist=NULL, mem.light=FALSE, upper.triangular=TRUE, full.matrix=FALSE, gapchar=c('-')){
	# assumes that the input alignment contains only bi-allelic sites ; other kind of variation will infringe assumptions for computations LD metrics, and likely yield abberant results (e.g. off [0;1] interval)
	# NB: gaps are natively treated as a distinct character ; if 'discard.gaps' is set to TRUE, gap-containg rows of the alignment at both considered sites will be discarded and computations will be done on the remaining rows.
	# return a numeric.matrix by default, but can rather return a list if 'mem.light' is set to TRUE, where every list element cover a row of the matrix, discarding all cells that were not computed because off the specified range (and would come as NA values in the matrix output)
	retD = (metric=="D")
	retDp = (metric=="D'")
	retr = (metric=="r")
	retr2 = (metric=="r2")
	retrho = (metric=="rho")
	retfi = (metric=="Fisher")
	retcompat = (metric=="compat")
	if (!any(retD, retDp, retr, retr2, retrho, retfi, retcompat)){ stop(paste('a wrong metric "', metric, '" was provided ; correct metrics are "D", "D\'", "r" "r2", "rho", "Fisher" or "compat"', sep='')) }
	M = dim(aln)[1]
	N = dim(aln)[2]
	if (full.matrix){
		if (is.null(max.dist)){ testjout = function(i,j, max.dist=NULL){ return(FALSE) }
		}else{ testjout = getIndexesOutDiagonalRibbon }	
	}else{
		if (is.null(max.dist)){ testjout = getIndexesOutUpperTriangular
		}else{ testjout = getIndexesOutUpperRibbon }
	}
	starttime = Sys.time()
	oneVsAllSiteLD = function(i){
		if (!quiet){ printProgressUpperMatrix(i, N, step=100, initclock=starttime) }
		if (discard.gaps){
			igaps = aln[,i] %in% gapchar
		}else{
			m = M
			si = aln[,i]
			ti = table(si)
			stopifnot(length(ti==2) & sum(ti)==m)	# asserts that there is 2 and only 2 alleles, which frequencies sum to the total
			ps = ti/m
			p1 = ps[1]
			p2 = ps[2]
			A1 = names(p1)
		}
		if (mem.light){ pair.range = which(!sapply(1:N, testjout, i=i, max.dist=max.dist)) 
		}else{pair.range = 1:N }
		sitecomp = sapply(pair.range, function(j){
			if (!mem.light && testjout(i,j, max.dist=max.dist)){ return(NA) }
			else{
				if (discard.gaps){
					jgaps = aln[,j] %in% gapchar
					ijnogaps = (!igaps & !jgaps)
					m = length(which(ijnogaps))
					if (m==0){ return(NA) }
#~ 					print(paste('m', m))
					si = aln[ijnogaps,i]
					ti = table(si)
#~ 					print(paste('ti', ti))
					stopifnot(length(ti==2) & sum(ti)==m)	# asserts that there is 2 and only 2 alleles, which frequencies sum to the total
					ps = ti/m
					p1 = ps[1]
					p2 = ps[2]
					A1 = names(p1)
					sj = aln[ijnogaps,j]
				}else{
					sj = aln[,j]
				}
				tj = table(sj)
				stopifnot(length(tj==2) & sum(tj)==m)	# asserts that there is 2 and only 2 alleles, which frequencies sum to the total
				qs = tj/m
				q1 = qs[1]
				q2 = qs[2]
				B1 = names(q1)
				ABs = apply(expand.grid(names(ps), names(qs)),1,paste, collapse=' ')
				h = factor(paste(si,sj), levels=ABs)
				xs = table(h)/m
				A1B1 = paste(A1, B1)
				x11 = xs[A1B1]
				if (retcompat){
					# four-gamete test: compatibility <=> maximum three of the four bi-allelic states are observed
					# return boolean stating sites are compatible (1) or not (0)
#~ 					return(as.integrer( 0 %in% xs ))
					return( 0 %in% xs ) # logical format, more compact
				}else{ if (retfi){
					tcont = matrix(table(h), 2, 2)
#~ 					print(tcont)
					Fpval = fisher.test(tcont)$p.value
					if (Fpval > 1+1e-12){
						print(tcont)
						print(fisher.test(tcont))
						stop('Error: p-value > 1')
					}
#~ 					print(Fpval)
					return(Fpval)
				}else{
				D = x11 - p1*q1
				if (verbose){ 
					print('', quote=F)
					print(c(i,j), quote=F)
					print(ps)
					print(qs)
					print(xs)
					print(c('D', D))
				}
				if (retD){ 
					return(D) 
				}else{ if (retr | retr2){
					r = D/sqrt(p1*q1*p2*q2)
					if (verbose){ 
						print(c('r', r))
					}
					if (retr){ return(r)
					}else{ return(r^2) }
				}else{ if (retDp){
					if (D < 0){ Dp = abs(D)/min(p1*q1, p2*q2)
					}else{ Dp = abs(D)/min(p1*q2, p2*q1) }
					return(Dp)
				}else{ if (retrho){
					rho = D/(p1*q2)
					return(rho)
				}else{ stop("undetermined return value")	
				}}}}}}
			}
		})
#~ 		if (is.list(sitecomp)){
#~ 			# inconsistent element lengths (some elemen have zero-length), correct it
#~ 			for (k in which(!as.logical(sapply(u, length)))){ sitecomp[[k]] = NA }
#~ 			sitecomp = unlist(sitecomp)
#~ 		}
#~ 		return(sitecomp)
	}
	
	if (multiproc > 1){
#~ 		lsiD = mclapply(1:N, oneVsAllSiteLD, mc.cores=multiproc, mc.preschedule=FALSE)
#~ 		msijD = do.call('rbind', lsiD)
		msijD = mclapply(1:N, oneVsAllSiteLD, mc.cores=multiproc, mc.preschedule=FALSE)
	}else{
		msijD = lapply(1:N, oneVsAllSiteLD)
	}
	if (!mem.light){ 
		msijD = simplify2array(msijD) 
		if (upper.triangular){ msijD = t(msijD) }
	}
	if (!quiet){ cat("\n") } 
	return(msijD)
}

as.matrix.LDlist = function(lsijD, quiet=TRUE){
	# return a lower triangular matrix
	starttime = Sys.time()
	N = length(lsijD)
	m = sapply(1:N, function(i){
		if (!quiet){ printProgressUpperMatrix(i, N, step=1000, initclock=starttime)  }
		c(rep(NA, i), as.numeric(lsijD[[i]]), rep(NA, N-i-length(lsijD[[i]])))
	})
	if (!quiet){ cat("\n") } 
	return(m)
}

unfold.triangular.LDlist = function(lsijD, maxLD=1){
	for (i in 1:length(lsijD)){
		js = sapply(1:i, function(j){ lsijD[[j]][i-j] })
		lsijD[[i]] = c(js, maxLD, lsijD[[i]])
	}
	return(NULL)
}

exponentialDecay = function(ldmat, pos, max.dist=1000, multiproc=1, quiet=FALSE, format='upper'){
	starttime = Sys.time()
	if ('list' %in% strsplit(format, '\\.')[[1]]){ N = length(ldmat)
	}else{ N = dim(ldmat)[1] }
	if (is.null(max.dist)){ k = Inf }else{ k = max.dist }
	if (format=='triangular.list'){ unfold.triangular.LDlist(ldmat, maxLD=1) }
#~ 	simplify2array(mclapply(1:N, function(i){
	sapply(1:N, function(i){
		if (!quiet){ printProgressUpperMatrix(i, N, step=100, initclock=starttime) }
		if ('list' %in% strsplit(format, '\\.')[[1]]){
			ldi = ldmat[[i]]
			pwld = ldi[max(1, floor((length(ldi)-k)/2)):min(length(ldi), floor((length(ldi)+k)/2))]
		}else{
			if (format=='upper'){
				# assumes the matrix to be upper triangular, as outputed by linkageDisequilibrium()
				if (i>1) v = ldmat[max(i-k,1):(i-1),i] else v = numeric()
				if (i<N) w = ldmat[i,(i+1):min(i+k,N)] else w = numeric()
			}else{ if (format=='lower'){
				# assumes the matrix to be lower triangular, as outputed by linkageDisequilibrium()
				if (i>1) v = ldmat[max(i-k,1):(i-1),i] else v = numeric()
				if (i<N) w = ldmat[i,(i+1):min(i+k,N)] else w = numeric()
			}else{stop('wrong format')}}
			pwld = c(v[!is.na(v)], 1, w[!is.na(w)])
	#~ 		pwld = c(v[!is.na(v)], w[!is.na(w)])
		}
		# mesure decay of a negative exponential fit
		di = abs(pos[i] - pos[(max(i-k,1):min(i+k,N))])
#~ 		di = abs(pos[i] - pos[setdiff((max(i-k,1):min(i+k,N)), i)])
		logdi = ifelse(di>0, log(di), NA)
		logld = ifelse(pwld>0, log(pwld), NA)
		print(c(length(logdi), length(logld)))
		lmlld = lm(logld ~ logdi)
		rate = lmlld$coefficients['logdi']	
#~ 	}, mc.cores=multiproc, mc.preschedule=FALSE))
	})
}

sampleLabels = function(index, nlabs=100){ ifelse(1:length(index) %in% round(seq(1, length(index), length.out=nlabs)), index, NA) }

plotLDHeatMaps = function(ldmat, full.aln.index, nlabs=100, main="correlation coefficient r^2 between sites", ignore.nonmat=F){
	labsample = sampleLabels(full.aln.index, nlabs=nlabs)
	if (class(ldmat)=='matrix'){
		gplots::heatmap.2(data.matrix(ldmat), Rowv=FALSE, Colv='Rowv', dendrogram="none", col=RColorBrewer::brewer.pal(9, 'Blues'), trace="none", cexRow=1.5 , cexCol=1.5, margins = c(12, 12),
		 labRow=labsample, labCol=labsample, main=main)
	}else{  if (!ignore.nonmat) stop("'ldmat' must be a numeric matrix") }
	return(NULL)
}

absentSeq = function(aln, propgap=1.0, relpropgap=0){
	alnpg = base.freq(aln, all=T)['-']
	sapply(1:dim(aln)[1], function(i){
		pg = base.freq(aln[i,], all=T)['-']
		return((pg >= propgap) & (pg/alnpg >= relpropgap))
	})
}

gapDensity = function(aln, freq=FALSE){ base.freq(aln, freq=freq, all=T)['-'] }

getMaxNbAlleleIndexes = function(aln, maxnall=1, as.logical=FALSE, multiproc=1){
	sites = simplify2array(mclapply(1:dim(aln)[2], function(j){
		statecounts = table(base.freq(aln[,j], freq=T, all=T))
		sum(statecounts > 0) <= maxnall
	}, mc.cores=multiproc, mc.preschedule=TRUE))
	if (as.logical){ return(sites)
	}else{ return(which(sites)) }
}

getGaplessIndexes = function(aln, as.logical=FALSE, maxgap=0, gapchar='-', multiproc=1){
	gaps = simplify2array(mclapply(1:dim(aln)[2], function(j){
		sum(base.freq(aln[,j], all=T, freq=T)[gapchar], na.rm=T) > maxgap
	}, mc.cores=multiproc, mc.preschedule=TRUE))
	if (as.logical){ return(!gaps)
	}else{ return(which(!gaps)) }
}

getSNPindexes = function(aln, as.logical=FALSE, multiproc=1){
	monomorphic = simplify2array(mclapply(1:dim(aln)[2], function(j){
		1 %in% base.freq(aln[,j])
	}, mc.cores=multiproc, mc.preschedule=TRUE))
	if (as.logical){ return(!monomorphic)
	}else{ return(which(!monomorphic)) }
}

getParsInformIndexes = function(aln, as.logical=FALSE, consider.gaps.as.variant=F, multiproc=1){
	nonsinglepolymorphic = simplify2array(mclapply(1:dim(aln)[2], function(j){
		statecounts = table(base.freq(aln[,j], freq=T, all=consider.gaps.as.variant))
		return(sum(statecounts[names(statecounts)>1])>1)
	}, mc.cores=multiproc, mc.preschedule=TRUE))
	if (as.logical){ return(nonsinglepolymorphic)
	}else{ return(which(nonsinglepolymorphic)) }
}

getMinorAlleleFreq = function(aln, consider.gaps.as.variant=FALSE, multiproc=1){ 
	mincount = simplify2array(mclapply(1:dim(aln)[2], function(j){
		statecounts = table(base.freq(aln[,j], freq=T, all=consider.gaps.as.variant))
		# first term captures biallelic site with sufficient minor allele frequency
		return(min(statecounts))
	}, mc.cores=multiproc, mc.preschedule=TRUE))
	return(mincount)
}

getBiAllelicIndexes = function(aln, as.logical=FALSE, minallelefreq=NULL, nonsingle=FALSE, consider.gaps.as.variant=FALSE, multiproc=1){
	if (nonsingle==TRUE | (!is.null(minallelefreq) && minallelefreq>1)){ print('not considering singletons/rare alleles') }
	if (is.null(minallelefreq)){ minoccur = as.integer(nonsingle) + 1 
	}else{ minoccur = minallelefreq }
	biallelic = simplify2array(mclapply(1:dim(aln)[2], function(j){
		statecounts = table(base.freq(aln[,j], freq=T, all=consider.gaps.as.variant))
		# first term captures biallelic site with sufficient minor allele frequency
		bials = sum(statecounts[names(statecounts)>=minoccur])==2
		# second make sure to filter sites where singletons/rare alleles would make third or further alternative alleles, what would mess up LD calculations!
		nothirdalt = sum(statecounts[names(statecounts)<minoccur & names(statecounts)>0])==0
#~ 		print(c(j, statecounts, nothirdalt))
		return(bials & nothirdalt)
	}, mc.cores=multiproc, mc.preschedule=TRUE))
	if (as.logical){ return(biallelic)
	}else{ return(which(biallelic)) }
}

getGapTracks = function(aln, gapchar='-'){
	gaptracks = lapply(rownames(aln), function(s){
		cseq = as.vector(as.character(aln[s,]))
		gt = data.frame(begin=integer(), end=integer())
		g = c()
		gap = F
		for (i in 1:length(cseq)){
			if (cseq[i]==gapchar){
				if (!gap){	
					# open gap
					g = i
					gap = T
				}
			}else{
				if (gap){
					# close gap
					gt = rbind(gt, c(g, i-1))
					gap = F
				}
			}
		}
		colnames(gt) = c('begin', 'end')
		gt$size = gt$end - gt$begin +1
		return(gt)
	})
	names(gaptracks) = rownames(aln)
	return(gaptracks)
}

toIndelPolymorphismAlignment = function(aln, code=c(0,1), collapse=TRUE){
	charaln = as.character(aln)
	codealn = sapply(1:dim(charaln)[2], function(j){ ifelse(charaln[,j]=='-', code[1], code[2]) })
	if (collapse){
		endaln = codealn[,1]
		m = n = 1
#~ 		ltracks = list() 
		mtracks = matrix(numeric(), 2, 0)
		while (n<=dim(charaln)[2]){
			while (all(codealn[,m]==codealn[,n]) & n<dim(codealn)[2] ){ n = n + 1 }
			if (n<dim(codealn)[2] | any(codealn[,m]!=codealn[,n])){
				# add column describing new homogeneous segment
				endaln = cbind(endaln, codealn[,n])
			}
			# store previous segment coordinates
#~ 			ltracks = c(ltracks, list(c(m, n-1)))
			mtracks = cbind(mtracks, c(m, n-1))
			# update counters to new segment coordinates
			m = n
			n = n + 1 
		}
	}else{
		endaln = codealn
#~ 		ltracks = lapply(1:dim(codealn)[2], function(j){ c(j,j) })
		mtracks = sapply(1:dim(codealn)[2], function(j){ c(j,j) })
	}
	result = list(endaln, mtracks)
	names(result) = c("indel.aln", "indel.coords")
	return(result)
}

convertTo01genotypes = function(aln, get.biallelic=FALSE, get.gapless=FALSE, multiproc=1){
	if (get.biallelic){ bial = getBiAllelicIndexes(aln, as.logical=TRUE) }else{ bial = TRUE }
	if (get.gapless){ nogap = getGaplessIndexes(aln, as.logical=TRUE) }else{ nogap = TRUE }
	caln = as.character(aln)[,(bial & nogap)]
	print(dim(caln))
	galn = simplify2array(mclapply(1:dim(caln)[2], function(i){
#~ 	galn = sapply(1:dim(caln)[2], function(i){
		si = caln[,i]
		ps = table(si)
		if (('-' %in% names(ps)) | (length(ps)!=2)){ stop("can only convert to 0/1 haploid genotypes on gapless biallelic site data")}
		ch = names(ps)
		major = ch[ps==max(ps)][1]
		minor = ch[ch!=major]
		tonum = c(0,1) ; names(tonum) = c(major, minor)
		return(as.vector(tonum[si]))
	}, mc.cores=multiproc, mc.preschedule=TRUE))
#~ 	})
	rownames(galn) = rownames(aln)
	return(galn)
}

exportToPED = function(aln, nfout, num.alleles=FALSE, sep=' ', get.biallelic=TRUE, get.gapless=TRUE, STRUCTUREfmt=FALSE, multiproc=1){
	tonum = 0:4 ; names(tonum) = c("N", "A", "C", "G", "T")
	stopifnot(is.matrix(aln))
	if (get.biallelic){ bial = getBiAllelicIndexes(aln, as.logical=TRUE) }else{ bial = TRUE }
	if (get.gapless){ nogap = getGaplessIndexes(aln, as.logical=TRUE) }else{ nogap = TRUE }
	caln = as.character(aln)
	if (is.vector(caln)){ caln = matrix(caln, dim(aln)[1], dim(aln)[2]) }
	caln = caln[,(bial & nogap)]
	rownames(caln) = rownames(aln)
	L = dim(caln)[2]+6
	if (num.alleles){ nogenotype = "0"
	}else{ nogenotype = "N" }
	
	fout = file(nfout, 'w')
	for (rn in rownames(caln)){
		if (num.alleles){ alleles = tonum[toupper(caln[rn,])]
		}else{ alleles = toupper(caln[rn,]) }
		if (STRUCTUREfmt){  write(c(rn, paste(alleles, alleles)), file=fout, sep=sep, ncolumns=L, append=TRUE)
		}else{ write(c(rn, rn, rep(0, 4), paste(alleles, alleles)), file=fout, sep=sep, ncolumns=L, append=TRUE) }
	}
	close(fout)
#~ 	pedtab = t(simplify2array(mclapply(rownames(caln) ,function(rn){
#~ 		if (num.alleles){ alleles = tonum[toupper(caln[rn,])]
#~ 		}else{ alleles = toupper(caln[rn,]) }
#~ 		if (STRUCTUREfmt){ return(c(rn, paste(alleles, alleles))) 
#~ 		}else{ return(c(rn, rn, rep(0, 4), paste(alleles, alleles))) }
#~ 	}, mc.cores=multiproc, mc.preschedule=TRUE)))
#~ 	write.table(pedtab, file=nfout, sep=' ', row.names=F, col.names=F, quote=F)
}

exportSNPpos = function(indexes, nfout, genpos=rep(0.0, length(indexes)), chr=rep(1, length(indexes)), EIGEINSTRATfmt=FALSE){
	if (EIGEINSTRATfmt){
		write.table(as.data.frame(t(sapply(1:length(indexes), function(k){ c(paste('s', k, sep=''), chr[k], genpos[k], indexes[k]) }))), file=nfout, sep='\t', row.names=F, col.names=F, quote=F)
	}else{
		write.table(as.data.frame(t(sapply(1:length(indexes), function(k){ c(chr[k], paste('s', k, sep=''), genpos[k], indexes[k]) }))), file=nfout, sep='\t', row.names=F, col.names=F, quote=F)
	}
}

haplotypeDiversity = function(aln){
	# expects a character matrix
	if (!is.character(aln)){ aln = as.character(aln) }
	N = dim(aln)[1]
	haplos = apply(aln, 1, paste, sep='', collapse='')
	fhaplos = table(haplos)/N
	return((1 - sum(fhaplos^2))*N/(N-1))
}

fixwidthWindowSubsampleSites = function(alnrange, siteindexes, maxsize, distribute=T){
	there = siteindexes %in% alnrange
	k = min(length(which(there)), maxsize)
#~ 	print('k, alnrange')
#~ 	print(list(k, summary(alnrange)))
	if (k>0){
		if (k==1 | !distribute){ return(sort(sample(siteindexes[there], k)))
		}else{
			# try to maximize the equal spacing of points in space
			refpoints = quantile(alnrange, p=(0:(k-1))/(k-1))
#~ 			print('refpoints')
#~ 			print(refpoints)
			# sample the sites that minimize the distance with quantiles of distance
			pool = siteindexes[there]
#~ 			print('pool')
#~ 			print(pool)
			gridsample = numeric()
#~ 			print('gridsample')
			for (ref in refpoints){
				dist2ref = abs(pool - ref)
				i = which(dist2ref==min(dist2ref))[1]
				gridsample = c(gridsample, pool[i])
				pool = pool[-i]
#~ 				print(gridsample)
			}
			return(sort(gridsample))
		}
	}else{
		return(numeric(0))
	}
}

callPhi = function(aln, alnname, phiwindow=10, phipackpath="~/Programs/PhiPack/Phi", tempdirpath="~/tmp"){ 
	nfaln = sprintf("%s/%s", tempdirpath, alnname)
	write.dna(aln, file=nfaln, format='fasta')
	phicall = sprintf("%s -w %d -f %s -v", phipackpath, phiwindow, nfaln)
	phiout = system(phicall, intern=TRUE)
	phistat =as.numeric(substring(phiout[length(phiout)-9], 15, 25))
	phipval = as.numeric(substring(phiout[length(phiout)-1], 15))
	unlink(nfaln)	# destroy temporary alignment file
	return(c(phistat, phipval))
}

reportsnpdens = function(alnrange){ length(alnrange) }

rollStats = function(aln, windowsize=1000, step=250, foci=NULL, fun=NULL, measures=c("nucdiv", "gapdens"), fun.usedist=list(), fun.userange=list(), fun.res=list(), dist.model='raw', rmAbsSeq=FALSE, propgap=0.35, relpropgap=10, subsample=NULL, tempdirpath="~/tmp", alnname='aln', ctrl.subalnsize=FALSE, multiproc=1, quiet=FALSE){	#, measures=c("distvar", "gapdens", "nucdiv")
	stopifnot(is.character(measures))
	starttime = Sys.time()
	for (measure in measures){
		if (is.null(fun.usedist[[measure]])){ fun.usedist[[measure]]=FALSE }
		if (is.null(fun.userange[[measure]])){ fun.userange[[measure]]=FALSE }
	}
	if (dist.model %in% c("pairwise", "percentage")){ distfun = dist.gene }else{ distfun = dist.dna }
	alnlen = dim(aln)[2]
	halfwin = windowsize/2
	if (is.null(foci)){
		foci = sort(setdiff(ceiling(halfwin/step):(alnlen%/%step-ceiling(halfwin/step)+1)*step, 0))
	}
	N = length(foci)
	rollfun = function(k){ 
		i = foci[k]
		if (!quiet){ printProgressUpperMatrix(k, N, step=100, initclock=starttime) }
		a = max((i-halfwin+1), 1)
		b = min(i+halfwin, alnlen)
		stopifnot(a<=b)
		alnrange = a:b
		if (!is.null(subsample)){
			# expects a list with component `$sites` (a vector of integers, indexes of sites in the full alignment) and `$maxsize` (an integer, maximum number of site to sample per window)
			subsize = subsample$maxsize
			alnrange = fixwidthWindowSubsampleSites(alnrange, subsample$sites, subsample$maxsize, distribute=T)
		}else{ subsize = windowsize }
		subaln = aln[, alnrange]
		subalnsize = dim(subaln)[2]
		if (rmAbsSeq){
			rmseq = which(absentSeq(subaln, propgap=propgap, relpropgap=relpropgap))
			if (length(rmseq)>0) subaln = subaln[-rmseq, ]
		}
		if (TRUE %in% fun.usedist | "distvar" %in% measures){ 
			dd = distfun(subaln, dist.model, pairwise.deletion=TRUE)
			d = dd[!is.na(dd)]
		}
		if (("chisq.spacing.sites" %in% measures) | ("wilcox.all.site.dist" %in% measures)){
			# expected spacing of (subsampled) sites under uniformity
			expspacing = windowsize/subsize
			# expected distribution of (subsampled) sites under uniformity
			expindexs = seq(0, windowsize, length.out=subsize)
			# expected distribution distances under uniform spacing of (subsampled) sites (upper triangular matrix of inter-site distances)
			expalldist = unlist(sapply(1:(subsize-1), function(i){ sapply((i+1):subsize, function(j){ expindexs[j] - expindexs[1] }) }))
		}
#~ 		print(subaln)
		results = sapply(measures, function(measure){
#~ 			print(measure)
			if (measure=="nucdiv"){         			result = pegas::nuc.div(subaln, pairwise.deletion=TRUE)
			}else{ if (measure=="gapdens"){ 			result = gapDensity(subaln)
			}else{ if (measure=="distvar"){ 			result = var(d, na.rm=TRUE)
			}else{ if (substring(measure, 1,2)=="LD"){  if ( ctrl.subalnsize & subalnsize==subsize ){
															result = mean(linkageDisequilibrium(as.character(subaln), metric=substring(measure, 3),
														 discard.gaps=FALSE, multiproc=1, quiet=TRUE), na.rm=TRUE)
														}else{ result = NA }
			}else{ if (measure=="PHI"){					if ( ctrl.subalnsize & subalnsize==subsize ){
															result = callPhi(subaln, sprintf("%s_%d-%d.fas", alnname, a, b), as.integer(subsize/10), tempdirpath=tempdirpath)
														}else{ result = rep(NA, 2) }
														names(result) = c('statistic', 'p.value')
			}else{ if (measure=="pscore"){				if ( ctrl.subalnsize & subalnsize==subsize ){
															pdaln = phyDat(subaln)
															# avoid printing the advance of parsimony search
															sink('/dev/null')
															ptree = pratchet(data=pdaln, maxit=100)
															# stop sinking
															sink(NULL)
															result = parsimony(ptree, pdaln)
														}else{ result = NA }
			}else{ if (measure=="chisq.spacing.sites"){
														# measure the departure of consecutive site spacing from uniformity
														if ( ctrl.subalnsize & !is.null(subsample) & subalnsize==subsize ){
															obsspacing = sapply(2:length(alnrange), function(i){alnrange[i] - alnrange[i-1]})
															chisqdeviates = sapply(obsspacing, function(x){ (x - expspacing)^2 / expspacing })
															result = sum(chisqdeviates) / subalnsize
														}else{ result = NA }
			}else{ if (measure=="wilcox.all.site.dist"){
														# test the difference of distribution of all inter-site distances between observed and expected under uniform (subsampled) site spacing 
														if ( ctrl.subalnsize & !is.null(subsample) & subalnsize==subsize ){
															centerindexes = alnrange - min(alnrange)
															obsalldist = unlist(sapply(1:(subalnsize-1), function(i){ sapply((i+1):subalnsize, function(j){ centerindexes[j] - centerindexes[1] }) }))
															result = unlist(wilcox.test(obsalldist, expalldist)[c('statistic', 'p.value')])
														}else{ result = rep(NA, 2) ; names(result) = c('statistic.W', 'p.value')}
			}else{ 
				if (!is.null(fun)){ 
					measfun = fun 
				}else{ 
#~ 					print(sys.frame())
					measfun = get(measure, inherits=TRUE)
					warnings()
					stopifnot(class(measfun)=="function")
				}
				if (fun.usedist[[measure]]){			indata = d 
				}else{ if (fun.userange[[measure]]){	indata = alnrange
				}else{									indata = subaln 
				}}
				res = measfun(indata)
				if (measure %in% names(fun.res)){ result = res[[fun.res[[measure]]]]
				}else{ result = res }
			}}}}}}}}
			return(result)			
		})
		if (is.list(results)){
			results = c(results, recursive=T)
			measnames = names(results)
#~ 			print(results)
#~ 			print('')
#~ 			measnames = character()
#~ 			for (i in length(results)){
#~ 				if ( !is.null(names(results[[i]])) ){
#~ 					measnames = c(measnames, paste(names(results)[i], names(results[[i]])))
#~ 				}else{
#~ 					measnames = c(measnames, names(results)[i])
#~ 				}
#~ 			}
		}else{ if (is.matrix(results) | is.data.frame(results)){
			stopifnot(1 %in% dim(results))
			d = which(dim(results)==1)[1]
			measnames = paste(colnames(results), rownames(results))
			if (d==1){ results = as.vector(results[1,])
			}else{ results = as.vector(results[,1]) }	
		}else{ measnames = measures }}
#~ 		print(results)
#~ 		print(measnames)
		stopifnot(is.vector(results), length(results)==length(measnames))
		names(results) = measnames
		return(results)
	}
	if (multiproc > 1){ 
		roll = simplify2array(mclapply(1:N, rollfun, mc.cores=multiproc, mc.preschedule=FALSE))
	}else{  
		roll = sapply(1:N, rollfun)
	}
	if (length(measures)>1){ roll = as.data.frame(t(roll))
	}else{ roll = as.data.frame(roll) ; colnames(roll) = measures[1] }
	roll = cbind(data.frame(foci=foci), roll)
	if (!quiet){ cat("\n") }
	return(roll)
}

getHyperVarClusterProfile = function(aln, windowsize=100, step=25, ldv.thresh=-3, dist.thresh=0.05, dist.model='raw', multiproc=1){
	stopifnot(class(dist.thresh) %in% c("numeric", "character", "function"))
	if (dist.model %in% c("pairwise", "percentage")){ distfun = dist.gene }else{ distfun = dist.dna }
	rdistvar = rollStats(aln, windowsize=windowsize, step=step, measures=c("distvar"), dist.model=dist.model, multiproc=multiproc)
	alnlen = dim(aln)[2]
	print(summary(rdistvar))
	print(paste("threshold: log10(distvar) > ", ldv.thresh, sep=''))
	foci = rdistvar$foci[(log(rdistvar$distvar, 10) > ldv.thresh) & !is.na(rdistvar$distvar)]
	N = length(foci)
	halfwin = windowsize/2
	if (!is.numeric(dist.thresh)){ 
		if (is.character(dist.thresh)){ 
			threshfun = get(dist.thresh) 
			stopifnot(class(threshfun)=="function")
		}else{ threshfun = dist.thresh }
	}else{ threshfun = NULL }
	roll = simplify2array(mclapply(1:N, function(k){ 
#~ 	roll = sapply(1:N, function(k){ 
		i = foci[k]
		printProgressUpperMatrix(k, N, step=100)
		alnrange = (i-halfwin+1):min(i+halfwin, alnlen)
		dd = distfun(aln[, alnrange], dist.model, pairwise.deletion=TRUE)
		dth = ifelse(is.null(threshfun), dist.thresh,  threshfun(dd[!is.na(dd)]))
#~ 		clustering = hclust(dd, method='single')
#~ 		cutree(clustering, h=dth)
		if (is.na(dth)){ print(table(dd, useNA='ifany')) }
		cl = singleLinkageClustersFromMatrix(dd, dist.thresh=dth)
		return(cl)
	}, mc.cores=multiproc, mc.preschedule=TRUE))
#~ 	})
	cat("\n")
	colnames(roll) = foci
	return(roll)
}

codeHyperVarClusterProfile = function(hvprofile, hvreg=NULL, coding="AA", write.full=NULL, write.reg=NULL, return.fmt="matrix"){
	aminoacids = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
	if (coding=="AA"){ code = aminoacids
	}else{ code = coding }
	if (is.list(hvprofile)){ 
		lprofreg = hvprofile 
		seqnames = unique(c(lapply(hvprofile, rownames), recursive=TRUE))
	}else{
		stopifnot(is.matrix(hvprofile))
		lprofreg = lapply(hvreg, function(reg){ hvprofile[, as.numeric(colnames(hvprofile)) %in% reg[1]:reg[2]] }) 
		names(lprofreg) = sapply(hvreg, paste, collapse='-')
		seqnames = rownames(hvprofile)
	}
	print(seqnames)
	hvprofile.coded = do.call(cbind, lapply(names(lprofreg), function(nprofreg){
		profreg = lprofreg[[nprofreg]]		
#~ 		codhv = apply(as.matrix(profreg), 2, function(v){ code[v] })
#~ 		codhv = ifelse(is.na(codhv), '-', codhv)
		lseqcodhv = lapply(seqnames, function(nsq){ 
			if (nsq %in% rownames(profreg)){
				sq = profreg[nsq,]
				return(ifelse(is.na(sq), '-', code[sq]))
			}else{ 
				return(rep('-', dim(profreg)[2]))
			}
		})
		if (!is.null(write.reg)){
			if (length(write.reg) > 1){ hvclustdir = write.reg[1] ; hvclustprefix = write.reg[2] 
			}else{ hvclustdir = '.' ; hvclustprefix = write.reg[1] }
			seqinr::write.fasta(lseqcodhv, names=rownames(profreg), file.out=file.path(hvclustdir, paste(hvclustprefix, nprofreg, 'fasta', sep='.')))
		}
		return(t(sapply(lseqcodhv, function(x){c(x, '-')})))
	}))
	rownames(hvprofile.coded) = seqnames
	hvprofile.coded.list = lapply(seqnames, function(nsq){ hvprofile.coded[nsq,] })
	if (!is.null(write.full)){
		if (length(write.full) > 1){ hvclustdir = write.full[1] ; hvclustprefix = write.full[2] 
		}else{ hvclustdir = '.' ; hvclustprefix = write.full[1] }
		print(paste(hvclustdir, paste(hvclustprefix, 'fasta', sep='.'),  sep='/'))
		seqinr::write.fasta(hvprofile.coded.list, names=seqnames, file.out=paste(hvclustdir, paste(hvclustprefix, 'fasta', sep='.'),  sep='/'))
	}
	if (return.fmt=="phyDat"){ 
		hvprofile.coded.pd = as.phyDat(hvprofile.coded, type="AA")
		names(hvprofile.coded.pd) = rownames(hvprofile)
		return(hvprofile.coded.pd) 
	}else{
		return(hvprofile.coded)
	}
}


fullAln2RefCoords = function(fullaln, reflabel=1, bijective=TRUE){
	k = 0
	map = c()
	for (j in 1:dim(fullaln)[2]){
		if (fullaln[reflabel,j]!=as.DNAbin('-')){
			k = k + 1
			map = c(map, k)
		}else{
			if (bijective){ map = c(map, NA)
			}else{ map = c(map, k) }
		}
	}
#~ 	imap = sapply(1:dim(fullaln)[2], function(j){
#~ 		if (fullaln[reflabel,j]!=as.DNAbin('-')){
#~ 			k = k + 1
#~ 			fwd = TRUE
#~ 		}else{
#~ 			fwd = FALSE
#~ 		}
#~ 		return(c(fwd, k))
#~ 	})
#~ 	map = imap[imap[,1]==TRUE,2]
#~ 	if (bijective){ map[map==0] = NA }
	return(map)
}

getAlnFmt = function(seq){
	wrongfmt = 'unvalid format, needs DNAbin or character matrix single-sequence alignment'
	if (!(is.matrix(seq) & dim(seq)[1]==1)){ stop(wrongfmt) }
	if (class(seq)=='DNAbin'){ asf = as.DNAbin
	}else{ if (class(refalnseq[1,1])=='character'){ asf = as.character
	}else{ stop(wrongfmt)
	}}
	return(asf)
}

translateAlnCoords = function(refalnseq, newcoords, asf=as.DNAbin){
		transalnseq = as.character(refalnseq)[,newcoords]
		transalnseq[is.na(transalnseq)] = '-'
		transalnseq = asf(transalnseq)
		return(transalnseq)
}

refAln2NewAlnCoords = function(refalnseq, newalnseq, returnTranslatedAln=F){
	# same fullAln2RefCoords(), but can handle an alredy aligned sequence as reference
	# it detects where gaps have been INTRODUCED (cannot deal with gap removal)

	as1 = getAlnFmt(refalnseq)
	as2 = getAlnFmt(newalnseq)
	
	k = 1
	map = c()
	for (j in 1:dim(newalnseq)[2]){
		if (newalnseq[1,j]==as2('-')){
			if (refalnseq[1,k]==as1('-')){
				map = c(map, k)
				k = k + 1
			}else{ 
				# gap was introduced in newalnseq
				map = c(map, NA)
			}
		}else{
			map = c(map, k)
			k = k + 1
		}
	}
	if (!returnTranslatedAln){
		return(map)
	}else{
		# for verification purpose
		transalnseq = translateAlnCoords(refalnseq, map, asf=as2)
		# so that: transalnseq == newalnseq
		return(transalnseq)
	}
}

regionList2indexes = function(reglist, concat=TRUE){
	regs = lapply(reglist, function(x){
		begreg = x[1]
		endreg = x[2]
		return(begreg:endreg)
	})
	if (concat){
		return(do.call('c', regs))
	}else{
		return(regs)
	}
}

regionIndexes2List = function(indexes, step=1, windowsize=0, alnlen=Inf){
	halfwin = windowsize/2
	lout = list()
	indexes = sort(indexes)
	curr = init = indexes[1]
	for (n in 2:length(indexes)){
		reg = list(c(max(1, init-halfwin+1), min(curr+halfwin, alnlen)))
		if (indexes[n] > curr+step+halfwin){
			lout = c(lout, reg)
			curr = init = indexes[n]
		}else{
			curr = indexes[n]
		}
	}
	lout = c(lout, reg)
	return(lout)
}

# for ploting region lists as from regionIndexes2List()
plotbound = function(pair, ybt=par()$usr[3:4], ...){
	rect(xleft=pair[1], xright=pair[2], ybottom=ybt[1], ytop=ybt[2], density=NA, ...)
}
plotbound0 = function(pair, ybt=par()$usr[3:4], coul='grey'){
	rect(xleft=pair[1], xright=pair[2], ybottom=ybt[1], ytop=ybt[2], density=NA, col=coul)
}

writeSiteLocs = function(sitelocs, lenaln, nflocs){
	write(c(length(sitelocs), lenaln, "L"), file=nflocs, ncolumns=3)
	write(sitelocs, file=nflocs, append=TRUE, ncolumns=1)
}

writeSubAln = function(reg, fullaln, fileprefix='', sites='all', exclude.seq=NULL, format='fasta'){
	begreg = reg[1]
	endreg = reg[2]
	regaln = fullaln[, begreg:endreg]
	if (!is.null(exclude.seq)){ regaln = regaln[-exclude.seq,] }
	if (sites!='all'){
		if (tolower(sites)=='snp'){ is = getSNPindexes(regaln) }
		if (tolower(sites)=='parsinform'){ is = getParsInformIndexes(regaln) }
		if (tolower(sites)=='biallelic'){ is = getBiAllelicIndexes(regaln) }
		if (length(is)>0){
			regaln = regaln[,is] 
			locs = is + begreg -1
		}else{
			regaln = NULL 
			locs = NULL
		}
	}else{
		locs = NULL
	}
	if (!is.null(regaln)){
		nfpref = paste(fileprefix, begreg, "-", endreg, sep="")
		nfaln = paste(nfpref, ifelse(format %in% c('interleaved', 'sequential'), 'phylip', format), sep=".")
		write.dna(regaln, nfaln, format=format)
		if (!is.null(locs)){
			nflocs = paste(nfpref, 'locs', sep=".")
			writeSiteLocs(locs, dim(fullaln)[2], nflocs)
		}
	}
	l = list(reg, regaln, locs)
	names(l) = c('coords', 'aln', 'locs')
	return(l)
}

#### dealing with binary coding of bipartitions as in MrBayes output

#~ str2set = function(s, sep=','){ sort(unique(strsplit(as.character(s), split=sep)[[1]])) }
#~ 
#~ distclades = function(cla1, cla2){
#~ 	s = str2set(cla1, ',')
#~ 	S = str2set(cla2, ',')
#~ 	return(max(length(s), length(S)) - length(intersect(s,S)))
#~ }
#~ 
#~ distbiparts = function(bip1, bip2){
#~ 	s = str2set(bip1, '|')
#~ 	S = str2set(bip2, '|')
#~ 	mdist = sapply(s, function(cla){
#~ 		sapply(S, function(Cla){
#~ 			distclades(cla, Cla)
#~ 		})
#~ 	})
#~ 	return(min(mdist))
#~ }
#~ 
#~ 
#~ fuzziness = function(bip1, bip2){
#~ 	if (is.na(cla1) | is.na(cla2)){ return(NA)}
#~ 	else{
#~ 		s = str2set(bip1, '|')
#~ 		S = str2set(bip2, '|')
#~ 		return(distbiparts(bip1, bip2)/min(length(s[1]), length(S[1])))
#~ 	}
#~ }

dcat = function(x, y){
	stopifnot(length(x)==length(y))
	sum(as.numeric(sapply(1:length(x), function(i){ x[i]!=y[i] })))
}
distcat = function(dmat){
	d = sapply(1:dim(dmat)[1], function(i){
		sapply(1:dim(dmat)[1], function(j){
			dcat(dmat[i,], dmat[j,])
		})
	})
	dimnames(d) = list(rownames(dmat), rownames(dmat))
	return(as.dist(d))
}

translateCladeToBipartBitaray = function(taxset, ltax){ paste(as.numeric(ltax %in% taxset), collapse='') }

identifyBipart = function(taxset, ltax, dbbiparts){  dbbiparts$profile==translateCladeToBipartBitaray(taxset, ltax) }

translateBipartToTaxSet = function(bip, ltax, value='bipart', as.vector=FALSE, name.dict=NULL, sep.strains=', ', sep.clades=' | '){
	prof = as.numeric(strsplit(bip, split='')[[1]])
	c1 = character()
	if (value=='bipart'){ 
		c2 = character()
		for (i in 1:length(prof)){
			t = ltax[i]
			if (!is.null(name.dict)){
				t = name.dict[name.dict[,1]==t,2] 
				stopifnot(length(t)>0)
			}
			if (prof[i]==1){ c1 = c(c1, t)
			}else{ if (prof[i]==0){ c2 = c(c2, t) 
			}else{ stop('unproper bipart coding, need 0 or 1') 
			}}
		}
		if (as.vector){ clades = list(c1, c2)
		}else{ clades = c(paste(c1, collapse=sep.strains), paste(c2, collapse=sep.strains)) }
	}else{ if (value=='smallclade'){ 
		for (i in 1:length(prof)){
			t = ltax[i]
			if (!is.null(name.dict)){
				t = name.dict[name.dict[,1]==t,2] 
				stopifnot(length(t)>0)
			}
			if (prof[i]==1) c1 = c(c1, t)
		}
		if (as.vector){ clades = c1
		}else{ clades = paste(c1, collapse=sep.strains) }
	}}
	if (as.vector){ return(clades) 
	}else{ return(paste(clades, collapse=sep.clades)) }
	
}

getSupportedClades = function(lsupbip, ltax, excl.strains=NULL){
	# use only a set of bipart with at least 0.5 support ; i.e. there are no conflicting bipart in the set
	lsupclades = do.call(c, lapply(lsupbip, translateBipartToTaxSet, ltax=ltax, value='bipart', as.vector=TRUE))
	lsupclades = lsupclades[order(sapply(lsupclades, length))]
	lnrsupclades = list(lsupclades[[1]])
	for (k in 2:length(lsupclades)){
		lsc = lsupclades[[k]]
		for (lnrsc in lnrsupclades){
			lsc = setdiff(lsc, lnrsc)
			lsc = setdiff(lsc, excl.strains)
		}
		if (length(lsc)>0){ lnrsupclades = c(lnrsupclades, list(lsc)) }
	}
	# put excluded strain in a separate clade, should be ignored in subsequent
	if (!is.null(excl.strains)){ lnrsupclades = c(lnrsupclades, list(excl.strains)) }
	return(lnrsupclades)
}

codeSequenceTypes = function(lnrsupclades, ltax, as.bin.matrix=FALSE, excl.strains=NULL, excl.as=NA){
	prof = sapply(ltax, function(tax){ which(sapply(lnrsupclades, function(lnrsc){ tax %in% lnrsc })) })
	if (!is.null(excl.strains)){
		# mark absent or excluded strains by attributing them e.g. NA or 0 as a sequence type
		prof[ltax %in% excl.strains] = excl.as
	}
	prof = unlist(prof)
	if (!as.bin.matrix){ return(prof)
	}else{ return(sapply(1:length(lnrsc), function(x){ prof==x })) }
}

testRVCorMat = function(lprofbinmat){
	lcorbinmat = sapply(names(lprofbinmat), function(i){
		lapply(names(lprofbinmat), function(j){
			rv = FactoMineR::coeffRV(lprofbinmat[[i]], lprofbinmat[[j]])
		})
	})
	return(lcorbinmat)
}

toBinProfile = function(matprofseqtypes){
	lprofbinmat = lapply(colnames(matprofseqtypes), function(n){
		prof = matprofseqtypes[,n]
		sapply(min(prof):max(prof), function(x){ as.numeric(prof==x) })
	})
	names(lprofbinmat) = colnames(matprofseqtypes)
	return(lprofbinmat)
}

testCorSeqTypes = function(matprofseqtypes){
	lprofbinmat = toBinProfile(matprofseqtypes)
	lcorbinmat = testRVCorMat(lprofbinmat)
	dimnames(lcorbinmat) = list(colnames(matprofseqtypes), colnames(matprofseqtypes))
	return(lcorbinmat)
}

exportSeqTypesAlignments = function(matprofseqtypes, laln, dir.out=NULL, gname.long.alias=NULL){
	if (!is.null(dir.out)){
		diroutseqs = file.path(dir.out, 'seqtype_sequences')
		dir.create(diroutseqs, showWarnings=F, recursive=T)
		diroutcons = file.path(dir.out, 'seqtype_consensus')
		dir.create(diroutcons, showWarnings=F, recursive=T)
	}
	lcons = sapply(colnames(matprofseqtypes), function(gname){
		profseqtypes = matprofseqtypes[,gname]
		lstcons = list() ; lstaln = list() ; nstcons = character()
		if (!is.null(gname.long.alias)){ genename = gname.long.alias[gname] }else{ genename = gname }
		for (st in levels(as.factor(profseqtypes))){
			if (as.numeric(st)!=0){
				staln = laln[[gname]][rownames(matprofseqtypes)[profseqtypes==st], ]
#~ 				print(staln)
				if (!is.null(dir.out)){ write.dna(staln, file=file.path(diroutseqs, paste(genename, '-ST', st, '.aln', sep='')), format='fasta') }
				cons = seqinr::consensus(ape::as.alignment(staln), method='majority')
				lstcons = c(lstcons, list(cons))
				lstaln = c(lstaln, list(staln))
#~ 				nstcons = c(nstcons, paste(genename, '-ST', st, sep=''))
				nstcons = c(nstcons, paste('ST', st, sep=''))
			}
		}
		if (!is.null(dir.out)){ seqinr::write.fasta(lstcons, names=nstcons, nbchar=max(sapply(lstcons, length), na.rm=T),
		 file.out=file.path(diroutcons, paste(genename, '-STconsensus', '.aln', sep=''))) }
		mstcons = as.DNAbin(t(simplify2array(lstcons)))
		rownames(mstcons) = nstcons
		names(lstaln) = nstcons
		return(list(mstcons, lstaln))
	})
	return(lcons)
}

compatibleBiparts = function(lbip, ltax, nbcores=6, quiet=FALSE){
	# build list of clades implied by all known bipartitons
	matbipclades = t(sapply(lbip, translateBipartToTaxSet, ltax=ltax, value='bipart', as.vector=TRUE))
	# matrix of compatible bipartitions
	# <=> at least one of the implied clades of the first bipart can be found included in one of the clades of the second
	N = dim(matbipclades)[1]
#~ 	matcomp = simplify2array(mclapply(1:N, function(i){
#~ 		if (i < N){
#~ 			upper = sapply((i+1):N, function(j){
#~ 				any(sapply(1:2, function(m){ sapply(1:2, function(n){ all(matbipclades[i,m][[1]] %in% matbipclades[j,n][[1]]) }) }))
#~ 				
#~ 			})
#~ 		}else{ upper = logical(0) }
#~ 		# complete upper triangle matrix with NA's in lower traingle and TRUE in the diagonal
#~ 		return(c(rep(NA, i-1), TRUE, upper) 
#~ 	}, mc.cores=nbcores, mc.preschedule=FALSE))
#~ 	# replace the NA's in the lower traingle by transpose value
#~ 	mclapply(1:N, function(j){sapply((j+1):N, function(i){ 
#~ 		matcomp[j,i] = matcomp[i,j]
#~ 	})}, mc.cores=nbcores, mc.preschedule=FALSE)
	matcomp = simplify2array(mclapply(1:N, function(i){
#~ 	matcomp = sapply(1:N, function(i){
		if (!quiet){ printProgressUpperMatrix(i, N, step=10) }
		sapply(1:N, function(j){
			any(sapply(1:2, function(m){ sapply(1:2, function(n){ all(matbipclades[i,m][[1]] %in% matbipclades[j,n][[1]]) }) }))
		})
	}, mc.cores=nbcores, mc.preschedule=TRUE))
#~ 	})

}


getSupportedBlocks = function(PPbiparts, gapsize=1, ppthresh=0.35, subset.bip=NULL, gene.names=NULL, min.report=1){
	biparts = numeric()
	firstgene = numeric()
	lastgene = numeric()
	sizetrack = numeric()
	if (!is.null(subset.bip)) bips = subset.bip
	else bips = 1:dim(PPbiparts)[1]
	for (i in bips){
		k = 0
		g = 0
		for(j in 1:dim(PPbiparts)[2]){
			if (PPbiparts[i, j]>=ppthresh){
				if (k==0){ firstgene = c(firstgene, j) }
#~ 				print(c('here', k, g))
				k = k + 1
				g = 0
			}else{ if (k>0){
				if (g<gapsize){
#~ 					print(c('la', k, g))
					g = g + 1
					k = k + 1
				}else{
#~ 					print(c('qui', k, g))
					biparts = c(biparts, i)
					lastgene = c(lastgene, j-g-1)
					sizetrack = c(sizetrack, k-g)
					k = 0
					g = 0
			}}}
		}
		if (k>0){
#~ 			print(c('aqui', k, g))
			biparts = c(biparts, i)
			lastgene = c(lastgene, j-g-1)
			sizetrack = c(sizetrack, k-g)
		}
#~ 		print(length(biparts))
#~ 		print(length(firstgene))
#~ 		print(length(lastgene))
#~ 		print(length(sizetrack))
	}
	if (!is.null(rownames(PPbiparts))){
		biparts = rownames(PPbiparts)[biparts]
	}
	if (!is.null(gene.names)){
		firstgene = gene.names[firstgene]
		lastgene = gene.names[lastgene]
	}
#~ 	print(length(biparts))
#~ 	print(length(firstgene))
#~ 	print(length(lastgene))
#~ 	print(length(sizetrack))
	dfbip = data.frame(bipart=biparts, first.gene=firstgene, last.gene=lastgene, track.size=sizetrack)
	return(dfbip[dfbip$track.size>=min.report,])
}
