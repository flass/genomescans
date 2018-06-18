#!/usr/bin/Rscript --vanilla
library('getopt')

# module to summarize LD scans obtained by resampling (simulations) of genomes in a test dataset as done with genome-wide_LD_scan.r script

tablename = 'LocalLD-subsampled.minalfrq1.biallelicsites.max1gaps'
metrics = c('meanldrsub', 'logcompldfisub') ; names(metrics) = c('r2', 'Fisher')
mergecols = c('foci', 'reference.position')


sample.size = function(x, na.rm=F){
	n = length(which(!is.na(x)))
	return(ifelse(n>0, n, NA))
}

loadLDrollTables = function(ldtabdir, file.rad=tablename, min.snp.density=20){	
	ldrolltables = lapply(names(metrics), function(LDmetric){
		if (file.exists(file.path(ldtabdir, sprintf('LD_%s.%s.RData', LDmetric, tablename)))){
			load(file.path(ldtabdir, sprintf('LD_%s.%s.RData', LDmetric, tablename)))
			return(list(ldrollsub=ldrollsub, rollsubsnpdens=rollsubsnpdens))
		}else{
			ldrollsub = read.table(file.path(ldtabdir, sprintf('LD_%s.%s.tab', LDmetric, tablename)), row.names=1, header=T)
			return(list(ldrollsub=ldrollsub))
		}
	})
	ldrollsub = merge(ldrolltables[[1]][['ldrollsub']], ldrolltables[[2]][['ldrollsub']], by=mergecols, all=T, sort=F)
	if (!is.null(ldrolltables[[1]][['rollsubsnpdens']])){
		ldrollsub = merge(ldrollsub, ldrolltables[[1]][['rollsubsnpdens']], by=mergecols, all=T, sort=F)
		ldrollsub = ldrollsub[ldrollsub$reportsnpdens >= min.snp.density, ]
	}
	return(ldrollsub)
}

getAverageAndCI = function(LDmetric, ldrollsubtables, n.sims){
	meanldrsubs = sapply(ldrollsubtables, function(x){ x[,metrics[LDmetric]] })
	ldstats = as.data.frame(sapply(list(mean=mean, median=median, sd=sd, mad=mad, sample.size=sample.size), function(fun){ apply(meanldrsubs, 1, fun, na.rm=T) }))
	meandev = 1.96*ldstats$sd/sqrt(ldstats$sample.size) ; ldstats$lower.mean.IC0.95 = ldstats$mean - meandev ; ldstats$upper.mean.IC0.95 = ldstats$mean + meandev
	mediandev = 1.96*1.4826*ldstats$mad/sqrt(n.sims) ; ldstats$lower.median.IC0.95 = ldstats$median - mediandev ; ldstats$upper.median.IC0.95 = ldstats$median + meandev
	colnames(ldstats) = paste(LDmetric, colnames(ldstats), sep='.')
	print(summary(ldstats))
	return(cbind(ldrollsubtables[[1]][,mergecols], ldstats))
}

loadAndCombineSimLDrollTables = function(resampldir, min.snp.density=20){
	sims = list.files(resampldir, full.names=T)
	simldrollsubtables = lapply(sims, loadLDrollTables, min.snp.density)
	simldmetables = lapply(names(metrics), getAverageAndCI, ldrollsubtables=simldrollsubtables, n.sims=length(sims))
	simldmetable = merge(simldmetables[[1]], simldmetables[[2]], by=mergecols, all=T, sort=F)
}

getCIbounds = function(ldmetable, averagefun='median', plot.CI=F, ...){
	lapply(c('lower', 'upper'), function(bound){ 
		CIbound = ldmetable[,sprintf('r2.%s.%s.IC0.95', bound, averagefun)] 
		if (plot.CI) lines(ldmetable$reference.position, CIbound, ...)
		return(CIbound)
	})
}

compareExp2Sim = function(exptable, simtable, averagefuns=c('mean', 'median')){
	# corresponding positions in the simulations
	expsimpos = sapply(exptable$reference.position, function(x){ 
		i = which(simtable$reference.position==x); if (length(i)==1){ return(i) }else{ return(NA) }
	})
	for (average in averagefuns){
		CIbounds = getCIbound(simtable, averagefun=averagefun)
		exptable[,sprintf('in.sim.%s.IC0.95', average)] = (exptable$meanldrsub < CIbounds[[1]][expsimpos] | exptable$meanldrsub > CIbounds[[2]][expsimpos])
	}
	invisible(exptable)
}

plotGenomicMapStatSummary = function(simldmetable, averagefuns=c('mean', 'median'), opt, n.sims, expldrollsub=NULL){
	for (average in averagefuns){
		par(mar=c(8,8,8,8))
		plot(simldmetable$reference.position, simldmetable[,sprintf('r2.%s', average)],
		 col='white',
		 ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", opt$windowsize), xlab="reference genome coordinates",
		 main=sprintf("%s (%d combined simulations) - LD scan with fixed-size windows", opt$datasetname, n.sims)
		)
		abline(h=1:10*(1/10), col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
		CIbounds = getCIbound(simldmetable, averagefun=averagefun, plot.CI=T, col='slategrey')
		signifpoints = simldmetable[,sprintf('Fisher.%s', average)] > opt$signifthresh
		lines(simldmetable$reference.position, simldmetable[,sprintf('r2.%s', average)], col='blue')
		points(simldmetable$reference.position[signifpoints], simldmetable[signifpoints, sprintf('r2.%s', average)], pch=20, col='purple')
		legtext = c(sprintf("%s of %d simulations", average, n.sims),
					"(when significantly higher than the the genomic basal level)", sprintf("95%% confidence interval around %s", average))
		legcol = c('blue', 'purple', 'slategrey')
		leglwd = c(1, 3, 1)
		legend('topleft', legend=legtext, col=legcol, lwd=leglwd, bg='white')
		if (!is.null(expldrollsub)){
			# plot small dots/open rounds when in/out the range of simulated data; colour them black/red when non-sig/significant
			plot(x=expldrollsub$reference.position, y=expldrollsub$meanldrsub,
			 ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", opt$windowsize), xlab="reference genome coordinates",
			  main=sprintf("%s - LD scan with fixed-size windows", basename(opt$expdir)), col='white')
			abline(h=1:10*(1/10), col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
			points(x=expldrollsub$reference.position, y=expldrollsub$meanldrsub,
			 col=ifelse(expldrollsub$logcompldfisub > opt$signifthresh, 'red', 'grey'),
			 pch=ifelse(expldrollsub[,sprintf('in.sim.IC0.95', average)], 1, 46),
			 cex=0.5)
			legtext = c(sprintf("out of 95%% confidence interval around simul %s", average),
						 sprintf("within 95%% confidence interval around simul %s", average),
						 "significantly higher than the genomic basal level",
						 "not significantly higher")
			legpch = c(1, 46, 20, 20)
			legcol = c('black', 'black', 'red', 'grey')
			legend('topleft', legend=legtext, pch=legpch, col=legcol, bg='white')
		}
	
	}
}

SummarizeExpAndSimulGenomeWideLD = function(opt, file.rad=tablename){

	simldmetable = loadAndCombineSimLDrollTables(opt$resampldir, file.rad=file.rad, min.snp.density=opt$nbsnp)
	if (!is.null(opt$outtab)){ 
		write.table(simldmetable, file=opt$outtab, row.names=F, col.names=T, quote=F, sep='\t')
	}
	if (!is.null(opt$expdir)){
		# additional data to be compared to the simulations
		expldrollsub = loadLDrollTables(opt$expdir, file.rad=file.rad, min.snp.density=opt$nbsnp)
		compareExp2Sim(expldrollsub, simldmetable, averagefun=average)
	}else{ expldrollsub = NULL }
	
	if (!is.null(opt$outpdf)){ 
		### genomic map plot of stat summary
		pdf(file=opt$outpdf, width=15, height=10,
		 title=sprintf("%s - local LD scan - windows %dbp", opt$datasetname, opt$windowsize))
		plotGenomicMapStatSummary(ldmetable, averagefuns=c('mean', 'median'), opt=opt, n.sims=length(sims), expldrollsub=expldrollsub)
		dev.off()
	}	
	invisible( list(simldmetable=simldmetable, expldrollsub=expldrollsub) )
}

print(environment())
if (identical(environment(), globalenv())){
	# if loaded as a module (source file for functions),
	# to avoid executing the main function, one has to be called with:
	# source('/path/to/summarize_simul_genome-wide_LD.r', local=new.env())

	spec = matrix(c(
	  'resampldir',   'S', 1, "character", "path of folder with simulated LD scan",
	  'expdir',       'E', 2, "character", "path of folder with experimental LD scan to compare simulation with [optional]",
	  'windowsize',   'w', 2, "integer",   "physical size (bp) of the sliding windows in which LD is evaluated [default: 3000]",
	  'step',         's', 2, "integer",   "step (bp) of the sliding windows in which LD is evaluated [default: 10]",
	  'nbsnp',        'm', 2, "integer",   "number of biallelic SNP within each window that are used for LD computation\n(windows with less than that are excluded from the report) [default: 20]",
	  'signifthresh', 't', 2, "double",    "threshold of Local LD Index (deviation form Fisher's exact test p-value from the genomic norm) above which the observed LD is significant [default: 5]",
	  'outtab',       'a', 2, "character", "path of ouput table file for summary of simulated values",
	  'outpdf',       'o', 2, "character", "path of ouput PDF file for plots"
	), byrow=TRUE, ncol=5);
	opt = getopt(spec, opt=commandArgs(trailingOnly=T))
	
	if ( is.null(opt$windowsize  ) ){ opt$windowsize   = 3000  }
	if ( is.null(opt$step        ) ){ opt$step         = 10    }
	if ( is.null(opt$nbsnp       ) ){ opt$nbsnp        = 20    }
	if ( is.null(opt$signifthresh) ){ opt$signifthresh = 5     }
	if ( is.null(opt$datasetname ) ){ opt$datasetname = basename(opt$resampldir)}
	if ( is.null(opt$outtab      ) ){ opt$outtab = sprintf("%s/%s_summary.tab", dirname(opt$resampldir), opt$datasetname)}
	if ( is.null(opt$outpdf      ) ){ opt$outpdf = sprintf("%s/%s_summary%s.pdf", dirname(opt$resampldir), opt$datasetname, ifelse(!is.null(opt$expdir), sprintf("_vs_%s", basename(opt$expdir)), ""))}

	SummarizeExpAndSimulGenomeWideLD(opt)
}
	
