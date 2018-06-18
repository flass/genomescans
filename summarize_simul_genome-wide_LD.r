#!/usr/bin/Rscript --vanilla
library('getopt')

tablename = 'LocalLD-subsampled.minalfrq1.biallelicsites.max1gaps'
metrics = c('meanldrsub', 'logcompldfisub') ; names(metrics) = c('r2', 'Fisher')
mergecols = c('foci', 'reference.position')

loadLDrollTables = function(ldtabdir, maxsize=20){	
	ldrolltables = lapply(names(metrics), function(LDmetric){
		if (file.exists(file.path(ldtabdir, sprintf('LD_%s.%s.RData', LDmetric, tablename)))){
			load(file.path(ldtabdir, sprintf('LD_%s.%s.RData', LDmetric, tablename)))
			return(list(ldrollsub=ldrollsub, rollsubsnpdens=rollsubsnpdens))
		}else{
			ldrollsub = read.table(file.path(ldtabdir, sprintf('LD_%s.%s.tab', LDmetric, tablename)), row.names=1, header=T)
			# mock values for biallelic SNP density, set at the threshold so to pass the filter
			rollsubsnpdens = data.frame(reportsnpdens=rep(maxsize, nrow(ldrollsub)), foci=ldrollsub$foci, reference.position=ldrollsub$reference.position)
			return(list(ldrollsub=ldrollsub, rollsubsnpdens=rollsubsnpdens))

		}
	})
	ldrollsub = merge(ldrolltables[[1]][['ldrollsub']], ldrolltables[[2]][['ldrollsub']], by=mergecols, all=T, sort=F)
	ldrollsub = merge(ldrollsub, ldrolltables[[1]][['rollsubsnpdens']], by=mergecols, all=T, sort=F)
	return(ldrollsub)
}

spec = matrix(c(
  'resampldir',   'S', 1, "character", "path of folder with simulated LD scan",
  'expdir',       'E', 2, "character", "path of folder with experimental LD scan to compare simulation with [optional]",
  'windowsize',   'w', 2, "integer",   "physical size (bp) of the sliding windows in which LD is evaluated [default: 3000]",
  'step',         's', 2, "integer",   "step (bp) of the sliding windows in which LD is evaluated [default: 10]",
  'maxsize',      'm', 2, "integer",   "maximum number of biallelic SNP within each window that are used for LD compuation\n(windows with less than that are excluded from the report) [default: 20]",
  'signifthresh', 't', 2, "double",    "threshold of Local LD Index (deviation form Fisher's exact test p-value from the genomic norm) above which the observed LD is significant [default: 5]",
  'outtab',       'a', 2, "character", "path of ouput table file for summary of simulated values",
  'outpdf',       'o', 2, "character", "path of ouput PDF file for plots"
), byrow=TRUE, ncol=5);
opt = getopt(spec, opt=commandArgs(trailingOnly=T))

if ( is.null(opt$windowsize  ) ){ opt$windowsize   = 3000  }
if ( is.null(opt$step        ) ){ opt$step         = 10    }
if ( is.null(opt$maxsize     ) ){ opt$maxsize      = 20    }
if ( is.null(opt$signifthresh) ){ opt$signifthresh = 5     }
if ( is.null(opt$datasetname ) ){ opt$datasetname = basename(opt$resampldir)}
if ( is.null(opt$outtab      ) ){ opt$outtab = sprintf("%s/%s_summary.tab", dirname(opt$resampldir), opt$datasetname)}
if ( is.null(opt$outpdf      ) ){ opt$outpdf = sprintf("%s/%s_summary%s.pdf", dirname(opt$resampldir), opt$datasetname, ifelse(!is.null(opt$expdir), sprintf("_vs_%s", basename(opt$expdir)), ""))}
sims = list.files(opt$resampldir)
ldrollsubtables = lapply(sims, function(sim){ loadLDrollTables(file.path(opt$resampldir, sim), maxsize=opt$maxsize) })
ldmetables = lapply(names(metrics), function(LDmetric){
	meanldrsubs = sapply(ldrollsubtables, function(x){ ifelse(x$reportsnpdens >= opt$maxsize, x[,metrics[LDmetric]], NA) })
	sample.size=function(x, na.rm=F){ n = length(which(!is.na(x))) ; ifelse(n>0, n, NA) }
	ldstats = as.data.frame(sapply(list(mean=mean, median=median, sd=sd, mad=mad, sample.size=sample.size), function(fun){ apply(meanldrsubs, 1, fun, na.rm=T) }))
	meandev = 1.96*ldstats$sd/sqrt(ldstats$sample.size) ; ldstats$lower.mean.IC0.95 = ldstats$mean - meandev ; ldstats$upper.mean.IC0.95 = ldstats$mean + meandev
	mediandev = 1.96*1.4826*ldstats$mad/sqrt(length(sims)) ; ldstats$lower.median.IC0.95 = ldstats$median - mediandev ; ldstats$upper.median.IC0.95 = ldstats$median + meandev
	colnames(ldstats) = paste(LDmetric, colnames(ldstats), sep='.')
	print(summary(ldstats))
	return(cbind(ldrollsubtables[[1]][,mergecols], ldstats))
})
ldmetable = merge(ldmetables[[1]], ldmetables[[2]], by=mergecols, all=T, sort=F)
write.table(ldmetable, file=opt$outtab, row.names=F, col.names=T, quote=F, sep='\t')

### genomic map plot of stat summary
pdf(file=opt$outpdf, width=15, height=10,
 title=sprintf("%s - local LD scan - windows %dbp", opt$datasetname, opt$windowsize))
for (average in c('mean', 'median')){
	par(mar=c(8,8,8,8))
	plot(ldmetable$reference.position, ldmetable[,sprintf('r2.%s', average)],
	 col='white',
	 ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", opt$windowsize), xlab="reference genome coordinates",
	 main=sprintf("%s (%d combined simulations) - LD scan with fixed-size windows", opt$datasetname, length(sims))
	)
	abline(h=1:10*(1/10), col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
	CIbounds = lapply(c('lower', 'upper'), function(bound){ 
		CIbound = ldmetable[,sprintf('r2.%s.%s.IC0.95', bound, average)] 
		lines(ldmetable$reference.position, CIbound, col='slategrey')
		return(CIbound)
	})
#~ 	CIboundnonas.i = lapply(CIbounds, function(x){ which(!is.na(x)[seq(1, length(x), 10)]) }) # subsample the points for graphic representation
#~ 	polygon(x=c(ldmetable$reference.position[CIboundnonas.i[[1]]], rev(ldmetable$reference.position[CIboundnonas.i[[2]]])),
#~ 	 y=c(CIbounds[[1]][CIboundnonas.i[[1]]], rev(CIbounds[[2]][CIboundnonas.i[[2]]])), col='slategrey', border=NA)
	signifpoints = ldmetable[,sprintf('Fisher.%s', average)] > opt$signifthresh
	lines(ldmetable$reference.position, ldmetable[,sprintf('r2.%s', average)], col='blue')
	points(ldmetable$reference.position[signifpoints], ldmetable[signifpoints, sprintf('r2.%s', average)], pch=20, col='purple')
	legtext = c(sprintf("%s of %d simulations", average, length(sims)),
	            "(when significantly higher than the the genomic basal level)", sprintf("95%% confidence interval around %s", average))
	legcol = c('blue', 'purple', 'slategrey')
	leglwd = c(1, 3, 1)
	legend('topleft', legend=legtext, col=legcol, lwd=leglwd, bg='white')
	
	if (!is.null(opt$expdir)){
		# additional data to be compared to the simulations
		ldrollsub = loadLDrollTables(opt$expdir, maxsize=opt$maxsize)
		# only plot those with enough SNP density
		maxdenspos = which(ldrollsub$reportsnpdens >= opt$maxsize)
		# corresponding positions in the simulations
		cormaxdenspos = sapply(ldrollsub$reference.position[maxdenspos], function(x){ 
			i = which(ldmetable$reference.position==x); if (length(i)==1){ return(i) }else{ return(NA) } })
		# plot small dots/open rounds when in/out the range of simulated data; colour them black/red when non-sig/significant
		plot(x=ldrollsub$reference.position[maxdenspos], y=ldrollsub$meanldrsub[maxdenspos],
		 ylab=sprintf("LD strength (r^2) in\n%dbp-wide windows", opt$windowsize), xlab="reference genome coordinates",
		  main=sprintf("%s - LD scan with fixed-size windows", basename(opt$expdir)), col='white')
		abline(h=1:10*(1/10), col=ifelse((1:10)%%5==0, 'grey', 'lightgrey'))
		points(x=ldrollsub$reference.position[maxdenspos], y=ldrollsub$meanldrsub[maxdenspos],
		 col=ifelse(ldrollsub$logcompldfisub[maxdenspos] > opt$signifthresh, 'red', 'grey'),
		 pch=ifelse((ldrollsub$meanldrsub[maxdenspos] < CIbounds[[1]][cormaxdenspos] | ldrollsub$meanldrsub[maxdenspos] > CIbounds[[2]][cormaxdenspos]), 1, 46),
		 cex=0.5)
		legtext = c(sprintf("out of 95%% confidence interval around simul %s", average),
		             sprintf("within 95%% confidence interval around simul %s", average),
		             "significantly higher than the genomic basal level",
		             "not significantly higher")
		legpch = c(1, 46, 20, 20)
		legcol = c('black', 'black', 'red', 'grey')
	}
	legend('topleft', legend=legtext, pch=legpch, col=legcol, bg='white')
}
dev.off()
