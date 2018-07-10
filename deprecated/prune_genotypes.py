#! /usr/bin/python
import os, sys, getopt
import tree2

options, args = getopt.getopt(sys.argv[1:], 'g:o:f:r:m:s:M:b:hv', ['gene.tree.list=', 'output.dir=', 'rel.var.thresh=', 'min.seq.within=', 'min.var.within=', 'min.bs.stem=', 'fix.number.cuts=', 'split.nth.underscore=', 'help', 'verbose']) # 
dopt = dict(options)

if ('-h' in dopt) or ('--help' in dopt):
	print "usage: python prune_genotypes.py -g /path/to/gene_tree_path_list -o /path/to/output_directory [-r X][-m K][-f C]"

nflnfgenetrees = dopt.get('-g', dopt.get('--gene.tree.list'))
dirout  = dopt.get('-o', dopt.get('--output.dir'))
splitlastunderscore = int(dopt.get('-s', dopt.get('--split.nth.underscore', 0)))
rvt = float(dopt.get('-r', dopt.get('--rel.var.thresh', 2.0)))
msw = int(dopt.get('-m', dopt.get('--min.seq.within', 5)))
mvw = float(dopt.get('-M', dopt.get('--min.var.within', 1e-5)))
mbs = float(dopt.get('-b', dopt.get('--min.bs.stem', 0.8)))
fnc = int(dopt.get('-f', dopt.get('--fix.number.cut', 0)))
verbose = ('-v' in dopt) or ('--verbose' in dopt)

with open(nflnfgenetrees, 'r') as flnfgenetrees:
	lngt = [line.rstrip('\n') for line in flnfgenetrees]

if not os.path.exists(dirout): os.makedirs(dirout)

nfoutgenelist = '%s/%s_gene_list'%(dirout, os.path.basename(dirout))
with open(nfoutgenelist, 'w') as foutgenelist:
	foutgenelist.write('\n'.join([os.path.basename(pngt).split('.')[0] for pngt in lngt])+'\n')

if verbose: print '\nrelvarthresh: %g\tminseqwithin: %d\nfixed number of tree cuts: %d'%(rvt, msw, fnc)
outdir = '%s/list_genotypes_rvt%g_msw%d_fnc%d'%(dirout, rvt, msw, fnc)
if not os.path.exists(outdir): os.makedirs(outdir)
for pngt in lngt:
	ngt = os.path.basename(pngt)
	if verbose: print ngt
	ng = ngt.split('.')[0]
	gt = tree2.Node(fic=pngt, leafNamesAsNum=True)
	gt.complete_internal_labels()
	agt = gt.prune_genotypes(relvarthresh=rvt, minseqwithin=msw, minvarwithin=1e-5, minbs=0.8, returnLabels=True, fixnbcut=fnc, silent=(not verbose))
	nfoutlist = '%s/%s.geno_labels'%(outdir, ng)
	with open(nfoutlist, 'w') as foutlist:
		for geno in agt:
			if splitlastunderscore>0:
				# split the gene name at the n-th rightmost '_' and keep the left part
				genos = [g.rsplit('_', splitlastunderscore)[0] for g in geno]
			else:
				genos = geno
			foutlist.write('\t'.join(genos)+'\n')
