#! /usr/bin/python
import os, sys, getopt
import tree2

options, args = getopt.getopt(sys.argv[1:], 'f:r:m:s:', ['rel.var.thresh=', 'min.seq.within=', 'fix.number.cuts=', 'split.nth.underscore=', 'help']) # 
dopt = dict(options)

if len(args) < 3:
	print "usage: python prune_genotypes.py /path/to/gene_tree_path_list /path/to/output_directory"

nflnfgenetrees = args[1]
dirout = args[2]
splitlastunderscore = int(dopt.get('-s', dopt.get('--split.nth.underscore', 0)))
fnc = int(dopt.get('-f', dopt.get('--rel.var.thresh', 7)))
rvt = int(dopt.get('-r', dopt.get('--min.seq.within', 2)))
msw = int(dopt.get('-m', dopt.get('--fix.number.cut', 5)))

with open(nflnfgenetrees, 'r') as flnfgenetrees:
	lngt = [line.rstrip('\n') for line in flnfgenetrees]

if not os.path.exists(dirout): os.makedirs(dirout)

nfoutgenelist = '%s/%s_gene_list'%(dirout, os.path.basename(dirout))
with open(nfoutgenelist, 'w') as foutgenelist:
	foutgenelist.write('\n'.join([os.path.basename(pngt).split('.')[0] for pngt in lngt])+'\n')

print '\nrelvarthresh: %g\tminseqwithin: %d\nfixed number of tree cuts: %d'%(rvt, msw, fnc)
outdir = '%s/list_genotypes_rvt%g_msw%d_fnc%d'%(dirout, rvt, msw, fnc)
if not os.path.exists(outdir): os.makedirs(outdir)
for pngt in lngt:
	ngt = os.path.basename(pngt)
	print ngt
	ng = ngt.split('.')[0]
	gt = tree2.Node(fic=pngt, leafNamesAsNum=True)
	gt.complete_internal_labels()
	agt = gt.prune_genotypes(relvarthresh=rvt, minseqwithin=msw, minvarwithin=1e-5, minbs=0.8, returnLabels=True, fixnbcut=fnc)
	nfoutlist = '%s/%s.geno_labels'%(outdir, ng)
	with open(nfoutlist, 'w') as foutlist:
		for geno in agt:
			if splitlastunderscore>0:
				# split the gene name at the n-th rightmost '_' and keep the left part
				genos = [g.rsplit('_', splitlastunderscore)[0] for g in geno]
			else:
				genos = geno
			foutlist.write('\t'.join(genos)+'\n')
