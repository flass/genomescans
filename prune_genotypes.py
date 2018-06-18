#! /usr/bin/python
import os, sys
import tree2

if len(sys.argv) < 3:
	print "usage: python prune_genotypes.py /path/to/gene_tree_path_list /path/to/output_directory"
nflnfgenetrees = sys.argv[1]
dirout = sys.argv[2]
if len(sys.argv)>3:
	splitlastunderscore = int(sys.argv[3])
else:
	splitlastunderscore = 0
with open(nflnfgenetrees, 'r') as flnfgenetrees:
	lngt = [line.rstrip('\n') for line in flnfgenetrees]

if not os.path.exists(dirout): os.makedirs(dirout)

nfoutgenelist = '%s/%s_gene_list'%(dirout, os.path.basename(dirout))
with open(nfoutgenelist, 'w') as foutgenelist:
	foutgenelist.write('\n'.join([os.path.basename(pngt).split('.')[0] for pngt in lngt])+'\n')

fnc = 7	
rvt = 2
msw = 5
print '\nrelvarthresh: %g\tminseqwithin: %d\nfixed number of tree cuts: %d'%(rvt, msw, fnc)
outdir = '%s/list_genotypes_rvt%g_msw%d_fnc%d'%(dirout, rvt, msw, fnc)
if not os.path.exists(outdir): os.makedirs(outdir)
for pngt in lngt:
	ngt = os.path.basename(pngt)
	print ngt
	ng = ngt.split('.')[0]
	gt = tree2.Node(fic=pngt)
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
