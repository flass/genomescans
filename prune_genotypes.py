#! /usr/bin/python
import os, sys
import tree2

if len(sys.argv) < 3:
	print "usage: python prune_genotypes.py /path/to/gene_tree_directory /path/to/output_directory"
genetreedir = sys.argv[1]
dirout = sys.argv[2]
lngt = os.listdir(genetreedir)
if not os.path.exists(dirout): os.makedirs(dirout)

nfoutgenelist = '%s/%s_gene_list'%(dirout, os.path.basename(dirout))
with open(nfoutgenelist, 'w') as foutgenelist:
	foutgenelist.write('\n'.join([ngt.split('.')[0] for ngt in lngt]))

fnc = 7	
rvt = 2
msw = 5
print '\nrelvarthresh: %g\tminseqwithin: %d\nfixed number of tree cuts: %d'%(rvt, msw, fnc)
outdir = '%s/list_genotypes_rvt%g_msw%d_fnc%d'%(dirout, rvt, msw, fnc)
if not os.path.exists(outdir): os.makedirs(outdir)
for ngt in lngt:
	print ngt
	ng = ngt.split('.')[0]
	gt = tree2.Node(fic='%s/%s'%(genetreedir, ngt))
	gt.complete_internal_labels()
	agt = gt.prune_genotypes(relvarthresh=rvt, minseqwithin=msw, minvarwithin=1e-5, minbs=0.8, returnLabels=True, fixnbcut=fnc)
	nfoutlist = '%s/%s.geno_labels'%(outdir, ng)
	with open(nfoutlist, 'w') as foutlist:
		for geno in agt: foutlist.write('\t'.join(geno)+'\n')
