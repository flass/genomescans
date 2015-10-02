#! /usr/bin/python
import tree2
import sys, os, getopt
import numpy
from bitarray import bitarray
from itertools import product
import cPickle

##### bipartition manipulation functions

### parser funtions

def getFileNameFromPat(lnf, pat, nbmatch=None):
	l = [nf for nf in lnf if (pat in nf)]
	if nbmatch!=None and len(l)!=nbmatch: raise IndexError, "too many files (%d) match the pattern '%s' :\n%s"%(len(l), pat, str(l))
	return l
	
def codeBipartProfile(rawpat, asStr=False, s='*', d='.'):	
	"""Binary profile is coded as bit array, or as a 01 string if asStr=True 
	with the smallest clade of the two defined by the bipartition coded as '1's,
	so that bip.count() gives the size of the smallest clade.
	(in case a size equality when the global taxon set size is even, the leftmost taxon will always be coded as a '1')
	"""
	nstar = rawpat.count(s)
	ndot = rawpat.count(d)
	if nstar > ndot: dcode = {s:'0', d:'1'}
	elif nstar < ndot: dcode = {s:'1', d:'0'}
	elif rawpat[0]==s: dcode = {s:'1', d:'0'}
	else: dcode = {s:'0', d:'1'}
	prof = ''.join(dcode[char] for char in rawpat)
	if asStr: return prof
	return bitarray(prof)

def parseMrBayesTstat(nftstat, bsthresh=0, decreasing=True):
	"""returns a dict reference number of bipart for this MrBayes analysis to posterior probabilities
	
	only retruns the bipats wich proba >= bsthresh.
	if decreasing=True, assumes that the bipartitions are listed by decresing probability order 
	(as in standard MrBayes output) and stops when proba < bsthresh is reached.
	"""
	dnump = {}
	with open(nftstat, 'r') as ftstat:
		header = ftstat.readline() #discard header
		if header.startswith('[ID: '): header = ftstat.readline()
		for line in ftstat:
			lsp = line.strip('\n').split()
			num = int(lsp[0])
			p = float(lsp[2])
			if p < bsthresh:
				#~ print num, p
				if decreasing: break
				else: continue
			dnump[num] = p
	return dnump

def parseMrBayesParts(nfparts, dnump={}, decreasing=True, returnNum=False, ltax=None, dtaxorder={}, compltaxset=[], remtaxset=[], s='*', d='.', verbose=False):
	"""returns a dict of bipart binary profiles to reference number for this MrBayes analysis
	if a dictionary of posterior probabililties (PP) of biparts is provided, 
	PP are set as values of the retruned dict instead of their reference number, 
	unless unforced with returnNum=True.
	"""
	dbipnum = {} # 
	with open(nfparts, 'r') as fparts:
		nline = 0
		nID = 0
		for line in fparts:
			nline += 1
			if line.startswith('[ID: ') or line.startswith('ID'):
				nID += 1
				if nID > 2 : break
				else: continue
			lsp = line.strip('\n').split()
			num = int(lsp[0])
			inbip = lsp[1]
			bipps = [inbip+''.join(p) for p in product([s, d], repeat=len(compltaxset))]
			nbipp = 0
			for bipp in bipps:
				nbipp += 1
				if remtaxset:
					for i in remtaxset:
						bipp = bipp[:i]+bipp[i+1:]
				if not dtaxorder:
					bip = bipp
				else:
					bip = ''
					pos = dtaxorder.keys()
					pos.sort()
					for i in pos:
						bip += bipp[dtaxorder[i]]
						
				bipat = codeBipartProfile(bip, asStr=True, s=s, d=d)

				if dnump:
					if bipat.count('1')>1:
						if num not in dnump:
							continue
						p = dnump[num]
					else:
						p = 1.0
					if not returnNum: dbipnum[bipat] = p
					else: dbipnum[bipat] = num
				else:
					dbipnum[bipat] = num
	return dbipnum

#### distance functions
# uner the hypothesis that there are a (in)finite(?) number of differentiated clusters of strains / sequence types that recombine
# among each other, converting one strain sequence type into another type at a gene / genome segment,
# this happenning much more frequently than the strain within cluster diverge,
# the distance between clades (or more practically bipartitions) observed in two different gene trees can be interpreted as the
# number of single-strain recombination events that disturbed the putative common history of these gene trees concerning the focal clades 

# under this model of highly frequent recombination, observing several (consecutive) genes sharing common clades / bipartitions
# is not likely without the concourse of a force maintaining a particular genotype at a locus in a particular set of strains, 
# e.g. because of purifying selection or because of sexual isolation from (relaltively impaired recombination with) external strains. 
# The likelihood of observing common clades under a null hypothesis involning neither selection nor biased recombination is
# - decreasing with the size of the set / tracks of genes sharing a common clade, relatively to the frequence distribution of such set / tracks for all possible clades
# - decreasing with the size of the clade : the more strains, the more likely one was involed in a recombination event
# - but then potentially re-increasing with the size of the clade for large clades as they might sign long-term population structure consistent across genes
# - decreasing with the total length of the subtree (patristic distance) of the clade, i.e. with increasing time of vertical history (as recombination events erase divergent evolution) ; this favouring the hypothesis of long-term population structure ; but beware of putative very divergent strains beingbrought together systematically by Long Branch Artefact
# - decreasing with the length of the branch marking the bipartiton (relative to the total length of the subtree), indicative of sequnce conservation within the clade potentially associated with diversifying selection compared to other strains
# + increasing with the information content of a gene tree ; onem must normalize by the frequency of occurence of supported bipartitons in the focal gene trees
# + increasing with the relatedness of strains in the clade ; must normalize the frequency of observing them together by the mean / median distance among them in all gene trees (to spot groups of strain that are just part of the same very recent clone) ; caution this criterion relates to that of the total subtree length


def bitwiseBipartDist(bip1, bip2):
	"""return the distance between two bipartitions coded as binary profile as the sum of the differences in the profile
	
	take as input bit-coded biparts (bitarray object or any sort of boolean sequence recognized by bitarray()). 
	"""
	b1 = bitarray(bip1)
	b2 = bitarray(bip2)
	return min((b1 ^ b2).count(1), (b1 ^ ~b2).count(1))
	
def bitwiseCladeDist(bip1, bip2, cla1=1, cla2=1, maxUnilateral=False):
	"""compares clades and return their distance
	
	as the total number of differences (~ total recombination evolution in both genes from an ancestral state of common history) 
	or, if maxUnilateral=True, as the maximum number of differences between each clade and their common part (intersection)
	(~ maximum recombination evolution rom the putative common ancestral state)
	
	take as input bit-coded biparts (bitarray object or any sort of boolean sequence recognized by bitarray()), 
	and select in each of them the '1' or '0' clade to be compared via the 'cla1' and 'cla2' arguments.	
	"""
	if cla1: b1 = bitarray(bip1)
	else: b1 = ~bitarray(bip1)
	if cla2: b2 = bitarray(bip2)
	else: b2 = ~bitarray(bip2)
	if maxUnilateral: return (max(b1.count(1), b2.count(1)) - (b1 & b2).count(1))
	else: return (b1.count(1) + b2.count(1) - 2*(b1 & b2).count(1))	
	
def weightedBipartDist(bip1, bip2, dbipartingenes, ngene1, ngene2):
	bwd = bitwiseBipartDist(bip1, bip2)
	Pb1g1 = dbipartingenes.get(ngene1, {}).get(bip1, 0)	# bipart b1 found in consensus tree of gene g1
	Pb1g2 = dbipartingenes.get(ngene2, {}).get(bip1, 0)
	Pb2g1 = dbipartingenes.get(ngene1, {}).get(bip2, 0)
	Pb2g2 = dbipartingenes.get(ngene2, {}).get(bip2, 0)	# bipart b2 found in consensus tree of gene g2
	wb1 = (Pb1g1 / Pb2g1) # decreasing when b2 was also possibly found in g1, and if b1 was not that sure in g1
	wb2 = (Pb2g2 / Pb1g2)

#### representation functions

sizecmp = lambda x,y: len(x)-len(y)

def getBipartFromNode(node, ltax):
	taxset = node.get_leaf_labels()
	taxset.sort()
	taxset = tuple(taxset)
	antitaxset = list(set(ltax) - set(taxset))
	antitaxset.sort()
	antitaxset = tuple(antitaxset)
	bipart = [taxset, antitaxset]
	bipart.sort(cmp=sizecmp)
	return tuple(bipart)
	
def translateBipartToTaxSet(bip, ltax):
	"""translate a bipartition coded as a binary profile into sets of taxa
	
	returned coding of the bipart as a 2-tuple of clades (both tuples of string labels) ordered by increasing size,
	i.e. that bipart[0] is the smallest clade of bipart
	"""
	b = bitarray(bip)
	cla1 = [ltax[i] for i in range(len(ltax)) if b[i]]
	cla1.sort()
	cla2 = [ltax[i] for i in range(len(ltax)) if not b[i]]
	cla2.sort()
	bipart = [tuple(cla1), tuple(cla2)]
	bipart.sort(cmp=sizecmp)
	return tuple(bipart)

def translateBipartToBitaray(taxsets, ltax, asStr=False):
	"""translate a bipartition coded as sets of taxa into a binary profile
	
	assumes that the input bipart is a 2-tuple of clades (both tuples of string labels) ordered by increasing size,
	i.e. that bipart[0] is the smallest clade of bipart
	"""
	prof = ''
	for tax in ltax:
		if tax in taxsets[0]: prof += '1'
		else: prof += '0'
	if asStr:return prof
	else: return bitarray(prof)
	
	
#####

def indexBipart(bip, num, dbipartsnum, dnumbiparts, verbose=False):
	# storage as a string as two different bitarray objects will always generate a different hash key, even though they're equal bitwise 
	# string storage is more eficient than bitarray for ntaxa < 75
	if isinstance(bip, bitarray): b = bip.to01()	
	else: b = bip
	#~ if b=='010100000011101111':
		#~ print 'here', num
	n = dbipartsnum.setdefault(b, num+1)
	if verbose:
		print b
		print n, num
	if n>num: num+=1
	#~ print n
	p = dnumbiparts.setdefault(n, b)
	if p!=b:
		print 'num', num
		print 'bip', b
		print 'pib', p
		raise IndexError
	return num, n
	
def indexConsGeneBipart(bip, pp, num, dbipartsnum, dnumbiparts, ngene, dbipartingenes, dgeneswithbipart, verbose=False):
	num, n = indexBipart(bip, num, dbipartsnum, dnumbiparts, verbose=verbose)
	dbipartingenes.setdefault(ngene, {})[n] = pp
	dgeneswithbipart.setdefault(n, set()).add(ngene)
	return num	
	
def indexConsGeneTreeBipart(bip, num, dbipartsnum, dnumbiparts, ngene, cit, dbipartintrees, dtreeswithbipart, foutgttable, verbose=False):
	num, n = indexBipart(bip, num, dbipartsnum, dnumbiparts, verbose=verbose)
	dbipartintrees.setdefault(ngene, {})[n] = cit
	dtreeswithbipart.setdefault(n, set()).add(ngene)
	foutgttable.write('\t'.join([str(tic) for tic in [n]+cit])+'\n')
	return num


def main(nflngenes, dirbayesresults, bsthreshrefbip, bsthreshsampbip, outdir, bipartdistthresh, ltaxall=None, aliases='', useallgenes=False, removetrailing=False, indata=None, verbose=False):
	
	if not indata:
		dbipartsnum = {}
		dnumbiparts = {}
		dbipartingenes = {}
		dbipartintrees = {}
		dgeneswithbipart = {}
		dtreeswithbipart = {}
		num = 0
		setbip = set()
	else:
		dbipartsnum , dnumbiparts, dbipartingenes, dbipartintrees, dgeneswithbipart, dtreeswithbipart, num, ltaxall, setbip = indata
		print 'init', num, 'bipartitions', len(lngene), 'gene trees'
	
	dalias = dict([al.split('=') for al in aliases.split(',') if al])
	print 'aliases:', dalias
	
	if ltaxall:
		if len(ltaxall)<=5: print "ltaxall: %d taxa, %s"%(len(ltaxall), repr(ltaxall))
		else: print "ltaxall: %d taxa, ['%s', ...]"%(len(ltaxall), "', '".join(ltaxall[:5]))
	else:
		print "no input taxon list"
	
	lnfbayesresults = os.listdir(dirbayesresults)
	with open(nflngenes, 'r') as flngenes:
		lngenes = [line.strip(' \n') for line in flngenes]
	print len(lngenes), 'genes to evaluate'
	#~ lngenes.sort()
	
	# prepare output
	foutltax = open("%s/taxlabels"%(outdir), 'w')
	foutbiparttable = open("%s/bipart_db.tab"%(outdir), 'w')
	foutgttable = open("%s/bipart_intrees.tab"%(outdir), 'w')
	foutbipartPP = open("%s/bipart_PostProbs.tab"%(outdir), 'w')
	foutbipartclust = open("%s/bipart_clusters.tab"%(outdir), 'w')
	foutgenelist = open("%s/screened_gene_list.txt"%(outdir), 'w')

	foutgttable.write('\t'.join(['bipart', 'gene_label', 'node_label', 'smallerclade_size', 'branch_support', 'branch_length', 'rel_branch_length', 'subtree_length', 'rel_subtree_length'])+'\n')

	excludedgenes = []
	for ngene in lngenes:
		if ngene in dbipartintrees: continue
		print ngene
		lnfgenebayesresults = getFileNameFromPat(lnfbayesresults, ngene)
		# get taxon corespondence of taxonomic profile an consensus gene tree
		nfconstree = getFileNameFromPat(lnfgenebayesresults, '.con.tre', nbmatch=1)[0]
		dnexconstree = tree2.read_nexus("%s/%s"%(dirbayesresults, nfconstree), returnDict=True, allLower=False)
		genetree = dnexconstree['tree']['con_50_majrule']
		#~ ltax = dnexconstree['taxlabels']
		if removetrailing:
			### remove trailing gene name in sequence names
			# of mrbayes reference
			ltax = [lab.split('_')[0] for lab in dnexconstree['taxlabels']]
			# of gene trees
			for leaf in genetree.get_leaves():
				leaf.edit_label(leaf.label().split('_')[0])
		else:
			ltax = dnexconstree['taxlabels']
		# replace aliased strain names
		ltax = [dalias.get(tax, tax) for tax in ltax]
		if verbose: print "ltax: %s"%repr(ltax)	
		
		dtaxorder = {}
		compltaxset = []
		remtaxset = []
		if not ltaxall:
			ltaxall = ltax
		elif set(ltax) != set(ltaxall):
			if useallgenes:
				# record the positions to ignore in the bipart patterns
				suptax = set(ltax) - set(ltaxall)
				remtaxset = [ltax.index(tax) for tax in suptax]
				# completes the list of taxa with missing taxa (at the end)
				mistax = set(ltaxall) - set(ltax)
				# record the number of position to add at the end of bipart patterns
				compltaxset = range(len(ltax), len(ltax)+len(mistax))
				# add them
				ltax = ltax + list(mistax)
				
			else:
				print "different taxonomic sampling between genes:\nref: %s (%d)\n%s:%s (%d)"%(str(ltaxall), len(ltaxall), ngene, str(ltax), len(ltax))
				print 'skip gene', ngene
				excludedgenes.append(ngene)
				continue
		if ltaxall != ltax:
			try:
				dtaxorder = dict(zip(range(len(ltaxall)), [ltax.index(tax) for tax in ltaxall])) #
				invtaxorder = dict([(b,a) for a,b in dtaxorder.items()])
				if verbose: print "dtaxorder: %s"%repr(dtaxorder)
			except ValueError, e:
				print e
				print "ref: %s (%d)\n%s:%s (%d)"%(str(ltaxall), len(ltaxall), ngene, str(ltax), len(ltax))
				print 'skip gene', ngene
				excludedgenes.append(ngene)
				continue
				
		if verbose:
			print "compltaxset: %s"%repr([ltaxall[dtaxorder[i]] for i in compltaxset])
			if len(compltaxset)>0:
				print "indexes of missing taxa in bipartition profile:", [dtaxorder[i] for i in compltaxset]
				print "missing taxa are by default set to either '0' or '1' in bipartition profile"
		
		# get post. prob of biparts
		nftstat = getFileNameFromPat(lnfgenebayesresults, '.tstat', nbmatch=1)[0]
		dnump = parseMrBayesTstat("%s/%s"%(dirbayesresults, nftstat), bsthreshsampbip)
		# get binary coding of taxonomic profile representation of bipartitions
		nfparts = getFileNameFromPat(lnfgenebayesresults, '.parts', nbmatch=1)[0]
		dbippp = parseMrBayesParts("%s/%s"%(dirbayesresults, nfparts), dnump=dnump, dtaxorder=dtaxorder, ltax=ltax, compltaxset=compltaxset, remtaxset=remtaxset, verbose=verbose) #, ltax=ltax)
		
		totaltreelg = genetree.treelength(excludeSelf=True)
		genetree.complete_internal_labels(force=True, excludeLeaves=True)
		
		print len(dbippp), 'recorded bipartions'
		dbipartingenes[ngene] = dbippp
		setbip |= set(dbippp.keys())
		
		# record what bipart was set in the consensus gene tree
		for node in genetree:
			if node.is_root(): continue
			if node.is_leaf(): node.set_bs(1.0)
			if node.bs()<bsthreshrefbip: continue
			subtreelg = node.treelength()
			taxsets = getBipartFromNode(node, ltaxall)
			bipart = translateBipartToBitaray(taxsets, ltaxall, asStr=True)
			# !!!! missing taxa will by default be set to '0' by translateBipartToBitaray() 
			cit = [ngene, node.label(), bipart.count('1'), node.bs(), node.lg(), node.lg()/totaltreelg, subtreelg, subtreelg/totaltreelg]
			num = indexConsGeneTreeBipart(bipart, num, dbipartsnum, dnumbiparts, ngene, cit, dbipartintrees, dtreeswithbipart, foutgttable)

	foutgttable.close()
			
	ncgtbip = num
	print 'found', ncgtbip, 'different bipartitions with support >=', bsthreshrefbip, 'in', len(lngenes), 'consensus gene trees'
	
	# export list of taxon names
	foutltax.write('\n'.join(ltaxall)+'\n')
	foutltax.close()
	
	# restrict the gene list 
	if excludedgenes: print "excluded %d / %d genes for difference of taxonomic coverage"%(len(excludedgenes), len(lngenes)), excludedgenes
	for gene in excludedgenes:
		lngenes.remove(gene)
	
	lngenes.sort()
	foutgenelist.write('\n'.join(lngenes)+'\n')
	foutgenelist.close()
	
	# indexes non-consensus biparts
	for ngene in lngenes:
		for bip in dbipartingenes[ngene]:
			#~ num = indexConsGeneBipart(bip, dbippp[bip], num, dbipartsnum, dnumbiparts, ngene, dbipartingenes, dgeneswithbipart)
			num, numb = indexBipart(bip, num, dbipartsnum, dnumbiparts)
			dgeneswithbipart.setdefault(bip, set()).add(ngene)
			
	print 'found', num, 'different bipartitions with support >=', bsthreshsampbip, 'in bayesian tree samples of at least 1 gene'
	consnobip = set(dbipartsnum.keys()) - setbip
	print 'found', len(consnobip), 'bipartition in consensus trees not present in bayesian sample'
	print ''
	print ' '.join(ltaxall)
		

	# scores the bipart given their frequency in trees and size
	foutbiparttable.write('\t'.join(['profile', 'bipart', 'size_smallpart', 'freq_constrees', 'freq_genesPPgt%.2f'%bsthreshsampbip])+'\n')
	for numb in dnumbiparts:
		bip = dnumbiparts[numb]
		foutbiparttable.write('%s\t%d\t%d\t%d\t%d\n'%( bip, numb, bip.count('1'), len(dtreeswithbipart.get(numb, [])), len(dgeneswithbipart.get(bip, [])) ))
	
	foutbiparttable.close()
	
	# write the posterior probability profiles of biparts along genes
	foutbipartPP.write('\t'.join(lngenes)+'\n')
	for numbref in range(1, ncgtbip+1):
		foutbipartPP.write('\t'.join([str(numbref)]+["%.3f"%(dbipartingenes[ngene].get(dnumbiparts[numbref], 0)) for ngene in lngenes])+'\n')
	foutbipartPP.close()

	## degenerate matching among biparts
	bipartsim = numpy.ones((len(dtreeswithbipart), len(dgeneswithbipart)), dtype=int)*(-1)  # array initialized with -1 distance
	dbipartsimbipart = {}
	for numb in dtreeswithbipart:
		for numB in dtreeswithbipart:
			b = dnumbiparts[numb]
			B = dnumbiparts[numB]
			# mapping of similar bipartitions via comparison of the taxon content of their constituent clades
			# use maximum difference from compared clades to their intersect as a distance
			dcla = bitwiseBipartDist(b, B)
			bipartsim[numb-1, numB-1] = dcla
			# indexes bipart pairs that would potentially cluster with their distance
			if dcla <= bipartdistthresh:
				dbipartsimbipart.setdefault(numb, {})[numB] = dcla

	# search for clusters of similar biparts occuring in different gene trees
	dclusters = {}
	# for each bipart in a tree, try to find a similar bipart in each other tree
	for numbref in dnumbiparts:
		#~ print numbref
		dclusters[numbref] = {}
		for ngene in lngenes:
			mind = None
			matchbip = 'NA'
			for numbg in dbipartintrees.setdefault(ngene, {}):
				if numbref in dbipartsimbipart:
					if numbg in dbipartsimbipart[numbref]:
						if mind==None or dbipartsimbipart[numbref]<mind:
							mind = dbipartsimbipart[numbref]
							matchbip = numbg
			dclusters[numbref][ngene] = matchbip
			
	# write the fuzzy bipart profiles and store its patterns in binary form
	dpatbip = {}
	dbipclust = {}
	foutbipartclust.write('\t'.join(lngenes)+'\n')
	for numbref in range(1, ncgtbip+1):
		lbipclust = [dclusters.get(numbref, {}).get(ngene) for ngene in lngenes]
		dbipclust[numbref] = lbipclust
		bitpat = bitarray([bool(bip) for bip in lbipclust])
		foutbipartclust.write('\t'.join([str(bip) for bip in [numbref]+lbipclust])+'\n')
		dpatbip.setdefault(bitpat.to01(), []).append(numbref)
	foutbipartclust.close()

	return [dbipartsnum , dnumbiparts, dbipartingenes, dbipartintrees, dgeneswithbipart, dtreeswithbipart, num, ltaxall, setbip]

def usage():
	s =  "python bayesbipartprofile.py gene_list mrbayes_result_directory output_directory\n"
	s += "Mandatory parameters:\n"
	s += "\tgene_list\tpath to text file containing the list of gene names (one per line) ordered as in the genome.\n"
	s += "\t\tNames should match the prefix of files produced by Mrbayes run\n"
	s += "\tmrbayes_result_directory\tpath to directory where files produced by Mrbayes are located\n"
	s += "\toutput_directory\tpath to directory where the output of this script will be written\n"
	s += "Opional parameters:\n"
	s += "--aliases=[string]\tlist of aliase to be used to consider different leaf names in different gene trees to be the same taxon.\n"
	s += "\t\tProvide as \"alias1=taxon1[,alias2=taxon2,[..]\"\n"
	s += "--ltax=[path]tpath to file containing the orderred list of taxon names (one per line) to be considered as the list of all possible ones\n"
	s += "\t\t(in case of heterogeneous presence/absence of genes in taxa). Defaults to the order in the first gene tree, incremented by new taxa incountered in following gene trees.\n"
	s += "--bs.thresh.ref.bip=[float]\tminimum branch support to report a bipartition as present in the reference (consensus) tree (default 0.35)\n"
	s += "Opional parameters:\n"
	s += "--bs.thresh.samp.bip=[float]\tminimum branch support to report a bipartition from the bayesian sample (default 0.01)\n"
	s += "--bipart.dist.thresh=[integer]\tmaximum distance between bipartitions (as the number of non-shared taxa) for reporting their similarity (default 2)\n"
	s += "--remove.trail.in.sequence.names\tthe names of sequences in gene trees will be shortened (cutting at the first occurence of '_' character) to match the reference taxon list\n"
	s += "-v\t--verbose\tverbose mode\n"
	s += "-a\t--allgenes\tdo not discard genes where taxa of the reference list are missing; including genes with to much missing taxa will lead to non-ending computations as\n"
	s += "\t\tas the observed bipartitions will be extended to all the combinations of bipartitions including or excluding each of the missing taxa" 
	s += "\t\thaving at most 3 missing taxa per gene compared to the ttotal set of observed taxa appear tractable"
	return s
	
if __name__ == "__main__":

	options, args = getopt.getopt(sys.argv[1:], 'av', ['bs.thresh.ref.bip=', 'bs.thresh.samp.bip=', 'bipart.dist.thresh=', 'aliases=', 'allgenes', 'ltax=', 'verbose', 'remove.trail.in.sequence.names']) # 
	dopt = dict(options)

	if (len(args) < 3) or ('-h' in dopt) or ('--help' in dopt):
		print usage()
		exit()
		
	nflngenes = args[0]
	dirbayesresults = args[1]
	outdir = args[2]
	bsthreshrefbip = float(dopt.get('--bs.thresh.ref.bip', 0.35))
	bsthreshsampbip = float(dopt.get('--bs.thresh.samp.bip', 0.01))
	bipartdistthresh = int(dopt.get('--bipart.dist.thresh', 2))
	aliases = dopt.get('--aliases', '')
	verbose = ('-v' in dopt) or ('--verbose' in dopt)
	useallgenes = ('-a' in dopt) or ('--allgenes' in dopt)
	nfltax = dopt.get('--ltax')
	removetrailing = '--remove.trail.in.sequence.names' in dopt
	if nfltax:
		with open(nfltax, 'r') as fltax:
			ltax = [line.rstrip('\n') for line in fltax]
	else:
		ltax = None

	dall = main(nflngenes, dirbayesresults, bsthreshrefbip, bsthreshsampbip, outdir, bipartdistthresh, ltaxall=ltax, aliases=aliases, useallgenes=useallgenes, removetrailing=removetrailing, verbose=verbose)
	with open("%s/allbipartdata.pickle"%outdir, 'w') as fpickle:
		cPickle.dump(dall, fpickle, 2)
