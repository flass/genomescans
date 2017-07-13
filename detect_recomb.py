#!/usr/bin/python
#-*- coding: utf-8 -*-

"""Routines for batch calls to recombination detection programs"""

__author__ = "Florent Lassalle <florent.lassalle@univ-lyon1.fr>"
__date__ = "20 February 2014"

import sys
import os, shutil
import subprocess
import time

# local paths to programs to custom
mpirunpathlocal = "/usr/bin/mpirun"
homedir = os.environ["HOME"]
progdir = "%s/Programs"%homedir
hyphypathlocal = "%s/hyphy"%progdir
phipackpath = "%s/PhiPack"%progdir
geneconvpath = "%s/GENECONV/unix.source"%progdir
nfgconvcfg = "%s/geneconv_config_recdetect.cfg"%geneconvpath

# generic definition of paths to MPI and HyPhy software using environment variables exported from job submission script 'detect_recomb.qsub'
# defaults to relevant paths when run on UCL CS cluster 
myhyphy = os.environ.get('myhyphy', "/share/apps/genomics/hyphympi")
mympi = os.environ.get('mympi', "/opt/SUNWhpc/HPC8.2.1/gnu/bin/mpirun --mca btl_tcp_if_exclude lo,eth0")

# when these paths do not exist, defaults to the local software defined by custom variables above
if os.path.exists(mympi.split()[0]): mpirunpath = mympi
else: mpirunpath = mpirunpathlocal
if os.path.exists(myhyphy): hyphypath = myhyphy
else: hyphypath = hyphypathlocal

hyphybasepath = "%s/lib/hyphy"%hyphypath
hyphyexecpath = "%s/bin"%hyphypath
hyphycmdleader = "%s/HYPHYMPI BASEPATH=%s %s/TemplateBatchFiles"%(hyphyexecpath, hyphybasepath, hyphybasepath)

nbpermut = 1000
nbpermutgconv = nbpermut*10
cwd = os.getcwd()

# creates Geneconv configuration file if it does not exist
gconvconfig = """#GCONV_CONFIG
-Seqtype=NUCL
%% Random number seed
-Startseed=12345
%% Disable search for outer fragments
-Noouter
%% Significance threshold for permutation tests on global fragments : none
%%-Nomaxsimpval
-Maxsimglobalpval=0.05
%% Number of permutations
-Numsims=%d
%% Mismatch penalty scaling factor
-GScale=1
%% Output format
-ExpFormat -Dumptab
"""%(nbpermutgconv)
if not os.access(nfgconvcfg, os.F_OK):
	fgconvcfg = open(nfgconvcfg, 'w')
	fgconvcfg.write(gconvconfig)
	fgconvcfg.close()

def toStr(l):
	ls = []
	for o in l:
		ls.append(str(o))
	return ls
	
def Phipack(task, tasktag, silent=True):
	""" detect recombination with PhiPack's tests 
	PHI (Bruen et al., Genetics 2006),
	Neighbour Similarity Score (NSS) (Jakobsen & Easteal, Comput Appl Biosci 1996)
	and Max Chi^2 (Maynard Smith, J Mol Evol 1992)
	"""	
	phi = 'NA'
	phi_pvalue = 'NA'
	nss = 'NA'
	nss_pvalue = 'NA'
	maxchi2 = 'NA'
	maxchi2_pvalue = 'NA'
	maxchi2_bkpt = 'NA'
	cmd = '%s/Phi -o -p %d -v -f %s'%(phipackpath, nbpermut, task)
	nflog = '%s/%s_Phi.log'%(cwd, tasktag)
	flog = open(nflog, 'w')
	if not silent: print "#call:\n%s"%cmd
	subprocess.call(cmd, shell=True, stdout=flog)
	flog.close()
	flog = open(nflog, 'r')
	lines = flog.readlines()
	flog.close()
	if not silent: print "#output:\n[...]\n%s"%''.join(lines[-20:])
	for line in lines:
		# PHI
		if line.startswith('Observed:'):
			phi = line.split()[2]
		if line.startswith('PHI (Permutation):'):
			phi_pvalue = line.split()[2]
		# NSS
		if line.startswith('The Neighbour Similarity score is'):
			nss = line.rstrip('\n').split()[-1]
		if line.startswith('NSS:'):
			nss_pvalue = line.split()[1]
		# Max Chi^2
		if line.startswith('Value of maximum breakpoint is:'):
			maxchi2 = line.rstrip('\n').split()[-1]
			if maxchi2=='inf': print line
		if line.startswith('Max Chi^2:'):
			maxchi2_pvalue = line.split()[2]
		if line.startswith('Coordinates of breakpoint with all sites'):
			maxchi2_bkpt = line.split('=')[1].split(', ')[1]
	return (phi, phi_pvalue, nss, nss_pvalue, maxchi2, maxchi2_pvalue, maxchi2_bkpt)
	
def Geneconv(task, tasktag, silent=True):
	"""detect recombination with Geneconv (Sawyer et al., Mol Biol Evol 1989)
	must disable log file creation (-nolog option) to avoid segfaults
	"""
	nfgcout = '%s/%s_gconv.tab'%(cwd, tasktag)
	cmd = '%s/geneconv -Seqfile=%s -Outfile=%s -Config=%s -nolog'%(geneconvpath, task, nfgcout, nfgconvcfg)
	subprocess.call(cmd, shell=True)
	fgcout=open(nfgcout, 'r')
	nb_gifrags = 0
	min_globsim_pvalue = 'NA'
	for line in fgcout:
		if line.startswith('GI'):
			nb_gifrags += 1
			lsp = line.split('\t')
			mgp = float(lsp[2])
			if min_globsim_pvalue=='NA' or mgp < min_globsim_pvalue:
				min_globsim_pvalue = mgp
	
	return(str(nb_gifrags), str(min_globsim_pvalue))


def SBP(task, tasktag, mpiruncmd, saveintermediate=None, silent=True):
	"""detect recombination with HyPhy's SBP (Kosakovsky Pond et al., Mol. Biol. Evol 2006)
	alignement file need to be in minimal FASTA format, i.e. no space-separated comment in header lines, 
	and no dot character in sequence names (beware of ACNUC-generated CDS names)
	"""
	if saveintermediate:
		if not os.path.isdir(saveintermediate): raise ValueError, "expected a directory to save intermediate files"
	sbp_deltaAICc = 'NA'
	sbp_meansplits = 'NA'
	sbp_SHtest_pvalue = 'NA'
	# create input file
	nfin = "%s/%s_sbp_in"%(cwd, tasktag)
	fin = open(nfin, "w")
	nfout = "%s/%s_sbp_out"%(cwd, tasktag)
	sbpin = '\n'.join(["1", task, "2", "1000", "HKY85", "2", nfout])+'\n'
	fin.write(sbpin)
	fin.close()
	if not silent: print "#sbp_in:\n%s"%sbpin
	# run analysis
	cmd= "%s %s/SingleBreakpointRecomb.bf < %s" %(mpiruncmd, hyphycmdleader, nfin)
	nflog = "%s/%s_SBP.log"%(cwd, tasktag)
	flog=open(nflog, 'w')
	if not silent: print "#call:\n%s"%cmd
	subprocess.call(cmd, shell=True, stdout=flog) #, stderr=flog)
	flog.close()
	flog = open(nflog, 'r')
	lines=flog.readlines()
	flog.close()
	if not silent: print "#output:\n[...]\n%s"%''.join(lines[-16:])
	lines.reverse() # results are at the end of log file
	foundbkpt = False
	for i in range(len(lines)):
		if lines[i] == "AIC-c\n":
			if lines[i-2].startswith("Best supported breakpoint is located at position"):
				sbp_deltaAICc = lines[i-3].split(": an improvement of ")[1].split("AIC points")[0]
				foundbkpt = True
				break
			elif lines[i-2].startswith("There seems to be NO recombination in this alignment"):
				sbp_deltaAICc = '0'
				break
	if foundbkpt:
		# if breakpoints were found, test their significancy as traces of recombination
		# create GARDProcessor input file
		nfin = "%s/%s_gproc_in"%(cwd, tasktag)
		fin = open(nfin, 'w')
		gpin = '\n'.join([task, nfout+"_cAIC.splits"])
		fin.write(gpin)
		fin.close()
		if not silent: print "#gproc_in:\n%s"%gpin
		# run GARDProcessor analysis
		cmd= "%s %s/GARDProcessor.bf < %s"%(mpiruncmd, hyphycmdleader, nfin)
		nflog = "%s/%s_SBPGProc.log"%(cwd, tasktag)
		flog=open(nflog, 'w')
		if not silent: print "#call:\n%s"%cmd
		subprocess.call(cmd, shell=True, stdout=flog)
		flog.close()
		flog = open(nflog, 'r')
		lines=flog.readlines()
		flog.close()
		if not silent: print "#output:\n[...]\n%s"%''.join(lines[-14:])
		lines.reverse() # results are at the end of log file
		for i in range(len(lines)):
			if lines[i].startswith("Mean splits identify:"):
				sbp_meansplits = lines[i].strip('\n').split()[-1]
			if lines[i].startswith("Breakpoint | LHS Raw p | LHS adjusted p | RHS Raw p | RHS adjusted p"):
				sbp_SHtest_pvalue = lines[i-1].strip('\n').split()[-1]
				break
	if saveintermediate and sbp_deltaAICc!='NA':
		lnf = os.listdir(cwd)
		for nf in lnf:
			for restag in ["_sbp_out", "_multi"]:
				if nf.startswith(tasktag+restag):
					ext = nf.partition(restag)[2]
					shutil.copy("%s/%s"%(cwd, nf), "%s/%s%s"%(saveintermediate, task.split('/')[-1], ext))
					break
	return (sbp_deltaAICc, sbp_meansplits, sbp_SHtest_pvalue)
	
def GARD(task, tasktag, mpiruncmd, saveintermediate=None, silent=True):
	"""detect recombination with HyPhy's GARD (Kosakovsky Pond et al., Mol. Biol. Evol 2006)
	alignement file need to be in minimal FASTA format, i.e. no spce-separated comment in header lines, 
	and no dot character in sequence names (beware of ACNUC-generated CDS names)
	"""
	if saveintermediate:
		if not os.path.isdir(saveintermediate): raise ValueError, "expected a directory to save intermediate files"
	gard_deltaAICc = 'NA'
	gard_maxnb_bkpt = 'NA'
	gard_meansplits = 'NA'
	gard_SHtest_pvalue = 'NA'
	# create GARD input file
	nfin = "%s/%s_gard_in"%(cwd, tasktag)
	fin = open(nfin, "w")
	nfout = "%s/%s_gard_out"%(cwd, tasktag)
	#~ gardin = '\n'.join(["", task, "010010", "2", "4", nfout])
	gardin = '\n'.join([task, "010010", "2", "4", nfout])+'\n'
	fin = open("%s/%s_gard_in"%(cwd, tasktag), "w")
	fin.write(gardin)
	fin.close()
	if not silent: print "#gard_in:\n%s"%gardin
	# run GARD analysis
	cmd= "%s %s/GARD.bf < %s"%(mpiruncmd, hyphycmdleader, nfin)
	nflog = "%s/%s_GARD.log"%(cwd, tasktag)
	flog=open(nflog, 'w')
	if not silent: print "#call:\n%s"%cmd
	subprocess.call(cmd, shell=True, stdout=flog)
	flog.close()
	flog=open(nflog, 'r')
	lines=flog.readlines()
	flog.close()
	if not silent: print "#output:\n[...]\n%s"%''.join(lines[-15:])
	lines.reverse() # results are at the end of log file
	foundbkpt = False
	for i in range(len(lines)):
		if lines[i].startswith("Breakpoints    c-AIC  Delta c-AIC"):
			# finds the n-breakpoints model for which the refinement of AIC-c values is the highest
			lAICc = [float(lines[i-1].split()[1])]
			k = 0
			# while (k==0) or (lAICc[k] < lAICc[k-1]):
				# k += 1
			for k in range(i):
				s = lines[i-k-1].split()[1]
				if not silent: print s,k
				try:
					f = float(s)
				except ValueError:
					break
				else:
					if f <= lAICc[-1]:
						lAICc.append(f)
			gard_deltaAICc = str(lAICc[0] - lAICc[-1])
			maxnb_bkpt = len(lAICc) - 1
			if maxnb_bkpt: foundbkpt = True
			gard_maxnb_bkpt = str(maxnb_bkpt)
			if not silent: print "gard_maxnb_bkpt %s, gard_deltaAICc %s"%(gard_maxnb_bkpt, gard_deltaAICc)
			break
	if saveintermediate and gard_deltaAICc!='NA':
		lnf = os.listdir(cwd)
		for nf in lnf:
			for restag in ["_gard_out", "_multi"]:
				if nf.startswith(tasktag+restag):
					ext = nf.partition(restag)[2]
					shutil.copy("%s/%s"%(cwd, nf), "%s/%s%s"%(saveintermediate, task.split('/')[-1], ext if ext!='' else restag))
					break
	flog.close()
			
	if foundbkpt:
		# if breakpoints were found, test their significancy as traces of recombination
		# create GARDProcessor input file
		nfin = "%s/%s_gproc_in"%(cwd, tasktag)
		fin = open(nfin, "w")
		gpin = '\n'.join([task, nfout+"_splits"])
		fin.write(gpin)
		fin.close()
		if not silent: print "#gproc_in:\n%s"%gpin
		# run GARDProcessor analysis
		cmd= "%s %s/GARDProcessor.bf < %s"%(mpiruncmd, hyphycmdleader, nfin)
		nflog = "%s/%s_GARDGProc.log"%(cwd, tasktag)
		flog=open(nflog, 'w')
		if not silent: print "#call:\n%s"%cmd
		subprocess.call(cmd, shell=True, stdout=flog)
		flog.close()
		flog=open(nflog, 'r')
		lines=flog.readlines()
		flog.close()
		if not silent: print "#output:\n[...]\n%s"%''.join(lines[-15:])
		if float(gard_deltaAICc)>0 and gard_deltaAICc!='NA':
			lines.reverse() # results are at the end of log file
			for i in range(len(lines)):
				if lines[i].startswith("Mean splits identify:"):
					gard_meansplits = lines[i].strip('\n').split()[-1]
				if lines[i].startswith("Breakpoint | LHS Raw p | LHS adjusted p | RHS Raw p | RHS adjusted p"):
					# LHS/RHS = left/right heand side
					# actualy a Kishino-Hasegawa (KH) test, not a Shimodaera-Hasegawa (SH) test
					KHres = lines[i-1].strip('\n').split('|')
					gard_SHtest_pvalue = str(min(float(KHres[2].strip(' ')), float(KHres[4].strip(' '))))
					break
	
	return (gard_deltaAICc, gard_maxnb_bkpt, gard_meansplits, gard_SHtest_pvalue)






def main(tasklist, recdetect, fresults, nbproc, tasktag, silent=False, outputmode='w', saveHyPhyintermediate=None):
	"""runs different programs for detection of recombination in a set of DNA sequence alignments
	
	takes as arguments:
	'tasklist' (list object): 	a list of alignment file paths
	'recdetect' (str): 			a comma-separated list of program(s) to run
	'fresults' (file object):	a writeable file for tabular output
	'nbproc' (int): 			the number of processors used for running HyPhy's programs
	'tasktag' (str):			a tag for temporary file (to avoid concurential read/write access to temporary files by several independent sinstances of this program)
	"""
	
	if outputmode == 'w':
		resheader = '\t'.join(["family", "phi", "phi_pvalue", "nss", "nss_pvalue", "maxchi2", "maxchi2_pvalue", "maxchi2_bkpt", "geneconv_nb_gifrags", "geneconv_min_globsim_pvalue", "sbp_deltaAICc", "sbp_meansplits", "sbp_SHtest_pvalue", "gard_deltaAICc", "gard_maxnb_bkpt", "gard_meansplits", "gard_SHtest_pvalue", "computation_time"])
		fresults.write(resheader+'\n')
		if not silent: print "#results header:\n%s"%resheader
	elif outputmode != 'a':
		raise IOError, "'fresults' file must be in 'w' or 'a' mode" 
	mpiruncmd = "%s -np %d"%(mpirunpath, nbproc)
	
	doallrecdetectprg = ("all" in recdetect)

	n = 0
	start = time.time()

	for taskline in tasklist:
		if not silent: print taskline
		if not silent: print "\n-----\n"	
		task = taskline.rstrip('\n')
		fam = os.path.basename(task).rsplit('.', 1)[0]
		n += 1
		if not silent: print "##task %d: %s"%(n, task)
		## measure of recombination
		# sets the default values for output of recombination tests for non-tested method or missed families 
		# (e.g. Phi cannot compute the permutation test with families that do not contain enough polymorphism signal)
		phi = 'NA'
		phi_pvalue = 'NA'
		nss = 'NA'
		nss_pvalue = 'NA'
		maxchi2 = 'NA'
		maxchi2_pvalue = 'NA'
		maxchi2_bkpt = 'NA'
		geneconv_nb_gifrags = 'NA'
		geneconv_min_globsim_pvalue = 'NA'
		sbp_deltaAICc = 'NA'
		sbp_meansplits = 'NA'
		sbp_SHtest_pvalue = 'NA'
		gard_deltaAICc = 'NA'
		gard_maxnb_bkpt = 'NA'
		gard_meansplits = 'NA'
		gard_SHtest_pvalue = 'NA'
		
		
		# run analysis		
		if ("Phi" in recdetect) or doallrecdetectprg:
			## detect recombination with PhiPack's tests 
			# PHI (Bruen et al., Genetics 2006),
			# Neighbour Similarity Score (NSS) (Jakobsen & Easteal, Comput Appl Biosci 1996)
			# and Max Chi^2 (Maynard Smith, J Mol Evol 1992)
			if not silent: print "#run PhiPack's Phi, NSS, Max Chi^2"
			phi, phi_pvalue, nss, nss_pvalue, maxchi2, maxchi2_pvalue, maxchi2_bkpt = Phipack(task, tasktag, silent=silent)
					
		if ("SBP" in recdetect) or doallrecdetectprg:
			## detect recombination with HyPhy's SBP (Kosakovsky Pond et al., Mol. Biol. Evol 2006)
			if saveHyPhyintermediate:
				saveintermediate = saveHyPhyintermediate+"_SBPdetails"
				if not os.access(saveintermediate, os.F_OK): os.mkdir(saveintermediate)
			if not silent: print "#run HyPhy's SBP"
			sbp_deltaAICc, sbp_meansplits, sbp_SHtest_pvalue = SBP(task, tasktag, mpiruncmd, silent=silent, saveintermediate=saveintermediate)
			
		if ("GARD" in recdetect) or doallrecdetectprg:
			## detect recombination with HyPhy's GARD (Kosakovsky Pond et al., Mol. Biol. Evol 2006)
			if saveHyPhyintermediate:
				saveintermediate = saveHyPhyintermediate+"_GARDdetails"
				if not os.access(saveintermediate, os.F_OK): os.mkdir(saveintermediate)
			if not silent: print "#run HyPhy's GARD"
			gard_deltaAICc, gard_maxnb_bkpt, gard_meansplits, gard_SHtest_pvalue = GARD(task, tasktag, mpiruncmd, silent=silent, saveintermediate=saveintermediate)
						
		if ("GENECONV" in recdetect) or doallrecdetectprg:
			## detect recombination with Geneconv (Sawyer et al., Mol Biol Evol 1989)
			if not silent: print "#run Geneconv"
			geneconv_nb_gifrags, geneconv_min_globsim_pvalue = Geneconv(task, tasktag, silent=silent)
			

		tcomp = time.time() - start
		resline = '\t'.join([fam, phi, phi_pvalue, nss, nss_pvalue, maxchi2, maxchi2_pvalue, maxchi2_bkpt, geneconv_nb_gifrags, geneconv_min_globsim_pvalue, sbp_deltaAICc, sbp_meansplits, sbp_SHtest_pvalue, gard_deltaAICc, gard_maxnb_bkpt, gard_meansplits, gard_SHtest_pvalue, str(tcomp)])
		fresults.write(resline+'\n')
		if not silent: print "#results:\n%s\n"%resline
		if not silent: print "#task %d took %s"%(n, time.strftime("%Hh %Mm %Ss", time.gmtime(tcomp)))
		start = time.time()
		

if __name__ == "__main__":
	if len(sys.argv) < 4:
		print "Usage: python detect_recomb.py task_list program(s)(p1{Phi|SBP|GARD|GENECONV|all}[,p2[,...]]) tabular_output_file [nb_MPI_processes]"
		sys.exit(2)
	nftasklist = os.path.expanduser(sys.argv[1])
	recdetect = sys.argv[2].split(',')
	nfresults = os.path.expanduser(sys.argv[3])
	if len(sys.argv) > 4: nbproc = int(sys.argv[4])
	else: nbproc = 1
	param2print = ['nftasklist', 'recdetect', 'nfresults', 'nbproc']
#	print "nftasklist = '%s';\nrecdetect = '%s';\nnfresults = '%s';\nnbproc = '%d'"%(nftasklist, ','.join(recdetect), nfresults, nbproc)
	if 'GARD' in recdetect:	param2print.append('mpirunpath')
	if ('SBP' in recdetect) or ('GARD' in recdetect): param2print.append('hyphypath')
	print '\n'.join(["%s = %s"%(par, repr(eval(par))) for par in param2print])
	ftasklist = open(nftasklist, 'r')
	tasklist = ftasklist.readlines()
	ftasklist.close()
	tasktag = os.path.basename(nftasklist)
	fresults = open(nfresults, 'w')
	main(tasklist, recdetect, fresults, nbproc, tasktag, saveHyPhyintermediate=nfresults)
	fresults.close()
	
