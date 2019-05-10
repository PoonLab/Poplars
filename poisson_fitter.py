"""
This a Python implementation of the Poisson-Fitter 2.0 algorithm described by 
Giorgi et al. in https://doi.org/10.1186/1471-2105-11-532.

The Poisson Distibution is fit to observed pairwise HD distribution using the 
Maximum Likelihood method described in by 
Lee et al. in https://doi.org/10.1016/j.jtbi.2009.07.038
"""

import re 

import argparse
from Bio import SeqIO
from Bio import Phylo 
import scipy.stats.poisson as poisson 
import matplotlib 
import hypermut 

# c_pattern = ".CONSENSUS"
# conseq = re.compile(c_pattern)


# def parse_args(): 
# 	parser = argparse.ArgumentParser(
# 		description = 'Estimates the time since infection by analyzing homogenous DNA sequences from a single '
# 					  'HIV patient early in infection and determining the best fitting Poission distribution.'
# 	)
# 	parser.add_argument('in', help="<input> path to aligned file")
# 	parser.add_argument('amb', help="<input> remove ambiguous IUPAC codes")
# 	parser.add_argument('rmut', help="<input> mutation rate")

# 	parser.add_argument('--lscale', help="<option> path to file obtained through deep sequencing")
# 	return parser.parse_args()



# def is_fasta(handle):
# """
# Checks that input file is fasta
# """
# 	fasta = SeqIO.parse(handle, "fasta")
# 	return any(fasta)


# def eq_ln(list):
# """
# Checks that all sequences between 2 consecutive consensus sequences are of equal length
# """	
# 	for header in list:
# 		for each seq in h: 
# 			if len(seq) = len(seq) of others 
# 				return True 
# 			else: 
# 				print "seqs cannot be aligned"
# 				return False


# def dna_alphabet(list):
# 	"""
# 	Checks for DNA sequences 
# 	"""
# 	for sequence in list:
# 		if seq[i] in "AGCT-":
# 			return True 
# 		else: 
# 			print "not DNA"
# 			return False 


# def is_aligned(list):
# 	"""
# 	Checks that sequences are aligned 
# 	"""
# 	for sequence in list : 
# 		for each base in seq:
				
# 			# for every A, there is an A or - and the len(seq) are all equal
# 			if eq_ln(list) and dna_alphabet(list): 
# 				return True
# 			else: 
# 				print "not aligned"
# 				return False 


# def rm_amb(handle): 
# 	"""
# 	Remove ambiguous codes 
# 	"""
# 	match any of (U, M, R, W, S, Y. K, V, H, D, B, X, N) in the sequence (tuple[1])
# 	if matches in seq:
# 		skip seq 


# def p_fitter(infile):
# 	with open(infile) as handle:
# 		# Read in as a list of tuples [(h1, s1), (h2, s2)....]
# 		aligned = [list(s) for s in SeqIO.FastaIO.SimpleFastaParser(handle)]

# 		# Check if sequences are aligned 
# 		if is_fasta(handle) and dna_alphabet(list) and eq_ln(list) and is_aligned(list): 
# 			...

# 	if args.lscale:
# 		check that name is of the correct form 
# 		(anystr).(unique identifer).(mul) OR (anystr).(mul)
# 			search for multiplicity of sequence at the end of the header name using "_[0-9][0-9]"
# 		# Check if sequences are aligned 

# 	if args.amb:
# 		rm_amb(handle)

# 		if APOBEC positions are to be removed:
# 			if all positions are to be removed: 
# 				run hypermut with input sequences, store pvals in list
# 					find all GRD motifs and replace the G with - 
# 					save seq_corrected to new file 

# 			else:
# 				find seqs with p_value < input_threshold
# 				find GRD motifs in these seqs
# 				replace G with -
# 				save corrected sequence to new file 

# 		elif sequences with APOBEC mutations are to be removed:
# 			run hypermut for each sequence
# 			if hypermut p-value < input p-value:
# 				skip sequence  
# 				create new alignment 

# 	if args.rmut:
# 		if rmut is invalid:
# 			rmut = 2.16e-05
# 		else: use input rmut 


# if __name__ =='__main__':
# 	args = parse_args()
# 	p_fitter(args.aligned)
