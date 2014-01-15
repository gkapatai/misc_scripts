#!/usr/bin/env Python
from argparse import (ArgumentParser, FileType)
import os, glob, re
import sys

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from itertools import groupby
from operator import itemgetter
from collections import Counter

# ============================================================================ #
# PROGRAM USAGE                                                                #
# ============================================================================ #

_ustr = """

        pileup2fasta.py

	%s [options]
	
	Input files: - mpileup
		     - reference.fa
	
	Output: <outdir>/consensus/sample.fa
	
	Example:
		Output in current directory:
		%s --input <path-to-file>/ERR00000.mpileup --reference <path-to.referencefile>/reference.fa --prefix ERR00000
		Define a different output directory
		%s --input <path-to-file>/ERR00000.mpileup --reference <path-to.referencefile>/reference.fa --outdir <output directory> --prefix ERR00000
	""" % (sys.argv[0], sys.argv[0], sys.argv[0])
	
def Parse_Commandline():
	"Parse the input arguments, use '-h' for help"
	
	global _args

	parser = ArgumentParser(usage=_ustr, description='pileup2fasta.py')
	parser.add_argument(
		'--input', type=str, required=True, 
		help='Input mpileup file')
	parser.add_argument(
		'--reference', type=str, required=True, 
		help='Input reference fasta file')
	parser.add_argument('--outdir', type=str, required=False, default = os.getcwd(),
		help='Output root directory')
	parser.add_argument('--prefix', type=str, required=True, 
		help='Prefix for output files')
	parser.add_argument(
		'--verbose', type=int, default=0, required=False, help='Verbosity')
	
	if (len(sys.argv) <= 1):
		parser.print_help()
		sys.exit()

	try:
	    _args = parser.parse_args()
	except:
	        sys.stderr.write('Error: Failed to parse command-line arguments')
	        sys.exit(1)
# ============================================================================ #
# Functions - Process pileup file                                              #
# ============================================================================ #

def read_pileup(pileup_file,ref_seq):
	with open(pileup_file) as pileup:
		
		hash_alignment = {}

		# Split all lines in the pileup by whitespace
		pileup_split = ( x.split() for x in pileup )
		# Group the split lines based on the first field (allele) 
		for allele, lines in groupby(pileup_split, itemgetter(0)):
			hash_alignment[allele] = []
			ref_pos = 0
			for fields in lines:
				alt_bps = {}
				nuc_num = int(fields[1]) # Actual position in ref allele
				nuc_ref = fields[2]
				orig_depth = int(fields[3])
				if nuc_num != ref_pos+1:
					for i in range(ref_pos, nuc_num-1):
						hash_alignment[allele].append([i+1, ref_seq[allele][i], 0, 0, 0, 0, 'None', alt_bps, []])
					ref_pos = nuc_num-1			# filter reads based on Phred scores - cutoff Q20 and calculate matches and mismatches
				if orig_depth != 0:
					orig_match, orig_mismatch, nuc_match,nuc_mismatch, total_indels, alt_bps, report_insertions= pileup_extract_information(nuc_ref, fields[4], fields[5])
					nuc_depth = nuc_match + nuc_mismatch
					if orig_match+orig_mismatch != orig_depth or nuc_depth > orig_depth:
						print "Attention required!"
						print "Line: {0}".format(fields)
					elif nuc_num != ref_pos+1:
						for i in range(ref_pos, nuc_num):
							hash_alignment[allele].append([i+1, ref_seq[allele][i], 0, 0, 0, 0, 'None', alt_bps, []])
						hash_alignment[allele].append([nuc_num, nuc_ref, orig_depth, nuc_match, nuc_mismatch, nuc_depth, total_indels, alt_bps, report_insertions])
						ref_pos = nuc_num
					else:
						# Hash for later processing in R
						hash_alignment[allele].append([nuc_num, nuc_ref, orig_depth, nuc_match, nuc_mismatch, nuc_depth, total_indels, alt_bps, report_insertions])
						ref_pos += 1
				else: 
					hash_alignment[allele].append([nuc_num, nuc_ref, 0, 0, 0, 0, 'None', alt_bps, []])
					ref_pos += 1
	return hash_alignment

def pileup_extract_information(ref_bp, align_bps, qualities):
	match = 0
	mismatch = 0
	probabilities = {}
	
	# remove all indels of format -1a and +1a since they do not have corresponding qualities and they refer to the following bp.
	# for example deletion  cat in positions 242-4 first appears in pos 241 with the actual deletions appearing as * in
	# the respective positions
	#ndh     241     C       79      ,$,,,,,,,,,,,,,,,,,,,,,-3catG,-3cat,,-3cat,,,,,,,,,,,,..,.,....,.........,.,.......,,...,,..... ;FFFFHGHIIHIJIIJJHII!)!*G)!!I!GGCD!!CCDDDDDDDDDCDBDDDDDDDDDCDDDJGJEDDJJJD@HFFFF
	#ndh     242     C       79      ,,,,,,,,,,,,,,,,,,,,-2at*A*,*,-2at,-2at,,-2at,,,,,-2at,-2at,,..,.,....,.........,.,.......,,...,,.....^~.       BCCFFFHIHBGIGGIIGHG!!!!J!!!H!GECC!!C>DDCDCBDDDDDDDDDBDDDDDCDDDJJJBCCJGJCDHFFFFC
	#ndh     243     A       79      ,,,,,,,,,,,,,,,,,,,**C*,***,*,,,,**,,..,.,....,.........,.,.......,,...,,...... @CCDDFHHGFGHGEHIFAH!!!!G!!!G!ACCC!!AADDDDCDDDDCDDDDDBDBDC?ABDDIIJD@CJIIACHHHGFC
	#ndh     244     T       79      ,$,$,$,,,,,,,,,,,,,,,,**C*,***,*,,,,**,,..,.,....,.........,.,.......,,...,,......      @B@DDFDEBDFEEFEI?CB!!!!H!!!B!@C>C!!:>DDDC@CDDDCDDDDDDDCCCC>DDDIGJF:>JJJ>CJHHHFC
	# The script below demonstrates that
	# >>> st1 = ',$,$.$,...,...,,.,,.,.,,,,.,,...,,.,,.,,.....,.,,,..,,.,.....,,......,,...,,,,.,.,,,,.,.,,.,,.,,,....,.,,.,.,,...+2AG,,,,,,,.,'
	# >>> st2 = '>>>CDB>HD<B@FBFCBFDJHDJBJJDDDJIDJI@H<DDDDDDDDHBDDDJDDDDDDDDDBDDDDDDDDDDDDDBDDDDDD>GDIDDJBDJDB<JIJJDIB@JDJBBJH8DBB<BDDCD'
	# >>> re.findall(r'[^,\.]', st1)
	# ['$', '$', '$', '+', '2', 'A', 'G']
	# >>> len(st1) - len(re.findall(r'[^,\.]', st1))
	# 119
	# >>> len(st2)
	# 119
	if align_bps.find("$]") != -1:
		pass
	indels = re.findall(r'[\+\-]\d{1,2}[acgtnACGTN]*', align_bps)
	# it rectifies cases where the deletion is followed by a mismatch base - see above at position 241 (-3catG)
	Indels = list(set(indels))
	for index, i in enumerate(Indels):
		if len(re.search(r'[a-zA-Z]+', i).group()) != int(re.search(r'\d+', i).group()):
			dif = len(re.search(r'[a-zA-Z]+', i).group()) - int(re.search(r'\d+', i).group())
			indels = [i[:-dif] if x == i else x for x in indels]  
	for e in list(set(indels)):
		align_bps = (align_bps.replace(e,''))
	indels = [x.upper() for x in indels]
	indels_freq = Counter(indels).most_common() if indels != [] else 'None'
	# find all matches of format ^~. or .$ and remove extra symbols leaving only . and ,
	match1 = re.findall(r'\^[0-9a-zA-z\!\ "#$%&\'()\*\+,\.\-\/:;<>\?@\[\]\\\^_`\{\}\|~]{1}[,\.]{1}', align_bps)
	for e in list(set(match1)):
		align_bps = align_bps.replace(e[:-1],'')
	match2 = re.findall(r'[,\.]{1}\$', align_bps)
	for e in list(set(match2)):
		align_bps = align_bps.replace(e[-1:],'')
		
	# look for possible mismatches and remove extra symbols
	if set(','.join(align_bps)) != set([",","."]):
		mm1= re.findall(r'\^[0-9a-zA-z\!\ "#$%&\'()\*\+,\.\-\/:;<>\?@\[\]\\\^_`\{\}\|~]{1}[acgtACGT]{1}', align_bps)
		for e in list(set(mm1)):
			align_bps = align_bps.replace(e[:-1], '')
		mm2 = re.findall(r'[acgtACGT]{1}\$', align_bps)
		for e in list(set(mm2)):
			align_bps = align_bps.replace(e[-1:], '')
	
	# calculate matches and mismatches
	match = align_bps.count('.') + align_bps.count(',')
	mismatch = sum([align_bps.upper().count(x) for x in ('*', 'A', 'C', 'G', 'T', 'N')])
	
	# filter positions based on Q score  - cutoff 20
	alt_bps = {}
	filtered_mismatch = 0
	filtered_match = 0
	match_filtered = []
	# calculate filtered match and mismatch positions
	if mismatch > 0:
		mismatch_pos = {}
		for x in (r'\*{1}', r'A{1}', r'C{1}', 'G{1}', 'T{1}', 'N{1}'):
			pt = re.compile(x)
			alt = re.search(r'[\*\+ACGTN]{1}', x).group()[0]
			try:
				pt.search(align_bps.upper()).group()
			except AttributeError:
				continue
			else:
				mismatch_pos[alt] = [m.start() for m in pt.finditer(align_bps.upper())]
		mismatch_filtered = {}
		for bp in mismatch_pos.keys():
			mismatch_filtered[bp] = []
			for m in mismatch_pos[bp]:
				Q=ord(qualities[m])-33
				if Q > 20: mismatch_filtered[bp].append(m)
				
	
	match_pos = [m.start() for m in re.finditer(r'[,\.]{1}', align_bps)]
	for m in match_pos:
		Q= ord(qualities[m])-33
		if Q > 20: match_filtered.append(m)
		
	for bp in ('A', 'C', 'G', 'T', 'N', '*'):
		if bp != ref_bp:
			alt_bps[bp] = len(mismatch_filtered[bp]) if mismatch > 0 and bp in mismatch_filtered.keys() else 0
			filtered_mismatch += alt_bps[bp] 
		else:
			alt_bps[bp] = len(match_filtered)
			filtered_match = alt_bps[bp]
			
	insertions = [x for x in filter(lambda x:x[0]=="+", indels)]
	report_insertions = []
	for ins in list(set(insertions)):
		if filtered_match > 4 and insertions.count(ins) > len(qualities)/2: # accepts a read if it occurs in more than half the reads
			report_insertions.append(ins)
	
			
	return match, mismatch, filtered_match, filtered_mismatch, indels_freq, alt_bps, report_insertions

# ============================================================================ #
# Functions - SNP/Genotype calling                                             #
# ============================================================================ #

def generate_consensus(allele_info):
	
	consensus = ""
	ambiguous = {"R":["A", "G"], "Y":["C","T"], "S":["G","C"], "W":["A","T"], "K":["G","T"], "M":["A","C"], "B":["C","G","T"], "D":["A","G","T"], "H":["A","C","T"], "V":["A","C","G"]}
	
	for nuc_info in allele_info:
		nuc_num, nuc_ref, orig_depth, nuc_match, nuc_mismatch, nuc_depth, total_indels, alt_bps, report_insertions= nuc_info
		genotype = ''
                if nuc_num == 2062725:
                    pass
		# convert the number values to frequency for final output
		if orig_depth == 0:
			consensus += "N"
			continue
		for bp in alt_bps.keys():
			#if bp == top_var:
			snp_frequency = round(alt_bps[bp]/float(nuc_depth),3) if nuc_depth > 0 else 0
			alt_bps[bp] = snp_frequency
		# find the most frequent allele and derive a genotype based on the frequency values		
		sorted_x = sorted(alt_bps.iteritems(), key=itemgetter(1), reverse=True)
		top_var = sorted_x[0][0]
		top_freq = sorted_x[0][1]
		if top_freq > 0.8:
			genotype = top_var+'/'+top_var
		else:
			if top_freq >0.2 and top_freq<0.8:
				if sorted_x[1][1] != sorted_x[2][1]:
					genotype = top_var+'/'+sorted_x[1][0]
				else:
					genotype = top_var+'/'+sorted_x[1][0]+'/'+sorted_x[2][0]
		if round(sum(alt_bps.values())) != round(1) and nuc_depth != 0:
			print "Attention required at position " + nuc_num
		if report_insertions != []:
			top_indel = total_indels[0]
			#indel = max(indel_freq.iteritems(), key=itemgetter(1))[0]
			if top_indel[0].startswith('+') and top_indel[0] in report_insertions:
				ins = report_insertions[0]
				ins_freq = top_indel[1]/nuc_depth
				if genotype == nuc_ref+'/'+nuc_ref:
					#ins = re.sub(r'\d', '', ins)
					if ins_freq > 0.8:
						genotype = '+/+'
					elif ins_freq > 0.2 and ins_freq < 0.8:
						genotype = '+/'+nuc_ref
		elif genotype.find('*') != -1:
			genotype = genotype.replace('*', '-')
		# add bp to consensus
		if genotype == '':
			consensus += "N"
		elif genotype.find("-") != -1 or genotype.find("+") != -1:
			if genotype.find("-") != -1:
				if genotype.count("-") == 2:
					consensus += "-"
				else:
					bp = re.search(r'[ACGTN]', genotype).group()
					if orig_depth < 5:
						consensus += bp.lower()
					else:
						consensus += bp
			else:
				consensus += nuc_ref + report_insertions[0][2:]
		elif  len(set(genotype.split("/"))) == 1:
			bp = genotype.split("/")[0]
			if orig_depth < 5:
				consensus += bp.lower()
			else:
				consensus += bp
		else:
			for item in ambiguous.items():
				if set(genotype.split("/")) == set(item[1]):
					if orig_depth < 5:
						consensus += item[0].lower()
					else:
						consensus += item[0]
					break
	return consensus

# ============================================================================ #
# Functions - Main                                                             #
# ============================================================================ #

def main():
	outdir = os.path.join(_args.outdir, "consensus")
	if not os.path.exists(outdir): 
		os.mkdir(outdir)

	pileup_file = os.path.basename(_args.input)
	sample =os.path.splitext(os.path.basename(pileup_file))[0]
        handle  = open(_args.reference, "rU")
        ref_seqs = {}
        for record in SeqIO.parse(handle, "fasta"):
            ref_seqs[record.id] = record.seq
	#size = len(ref_seq.seq)
	hash_alignment = read_pileup(pileup_file,ref_seqs)
	records =[]
	for allele in hash_alignment.keys():
		consensus = generate_consensus(hash_alignment[allele])
		record = SeqRecord(Seq(consensus, generic_dna), id="{0}".format(sample), description = "{0}; Consensus from {1}".format(allele, pileup_file))
		records.append(record)
	SeqIO.write(records, '{outdir}/{sample}.fa'.format(outdir=outdir, sample=sample), "fasta")


if __name__ == "__main__":
	Parse_Commandline()
	main()

#pileup2fasta.py
