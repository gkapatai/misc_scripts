import sys
import os
from argparse import (ArgumentParser, FileType)

# ============================================================================ #
# PROGRAM USAGE                                                                #
# ============================================================================ #

_ustr = """

       alignment.py

	%s [options]
	
	Input files: - fasta_file_1
		     - fasta_file_2
	
	Example:

	%s --input <path-to-file>/fasta1 <path-to-file>/fasta2 --module <path-to-module>/alignment.py 

	""" % (sys.argv[0], sys.argv[0])
	
def Parse_Commandline():
	"Parse the input arguments, use '-h' for help"
	
	global _args

	parser = ArgumentParser(usage=_ustr, description='pileup2fasta.py')
	parser.add_argument(
		'--input', nargs=2, type=str, required=True, 
		help='Input mpileup file')
	parser.add_argument(
		'--modulepath', type=str, required=True, 
		help='Input path to module alignment.py')
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
# IMPORT MODULES                                                               #
# ============================================================================ #

Parse_Commandline()
module = _args.modulepath
module_folder = os.path.dirname(module)
if module_folder not in sys.path:
    sys.path.append(module_folder)

import alignment
from Bio import SeqIO

# ============================================================================ #
# MAIN PROGRAM                                                                 #
# ============================================================================ #

handle1 = open(_args.input[0], "rU")
records1 = {}
for record in SeqIO.parse(handle1, "fasta"):
    records1[record.id] = record.seq

handle2 = open(_args.input[1], "rU")
records2 = {}
for record in SeqIO.parse(handle2, "fasta"):
    records2[record.id]=record.seq

for allele in records1.keys():
    print allele
    alignment.needle(records1[allele], records2[allele])
    
