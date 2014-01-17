misc_scripts
============

Scripts written just to make my life easier

1. pileup2fasta.py

Run:

python pileup2fasta.py --input <path-to-file>/sample.mpileup --reference <path-to-ref-file>/reference.fa --prefix sample

By default the program creates a consensus folder within the current working directory and places the output file (sample.fa) there. 

Users can define a different output directory by using option --outdir, i.e.

python pileup2fasta.py --input <path-to-file>/sample.mpileup --reference <path-to-ref-file>/reference.fa --outdir <path-to-outdir>/outdir --prefix sample

2. run_alignment.py

Run

python run_alignment.py --input <path-to-file>/fasta1 <path-to-file>/fasta2 --module <path-to-module>/alignment.py

The script uses the alignment.py module, provided in the align2sequences folder. The alignment.py was adapted by version available at https://github.com/alevchuk/pairwise-alignment-in-python. This software implements the Smith-Waterman and Needleman-Wunsch algorithms to do local and global sequence analysis, respectively. 
The results for the global alignment were tested against results from
http://www.ebi.ac.uk/Tools/psa/emboss_needle/.
