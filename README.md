misc_scripts
============

Scripts written just to make my life easier

1. pileup2fasta.py

Run:

python pileup2fasta.py --input <path-to-file>/sample.mpileup --reference <path-to-ref-file>/reference.fa --prefix sample

By default the program creates a consensus folder within the current working directory and places the output file (sample.fa) there. 

Users can define a different output directory by using option --outdir, i.e.

python pileup2fasta.py --input <path-to-file>/sample.mpileup --reference <path-to-ref-file>/reference.fa --outdir <path-to-outdir>/outdir --prefix sample
