#! /usr/bin/env python

import sys


Usage = """
RELICmatch-align.py - version 1.0

This program visualize RELICmatch peptide alignment on a protein sequence.
The input file must be produced by RELICmatch-parse.py

Usage:
RELICmatch-align.py <path_to_protein_file> <path_to_parsed_RELICmatch_file>

The file "RELIC_ALIGNMENT.txt" will be saved to the current folder.

If you have any issues or questions on this program, please contact me:
vinni(at)hawaii.edu

Link to RELIC-match software:
http://www.northeastern.edu/xray/downloads/match-program/

Kirill Vinnikov 

"""

if len(sys.argv) <3 or sys.argv[1]=="-help" or sys.argv[1]=="help" or sys.argv[1]=="-h":
    print Usage
elif len(sys.argv) > 3:
    print "The number of arguments cannot be more than 2!\n"
    print Usage
else:
    prot_path = sys.argv[1]
    pep_path = sys.argv[2]
    
    prot_fasta = open(prot_path, 'r')

    protein = ""
    
    for line in prot_fasta:
        line = line.strip()
        if line.startswith(">"):
            continue
        else:
            protein = protein + line
    
    prot_fasta.close()


    a = len(protein)

    lines = [protein, " "*a]

    pep_file = open(pep_path, 'r')

    for line in pep_file:
        line = line.strip()
        AB = line.split()
        start = int(AB[1])-1
        end = int(AB[2])
        check = True
        for i in range(1,len(lines)):
            if lines[i][start:end]==" "*(end-start):
                lines[i] = lines[i][:start]+AB[0]+lines[i][end:]
                if len(lines[i]) != a:
                    print "Error!!!"
                check = False
                break
        if check:
            lines.append(" "*a)
            lines[len(lines)-1] = lines[len(lines)-1][:start]+AB[0]+lines[len(lines)-1][end:] 
    pep_file.close()

    outfile = open("RELIC_ALIGNMENT.txt","w")        
    for line in lines:
        outfile.write(line+'\n')
    outfile.close()
