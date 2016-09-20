#! /usr/bin/env python

import sys
import re
import os

Usage = """
RELICmatch-parse.py - version 1.1

This program converts RELIC-match output to non-interleaved format
and prints out the list of aligned peptides with their start and end positions.

Usage:
RELICmatch-parse.py <path_to_the_text_file_with_RELIC_match_alignment> [option]

Option:

<path_to_the_output_file_with_non-interleaved_alignment.txt>

If optional output path is provided, the program will save non-interleaved RELIC-match alignment to the text file.

If you have any issues or questions on this program, please contact me:
vinni(at)hawaii.edu

Link to RELIC-match software:
http://www.northeastern.edu/xray/downloads/match-program/

Kirill Vinnikov 

"""

if len(sys.argv) <2 or sys.argv[1]=="-help" or sys.argv[1]=="help" or sys.argv[1]=="-h":
    print Usage
elif len(sys.argv) > 3:
    print "The number of arguments cannot be more than 2!\n"
    print Usage
else:
    inpath = sys.argv[1]
    if len(sys.argv) == 3:
        outpath = sys.argv[2]
    else:
        outpath = False

    ### list of variables
    seq_count =0   # counts protein fragments
    depth = 0      # counts rows of aligned peptides
    fragment = 0   # counts end position for each protein fragment
    AB = []        # intermediate list of peptides per fragment
    protein = []   # list of protein fragments
    align = {}     # dictionary of protein fragments with peptide matches
    FINAL = []     # final list of all peptides

    ### open file and read all lines
    try: 
        infile = open(inpath, 'r')
        all_lines = infile.readlines()
    except IOError:
        print "\nPlese check the path for your input file. I can't find it.\n"
        sys.exit()
    

    ### filter lines and count peptides for each protein fragment
    for i in range(len(all_lines)):
        line = all_lines[i]
        if line.strip().endswith("0"):
            if len(AB) > 0:
                align[protein[seq_count-1]] = AB
            seq_count +=1
            protein.append(line.strip().split()[0])
            fragment = fragment + len(protein[seq_count-1])
            depth = 0 
            AB = []
        elif seq_count > 0:
            start = fragment - len(protein[seq_count-1])
            depth +=1
            ABs = line.strip().split()

            ABs_uniq = []
            for k in ABs:
                if k not in ABs_uniq:
                    ABs_uniq.append(k)
        
            for x in ABs_uniq:

                if len(re.findall("\s"+x+"\s", line)) < 1:
                    print "ERROR in the match procedure. It is probably the program bug. Please contact me immediately."
                    print x
                    print line 
                    sys.exit()
                else:
                    match = re.finditer("\s"+x+"\s",line)
                    for m in match:
                        AB.append([x, depth, start+m.start()+1,start+m.start()+len(x)]) 
                                       
        else:
            if line != " \n":
                print "Please check your input file. \nIt should be exact output text file with RELIC-match alignment."
                sys.exit()

    align[protein[seq_count-1]] = AB        
        
    infile.close()

    ### lines for output file
    lines = []

    ### constructing total protein sequence
    total_seq = "".join(protein)
    lines.append(total_seq+"\n")

    ### finding max depth
    max_depth = 0
    for fragment in protein:
        for AB in align[fragment]:
            if AB[1] > max_depth:
                max_depth = AB[1]

    ### Adding peptides to lines
    for line in range(max_depth):
        lines.append(" "*len(total_seq)+"\n")
    for fragment in protein:
        for AB in align[fragment]:
            lines[AB[1]] = lines[AB[1]][0:AB[2]-1]+AB[0]+lines[AB[1]][AB[3]:]

    ### printing the parsed output
    for i in range(1,len(lines)):
        line=lines[i].strip().split()

        for w in line:
            if len(w)%12 != 0:
                print "ERROR in peptide substring procedure. It is probably the program bug. Please contact me immediately."
                print w
                sys.exit()
            a = len(w)/12
            if a > 1:
                for j in range(a):
                    AB = w[j*12:12*(j+1)]
                    line.append(AB)

        line_uniq = []
        for l in line:
            if len(l)==12 and l not in line_uniq:
                line_uniq.append(l)

        for x in line_uniq:
            match = re.finditer(x, lines[i])
            for m in match:
                FINAL.append([x,m.start()+1,m.start()+12])

    ### writing optional output file
    if outpath:
        if len(outpath.split("/")) == 1:
            outpath = os.getcwd() +"/"+ outpath
        outfile = open(outpath,"w")        
        for line in lines:
            outfile.write(line)
        outfile.close()

    for x in FINAL:
        print x[0]+"\t"+str(x[1])+"\t"+str(x[2])
