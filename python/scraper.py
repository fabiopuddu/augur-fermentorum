#! /usr/bin/python
from Bio import Entrez
import random
import time
import os
import sys

Entrez.email = "f.puddu@gurdon.cam.ac.uk"     # Always tell NCBI who you are
aminoacids = ['', 'A', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'F', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'G']
gene='SLN1'
systematic_gene='YIL147C'

gene=sys.argv[1]
mutation=sys.argv[2]
systematic_gene=os.popen("cat "+sys.path[0]+"/all_yeast_genes.txt | grep "+gene+ " | cut -f1 | tr -d \"\n\"").read()
filename = sys.path[0]+"/orf_trans.fasta"
protein_database=open(filename, 'r')
found = 0
protein=''
for line in protein_database:
  if systematic_gene in line and not found:
    found = 1    
    continue
  if (found == 1):	
    if ('>' not in line):
      line = line.rstrip()	
      protein = protein+line
    else:
  	  found=0
#print "Protein\tREF\t#\tALT\tPMC ID\n";
AA_N = 227
protein = [protein[index] for index in range(AA_N, len(protein)) ]

search_string = '(yeast[Abstract] OR cerevisiae[Abstract]) AND '+gene+'[Body - All Words] AND '+mutation+'[Body - All Words]'
handle = Entrez.esearch(db="pmc", term=search_string)
record = Entrez.read(handle)
PMCid = record["IdList"]
if PMCid != []:
	print "\n".join(PMCid) 

# for AA in list(protein):
#   if ('*' not in AA and AA_N>1):
#     NOTFOUND = 0
#     for aa in aminoacids:
# 			search_string = '(yeast[Abstract] OR cerevisiae[Abstract]) AND '+gene+'[Body - All Words] AND '+AA+str(AA_N)+aa+'[Body - All Words]'
# 			handle = Entrez.esearch(db="pmc", term=search_string)
# 			record = Entrez.read(handle)
# 			PMCid = record["IdList"]
# 			print gene+"\t"+AA+"\t"+str(AA_N)+"\t"+aa+"\t"+ ', '.join(PMCid)
# 			time.sleep(1)
#   AA_N +=1





