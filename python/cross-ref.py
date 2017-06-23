#!/usr/bin/python
import csv
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn2
import subprocess
#read the data in
results_directory='./'
gene_lists={}
for rep in ('rDNA','CUP1','TEL','T1', 'T2'):
	for condition in ('+','-'):
		with open(results_directory+str(rep+'_'+condition+'.txt'))as f:
			reader=csv.reader(f,delimiter="\t")
			lst=list()
			for line in reader:
				gene=line[0].split('_')[1]	
				lst.append(gene)
			gene_lists[str(rep+condition)]=set(lst)
	gene_lists[str(rep)]=gene_lists[str(rep)+'-'].union(gene_lists[str(rep)+'+'])
telomere_lists={}
with open('telomere_length_annotations.txt') as f:
	reader=csv.reader(f,delimiter="\t")
	plus_list=list()
	minus_list=list()
	other_list=list()
	for line in reader:
		if '#' in line[0]:
			continue
		if 'increased' in line[2]:
			plus_list.append(line[0])
		elif 'decreased' in line[2]:
			minus_list.append(line[0])
		else:
			other_list.append(line[0])
	telomere_lists['TEL+']=set(plus_list)
	telomere_lists['TEL-']=set(minus_list)
	other_list=set(other_list)
	telomere_lists['TEL']=telomere_lists['TEL+'].union(telomere_lists['TEL-']).union(other_list)

v = venn2 ( [ gene_lists['TEL'], telomere_lists['TEL']], set_labels=('our', 'published'))
plt.title("Hits Overlap")
plt.savefig("test.png")

only_us=gene_lists['TEL'].difference(telomere_lists['TEL'])	
print "\n".join(list(only_us))
