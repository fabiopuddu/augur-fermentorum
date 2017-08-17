#!/usr/bin/python
import csv
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn2
import subprocess
import sys
from scipy.stats import hypergeom

def FormatOutput(input_list,x):
	return ['\t'.join(input_list[x * i: x * i + x]) for i in range(0, len(input_list)/x+1)]

#Read what we want to analyse from parameters
p1=sys.argv[1]
p2=sys.argv[2]
pop=4454
#read the data in
matplotlib.rcParams.update({'font.size': 26})
results_directory='./'
gene_lists={}
for rep in ('rDNA','CUP1','M', 'TEL','T1', 'T2', 'T3', 'T4', 'T5'):
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


intrsct=gene_lists[p1].intersection(gene_lists[p2]) 
v = venn2 ( [ gene_lists[p1], gene_lists[p2]], set_labels=(p1, p2))
v2 = venn2 ( [ gene_lists[p1], gene_lists[p2]], set_labels=('our', 'published'))
plt.title("Hits Overlap")
plt.savefig("test.png", dpi=300)

only_us=gene_lists['TEL'].difference(telomere_lists['TEL'])	
already_known=gene_lists['rDNA+'].intersection(gene_lists['CUP1+'])     



###Calculate statistics
#M is the total number of objects, n is total number of Type I objects. 
#The random variate represents the number of Type I objects in N drawn without replacement from the total population.
#hpd = hypergeom(M, n, N)
prob = hypergeom.sf(len(intrsct)-1, pop, len(gene_lists[p1]), len(gene_lists[p2]))

#Find the number of random chance matches.
i=1
x=0
j=0
flag=0
while x<0.95:
	x=hypergeom.cdf(i-1, pop, len(gene_lists[p1]), len(gene_lists[p2]))
	if  x<0.05:
		j=i+1
	print '{} {}'.format(i, x)
	i+=1
i-=1
print "STATISTICS"
print "Total number of genes analysed: "+str(pop)
print "Number of genes in dataset "+p1+": "+str(len(gene_lists[p1]))
print "Number of genes in dataset "+p2+": "+str(len(gene_lists[p2]))
print "Number of genes in both datasets: "+str(len(intrsct))
if prob<0.05:
	print '{} {} {} {:.2E}'.format('Probabilty of having', len(intrsct) ,'or more hits by chance: ',prob)
else:
	print '{} {} {} \033[91m {:.2E} \033[97m'.format('Probabilty of having', len(intrsct) ,'or more hits by chance:',prob)
print 'Random overlap would result in {} to {} genes being present in both sets (p>0.95)'.format( j , i)
print 'Genes present in both datasets:'
print "\n".join(FormatOutput(list(intrsct),10))

