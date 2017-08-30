#!/usr/bin/python
import csv
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn2
import subprocess
import sys
from scipy.stats import hypergeom
import pandas 
import seaborn as sns
import numpy as np
from optparse import OptionParser
#This function returns a list formed by joining every x element of the input list with tabs
def FormatOutput(input_list,x):
	return ['\t'.join(input_list[x * i: x * i + x]) for i in range(0, len(input_list)/x+1)]

def ReturnGlobalGraph():
	print "Running global analysis"
	Rows= ['rDNA-','CUP1-','M-', 'TEL-','T1-', 'T2-', 'T3-', 'T4-', 'T5-', 'rDNA+','CUP1+','M+', 'TEL+','T1+', 'T2+', 'T3+', 'T4+', 'T5+']
	Cols = Rows 
	#Initialise a matrix
	global_prob=np.zeros((len(Rows), len(Cols)))
	for n1 in range(0,18):
        	for n2 in range(0,18):
                	p1=Rows[n1]
                	p2=Cols[n2]
                	intrsct=gene_lists[p1].intersection(gene_lists[p2]) 
                	global_prob[n1,n2] = hypergeom.sf(len(intrsct)-1, pop, len(gene_lists[p1]), len(gene_lists[p2]))
                	#Set lower limit to 10 to the -9
			if global_prob[n1,n2] <  0.000000001:
                        	global_prob[n1,n2]= 0.000000001
			#When either set is empty, assign NaN
                	if len(gene_lists[p2]) == 0 or len( gene_lists[p1]) == 0 :
                        	global_prob[n1,n2]=float('NaN');
	#Mask all NaN values
	global_prob = np.ma.masked_invalid(global_prob)
	#Start plotting
	plt.clf()
	fig, ax = plt.subplots()
	heatmap=ax.pcolor(global_prob,cmap='plasma', norm=matplotlib.colors.LogNorm(vmin=0.000000001, vmax=1))
	fig.colorbar(heatmap)
	ax.patch.set(hatch='x', edgecolor='black')
	#Autoformat axes
	fig.autofmt_xdate()
	plt.yticks(np.arange(0.5, len(Rows), 1), Rows)
	plt.xticks(np.arange(0.5, len(Cols), 1), Cols)
	plt.savefig("all-against-all.png", dpi=300)
	return()

def PrintPairwiseGraph():
	plt.clf()
	matplotlib.rcParams.update({'font.size': 26})
	v = venn2 ( [ gene_lists[p1], gene_lists[p2]], set_labels=(p1, p2)) 
	plt.title("Hits Overlap")
	plt.savefig("pairwise.png", dpi=300)

#Read what we want to analyse from parameters
parser = OptionParser()
parser.add_option("-g", "--global", action="store_true", dest="glob",
                  help="Run all-against-all analys")
parser.add_option("-p", "--plot", action="store_true", dest="plot",
                  help="Plot Venn graph of the pairwise comparison")
(options, args) = parser.parse_args()
p1=args[0]
p2=args[1]
#read the list of all genes in the dataset
ALL=list()
with open ('ALL.txt') as all:
	reader=csv.reader(all,delimiter="\t")
	for line in reader:
		ALL.append(line[0])

pop=len(ALL)
ALL=set(ALL)
#read the data in
results_directory='./'
gene_lists={}
#Read in data from screening
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

#MAKE INTERSECTIONS
intrsct=gene_lists[p1].intersection(gene_lists[p2]) 
if options.plot:
        PrintPairwiseGraph()

if  (p1 == p2 and 'TEL' in p1):
	telomere_lists={}
	#read in data from previously annotated genes
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
		telomere_lists['TEL+']=telomere_lists['TEL+'].intersection(ALL)	#retain only the genes that we actually have screened
		telomere_lists['TEL-']=set(minus_list)
		telomere_lists['TEL-']=telomere_lists['TEL-'].intersection(ALL)	#retain only the genes that we actually have screened
		other_list=set(other_list)
		telomere_lists['TEL']=telomere_lists['TEL+'].union(telomere_lists['TEL-']).union(other_list)
		plt.clf()
		matplotlib.rcParams.update({'font.size': 26})
		v1 = venn2 ([ gene_lists[p1], telomere_lists[p1]], set_labels=('our', 'published'))
		plt.savefig("literature_reference.png", dpi=300)
		only_us=gene_lists[p1].difference(telomere_lists[p1])	
		#already_known=gene_lists['rDNA+'].intersection(gene_lists['CUP1+'])     
		print  "List of newly identified genes:"
		print  "\n".join(FormatOutput(list(only_us),10))

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
#	print '{} {}'.format(i, x)
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

if options.glob:
	ReturnGlobalGraph()

