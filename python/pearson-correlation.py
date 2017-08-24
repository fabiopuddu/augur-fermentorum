#!/usr/bin/python
import csv
from scipy.stats.stats import pearsonr
with open('rep_adjusted_all.txt') as f:
	reader=csv.reader(f,delimiter="\t")
	measures={}
	measures['rDNA']={}
	measures['TEL']={}
	measures['CUP1']={}
	measures['M']={}
	measures['T1']={}
	measures['T2']={}
	measures['T3']={}
	measures['T4']={}
	measures['T5']={}
	measures['GWM']={}
	SDlist=list()
	for line in reader:
		#print line
		SDname=line[0]
		SDlist.append(SDname)
		measures['rDNA'][SDname]=line[1]
		measures['CUP1'][SDname]=line[2] 
		measures['M'][SDname]=line[3] 
		measures['T1'][SDname]=line[4] 
		measures['T2'][SDname]=line[5]
		measures['T3'][SDname]=line[6]
		measures['T4'][SDname]=line[7]
		measures['T5'][SDname]=line[8]
		measures['GWM'][SDname]=line[9]  
		measures['TEL'][SDname]=line[10] 


features=('rDNA', 'CUP1', 'TEL', 'M', 'T1', 'T2', 'T3', 'T4', 'T5', 'GWM')
#print measures
maxrsquare=0
print "\t{}".format("\t".join(features))
for feature1 in features:
	rsquarelist=list()
	for feature2 in features:
		x=list()
		y=list()
		for k in SDlist:
			x.append(float(measures[feature1][k]))
			y.append(float(measures[feature2][k]))
		rsquare=pearsonr(x, y)[0]**2
		if rsquare>maxrsquare and rsquare<1:
			maxrsquare=rsquare
			maxfeature1=feature1
			maxfeature2=feature2
		rsquarelist.append(str("{:.4f}".format(rsquare)))
	print "{}\t{}".format(feature1, "\t".join(rsquarelist))
print "{} {} {}".format(maxrsquare, maxfeature1, maxfeature2)

