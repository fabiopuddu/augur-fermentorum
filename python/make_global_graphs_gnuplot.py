#!/usr/bin/env python
import sys
import os
import re
import subprocess
import math

def deconvolute(table, column):
    output={}
    for SDno in table.keys():
        del_key=table[SDno][14]
        output[del_key]=0
    for del_str in output.keys():
        col_sum=0
        found=0
        for SDno in table.keys():
            if table[SDno][14] == del_str:
                col_sum += float(table[SDno][column])
                found +=1
        if found >0:
            #print del_str
            #print str(col_sum)+"\t"+str(found)
            col_ave= col_sum/float(found)
            #print col_ave
        else:
            col_ave=0
        output[del_str]=col_ave

    for SDno in table.keys():
        del_key=table[SDno][14]
        table[SDno][34]=output[del_key]
    return table

def deconvolute_stdev(table, column):
    output={}
    for SDno in table.keys():
        del_key=table[SDno][14]
        output[del_key]=0
    for del_str in output.keys():
        col_sum=0
        col_sum_sq=0
        found=0
        for SDno in table.keys():
            if table[SDno][14] == del_str:
                col_sum += float(table[SDno][column])
                col_sum_sq += float(table[SDno][column])**2
                found +=1
        if found >1:
            col_SD = math.sqrt((col_sum_sq/found)- (col_sum/found)**2)
            #print col_SD;
            #print col_ave
        else:
            col_SD=0
        output[del_str]=col_SD
    return output
#######################
# Define what is being analysed

pth = os.path.dirname(sys.argv[0])
print (pth)

results_file_name={1:'rDNA',\
                   2:'CUP1',\
                   3:'mtDNA',\
                   4:'twomicron',\
                   5:'Ty1',\
                   6:'Ty2',\
                   7:'Ty3',\
                   8:'Ty4',\
                   9:'Ty5',\
                   10:'GWM',\
                   12:'Telomeres',\
                   32:'Aneuploidy',\
                   33:'CR'}


#object=int(sys.argv[1])
input_file_name=sys.argv[1]

fh=open(input_file_name, 'r')
data={}
for row in fh:
    if '#' in row:
        header = row.split("\t")
        #print header
        ncol= len(header)
        #print ncol
        continue
    row=row.rstrip()
    riga=row.split("\t")
    key=riga[0]
    data[key]={}
    for col in range(0, ncol):
        try:
            data[key][col] = riga[col]
        except:
            data[key][col] = 0
fh.close()
dec_stdev_data=dict()
os.makedirs("DeconvolvedData", exist_ok=True)
os.makedirs("Plots", exist_ok=True)
for object in (1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 32, 33):
    x=list() #x will contain the Del1111_YFG1 strain name
    y={} #Y will be nested dictionary, with keys Del1111_YFG1 and b,b2,b3,b4,b22,b12, etc and value the value of the measure
    dec_data=deconvolute(data,object)
    dec_keys=sorted(dec_data, key=lambda k: (dec_data[k][34],dec_data[k][14]))
    counter=0
    previous_strain_DELNAME=''
    with open('DeconvolvedData/'+results_file_name[object]+'.txt','w') as outfh:
        for line in dec_keys:
            rep_data=data[line]
            #print rep_data
            new_strain_DELNAME=rep_data[14]
            if new_strain_DELNAME != previous_strain_DELNAME:
                counter += 1
            for k in sorted(data[line].keys()):
                outfh.write(str(data[line][k])+"\t")
            outfh.write(str(counter)+"\n")
            previous_strain_DELNAME = new_strain_DELNAME


    x=list() #x will contain the Del1111_YFG1 strain name
    y={} #Y will be nested dictionary, with keys Del1111_YFG1 and b,b2,b3,b4,b22,b12, etc and value the value of the measure
    dec_stdev_data[object]=deconvolute_stdev(data,object)
    dec_keys_stdev=sorted(dec_stdev_data[object], key=lambda k: dec_stdev_data[object][k])
    counter=1
    with open('DeconvolvedData/'+results_file_name[object]+'.stdev.txt','w')as outfh:
        for line in dec_keys_stdev:
            outfh.write(str(counter)+"\t")
            outfh.write(str(line)+"\t")
            outfh.write(str(dec_stdev_data[object][line])+"\n")
            counter += 1;

CUPRDN=0
cuprdn=0
CUPrdn=0
cupRDN=0
with open('DeconvolvedData/rDNAvsCUP1.txt','w')as outfh:
    for line in data.keys():
        outfh.write(str(line)+"\t")
        outfh.write(data[line][14]+"\t")
        outfh.write(data[line][1]+"\t") #rDNA
        outfh.write(data[line][2]+"\t") #CUP1
        if (float(data[line][1])>135.06776 or float(data[line][1])<96.32024) and (float(data[line][2])<11.607038 or float(data[line][2])>14.438762):
            outfh.write("1\n")
        else:
            outfh.write("0\n")

with open('DeconvolvedData/rDNAvsCUP1.stdev.txt','w')as outfh:
        for line in dec_stdev_data[1].keys():
                outfh.write(str(line)+"\t")
                outfh.write(str(dec_stdev_data[1][line])+"\t") #rDNA
                outfh.write(str(dec_stdev_data[2][line])+"\n") #CUP1
               # if (float(data[line][1])>135.06776 or float(data[line][1])<96.32024) and (float(data[line][2])<11.607038 or float(data[line][2])>14.438762):
               #         outfh.write("1\n")
               # else:
               #         outfh.write("0\n")

gnuplot_command='gnuplot -e "f=\''+input_file_name+'\';p=\''+pth+'\'" '+pth+'/../gnuplot/wide-plots.gpl'
print(gnuplot_command)
os.system(gnuplot_command)
