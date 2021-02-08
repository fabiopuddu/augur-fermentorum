def get_metadata(experiment=0, force=False):
    import mysql.connector as mydb
    import os
    from dbconfig import credentials
    output=list()
    d=credentials()
    conn = mydb.connect(host=d['host'], user=d['user'], password=d['pass'], database=d['db'])
    mycursor = conn.cursor()
    query="""SELECT 	gDNA_samples.TubeBarcodeID,
                        `DelName`,
                        `ActualDelName`,
                        `FinalDelCheck`,
                        `SequencingPlateID`,
                        AlsoKnownAs,
                        Filename,
                        gDNA_samples.ERSNumber,
                        `Ploidy_Expected`,
                        submitted_ftp,
                        is_control,
                        fastq_ftp,
                        gDNA_samples.SampleName,
                       Ploidy_Observed
                    FROM gDNA_samples
                    LEFT JOIN  `gDNA_sequencing` ON `gDNA_sequencing`.`TubeBarcodeID` = `gDNA_samples`.`TubeBarcodeID`
                    LEFT JOIN gDNA_ENA on ERR = `run_accession`
                    WHERE Exper = """+str(experiment)+""" AND Filename IS NOT NULL
                    ORDER BY gDNA_samples.SequencingPlateID, gDNA_samples.DelName, gDNA_sequencing.TubeBarcodeID"""
    mycursor.execute(query)
    myresult = mycursor.fetchall()
    for x in myresult:
        if x[3] != "X" or force: #if the sample has not been blacklisted
            o=dict()
            o['bcID']=  x[0]
            if  x[2] == None or x[2] == "":  # If the strain is confirmed or unknown
                o['dname']= x[1]
            else:                                           # If the strain is not confirmed
                o['dname']= x[2]    #Use the ActualDelName Field when this has been set
            o['plate']= x[4]
            if x[5] != None:
                o['aka']=   x[5]
            else:
                o['aka']=   x[12]
            o['fname']= x[6]
            o['ers']=   x[7]
            o['ploidy']=x[8] if x[13]=="" or  x[13]==None else x[13]
            if x[9]  == None or x[9] == "":
                o['url']=   x[11]
            else:
                o['url']=   x[9]
            o['ctrl']=  x[10]
            o['SD']=x[12]
            o['delcheck']=x[3]
            output.append(o)
    return output

def write_nc(nameconversion,experiment_number):
    NC = get_metadata(experiment_number)
    with open(nameconversion, 'w') as name_conversion:
        for row in NC:
            line="\t".join(map(str,[ row['bcID'] , row['dname'] , row['plate'] , row['aka'] , row['fname'] , row['ers'] , row['ploidy'] , row['url'] , row['ctrl'], row['SD'], row['delcheck']]))
            name_conversion.write(line+"\n")

def read_nc(nameconversion):
    import urllib.parse
    ploidy=dict()
    delname=dict()
    ers=dict()
    sample_name=dict()
    bcode=dict()
    urls=dict()
    bcIDs=dict()
    plate_lib=dict()
    prev_TubeBarcodeID=""
    temp_l=list()
    aka=dict()
    fln_for_merge=dict()
    control=dict()
    delcheck=dict()
    with open (nameconversion) as name_conversion:
        for line in name_conversion:
            line=line.rstrip()
            TubeBarcodeID,dname,plate,gene,fname,ERSnumber,pl,url,ctrl,SD,dc = line.split("\t")
            #if (dname == cwd ):
            if ctrl != "0":
                try:
                    control[dname].append(TubeBarcodeID)
                except:
                    control[dname]=list()
                    control[dname].append(TubeBarcodeID)
            try:
                bcIDs[dname].append(TubeBarcodeID)
            except:
                bcIDs[dname]=list()
                bcIDs[dname].append(TubeBarcodeID)
            ploidy[TubeBarcodeID]=pl
            delname[TubeBarcodeID]=dname
            ers[TubeBarcodeID]=ERSnumber
            plate_lib[TubeBarcodeID]=plate
            sample_name[TubeBarcodeID]=SD
            bcode[fname]=TubeBarcodeID
            aka[TubeBarcodeID]=gene
            delcheck[TubeBarcodeID]=dc
            urls[fname]=urllib.parse.quote(url, safe="/")
            if TubeBarcodeID == prev_TubeBarcodeID:
                temp_l.append(fname)
            else:
                if len(temp_l)>0:
                    fln_for_merge[prev_TubeBarcodeID] = temp_l
                temp_l=list()
                temp_l.append(fname)
            prev_TubeBarcodeID = TubeBarcodeID
    fln_for_merge[prev_TubeBarcodeID] = temp_l
    for s in bcIDs.keys():
        bcIDs[s]=list(set(bcIDs[s]))
    return ploidy,delname,plate_lib,ers,bcode,urls,bcIDs,fln_for_merge,control,aka,sample_name,delcheck

def update_del_database(bcID,most_likely_candidate,exper):
    import mysql.connector as mydb
    from dbconfig import credentials
    d=credentials()   
    conn = mydb.connect(host=d['host'], port=d['write_port'], user=d['user'], password=d['pass'], database=d['db'])
    #THIS BIT SETS THE ACTUAL DEL NAME TO
                #"+" IF THE CORRECTED DELETION WAS DETECTED
                #"X" IF THE GENE WAS NOT DELETED AND INFORMATION ON THE DELETION CANNOT BE RETRIEVED FROM BARCODES
    if most_likely_candidate == "X" or most_likely_candidate == "+":
        mycursor = conn.cursor()
        query1="""UPDATE gDNA_samples_to_experiment
                  SET `FinalDelCheck` = '""" + most_likely_candidate + """'
                  WHERE `bcID` = '""" + bcID + """' AND Exper = """ + str(exper)
    #    print (query1)
        mycursor.execute(query1)
        #print (query1)
        conn.commit()
    else:
        #WHEN THE EXPECTED GENE IS NOT DELETED, AND WE HAVE INFO FROM THE BARCODES
        #THIS BIT QUERIES THE DATABASE TO UNDERSTAND WHAT DELETION NUMBER SHOULD BE GIVEN TO THE GENE CORRESPONDING TO THE DETECTED BARCODES
        mycursor = conn.cursor()
        query0= """SELECT DISTINCT IF(ActualDelName IS NOT NULL,
                    ActualDelName,
                    DelName) AS d_name
                    from `gDNA_samples`
                    WHERE `Exper`= """ + str(exper) + """
                    HAVING d_name LIKE '%"""+most_likely_candidate+"""'"""
        #print(query0)
        mycursor.execute(query0)
        #print (query0)
        myresult = mycursor.fetchall()
        newdelID=None
        for x in myresult:
            newdelID=x[0]
        #IF NOTHING IS FOUND, IT MEANS THAT THE GENE CORRESPONDING TO THE DETECTED BARCODES WAS NOT PRESENT IN THE EXPERIMENT TO START WITH
        #SO WE NEED TO ASSIGN IT A DELETION NUMBER CORRESPONDING TO THE HIGHEST IN THE DATABASE PLUS ONE
        if not newdelID or newdelID == 'NULL':
            mycursor = conn.cursor()
            query2=""" SELECT DISTINCT IF(ActualDelName IS NOT NULL AND ActualDelName != '+' AND ActualDelName != 'X',
			           SUBSTRING_INDEX(SUBSTRING_INDEX(ActualDelName,'_', 1), 'Del', -1),
			           SUBSTRING_INDEX(SUBSTRING_INDEX(DelName,'_', 1), 'Del', -1)) AS d_name
                       FROM `gDNA_samples`
                       WHERE exper= """+str(exper)+"""
                       ORDER BY CAST(d_name AS INT) DESC LIMIT 1
                       """
            mycursor.execute(query2)
            #print (query2)
            myresult = mycursor.fetchall()
            for x in myresult:
                highest_del_number=int(float(x[0]))
            newdelID = 'Del'+str(highest_del_number+1)+"_"+most_likely_candidate
        #WE NEED TO CHEC WHETHER THE DATABSE ALREADY CONTAINS SOMETHING
        mycursor = conn.cursor()
        query1="""SELECT ActualDelName
         	  FROM `gDNA_samples`
                  WHERE `TubeBarcodeID` = '""" + bcID + "'"
        mycursor.execute(query1)
        myresult = mycursor.fetchall()
        #print (query1)
        AdN=None
        for x in myresult:
            AdN=x[0]
        #print("ADN:",AdN)
        #AND NOW FINALLY WE NEED TO UPDATE THE DATABASE
        if AdN and AdN !='NULL':
            #IF THERE WAS SOMETHING IN THE ActualDelName THEN IT MEANS THAT WE HAVE ALREADY REASSIGNED THIS SAMPLE AND IF IT COMES UP TO HERE
            #WE NEED TO MARK IT AS PERMANENTLY FAILED
            mycursor = conn.cursor()
            query1="""UPDATE gDNA_samples_to_experiment
                      SET `FinalDelCheck` = 'X'
                      WHERE `bcID` = '""" + bcID + """' AND Exper = """+str(exper)
            #print (query1)
            mycursor.execute(query1)
            conn.commit()
        else:
            #OTHERWISE IT MEANS THAT IT'S THE FIRST TIME WE ARE LOOKING AT THIS SAMPLE
            #AND WE SHOULD GIVE IT ANOTHER OPPORTUNITY
            mycursor = conn.cursor()
            query1="""UPDATE gDNA_samples_to_experiment
                      SET `ActualDelName` = '""" + newdelID + """'
                      WHERE `bcID` = '""" + bcID + """' AND `ActualDelName` IS NULL AND Exper="""+str(exper)
            #print (query1)
            mycursor.execute(query1)
            conn.commit()



def update_repDNA_db(fields,numbers,exper):
    import mysql.connector as mydb
    from dbconfig import credentials
    f=fields.copy()
    f.append('Mat')
    f.append('Exper')
    if float(numbers[11]) <= -0.35:
        numbers.append("alpha")
    elif float(numbers[11]) >=  0.35:
        numbers.append("a")
    elif float(numbers[11]) > -0.35 and float(numbers[11]) < 0.35:
        numbers.append("a/alpha")
    else:
        exit(1)
    numbers.append(str(exper))
    f=list(map(lambda a: "`"+a+"`",f))
    numbers=list(map(lambda a: "\""+a+"\"",numbers))
    updateLine=", ".join([str(x[0])+"="+str(x[1]) for x in list(zip(f,numbers))])
    d=credentials()
    conn = mydb.connect(host=d['host'], port=d['write_port'], user=d['user'], password=d['pass'], database=d['db'])
    query="""INSERT INTO gDNA_raw_repDNA ("""+ ','.join(f) +""")
             VALUES ("""+ ','.join(numbers) + """) ON DUPLICATE KEY UPDATE """ + updateLine
    #print (query)
    mycursor = conn.cursor()
    mycursor.execute(query)
    conn.commit()
def update_suppressor_db(sample,valid_sup="", proposed_sup="", all_suppressors=""):
    import mysql.connector as mydb
    from dbconfig import credentials
    f=['TubeBarcodeID', 'ValidatedSuppress', 'PutativeSuppress', 'full_genotype']
    val=[sample,valid_sup,proposed_sup,all_suppressors]
    val=list(map(lambda x:"\""+x+"\"",val))
    updateLine=", ".join([str(x[0])+"="+str(x[1]) for x in list(zip(f,val))])
    d=credentials()
    conn = mydb.connect(host=d['host'], port=d['write_port'], user=d['user'], password=d['pass'], database=d['db'])
    query="""INSERT INTO gDNA_suppressors ("""+ ','.join(f) +""")
             VALUES ("""+ ','.join(val) + ") ON DUPLICATE KEY UPDATE " + updateLine
    #print (query)
    mycursor = conn.cursor()
    mycursor.execute(query)
    conn.commit()

def update_repl_asymDB(decile,exper,delname,header,riga):
    import mysql.connector as mydb
    from dbconfig import credentials
    h=header.copy()
    h.append("Exper")
    h.append("DelName")  
    h.append("decile")
    riga.append(exper)
    riga.append(delname)
    riga.append(decile)
    f=list(map(lambda x:"`"+x+"`",h))
    val=list(map(lambda x:"\""+str(x)+"\"",riga))
    updateLine=", ".join([str(x[0])+"="+str(x[1]) for x in list(zip(f,val))])
    d=credentials()
    conn = mydb.connect(host=d['host'], port=d['write_port'], user=d['user'], password=d['pass'], database=d['db'])
    query="""INSERT INTO gDNA_SNPS_repl_bias ("""+ ','.join(f) +""")
             VALUES ("""+ ','.join(val) + ") ON DUPLICATE KEY UPDATE " + updateLine
    #print (query)
    mycursor = conn.cursor()
    mycursor.execute(query)
    conn.commit()    

def update_mutDB(bcid,data,exper):
    import mysql.connector as mydb
    from dbconfig import credentials
    f=sorted(list(data.keys()))
    val=[data[k] for k in f]
    f.append('bcID')
    f.append('Exper')
    val.append(bcid)
    val.append(str(exper))
    f=list(map(lambda x:"`"+x+"`",f))
    val=list(map(lambda x:"\""+x+"\"",val))
    updateLine=", ".join([str(x[0])+"="+str(x[1]) for x in list(zip(f,val))])
    d=credentials()
    conn = mydb.connect(host=d['host'], port=d['write_port'], user=d['user'], password=d['pass'], database=d['db'])
    query="""INSERT INTO gDNA_mutations ("""+ ','.join(f) +""")
             VALUES ("""+ ','.join(val) + ") ON DUPLICATE KEY UPDATE " + updateLine
    #print (query)
    mycursor = conn.cursor()
    mycursor.execute(query)
    conn.commit()

def update_vcf_db(exper,sample,CHROM,POS,ID,REF,TYPE,ALT,QUAL,FILTER,CSQ,GT,DP,DV,GQ):
    import mysql.connector as mydb
    from dbconfig import credentials
    h=["Exper","TubeBarcodeID","CHROM","POS","ID","REF","TYPE","ALT","QUAL","FILTER","CSQ","GT","DP","DV","GQ"]
    riga=[exper,sample,CHROM,POS,ID,REF,TYPE,ALT,QUAL,FILTER,CSQ,GT,DP,DV,GQ]
    f=list(map(lambda x:"`"+x+"`",h))
    val=list(map(lambda x:"\""+str(x)+"\"",riga))
    updateLine=", ".join([str(x[0])+"="+str(x[1]) for x in list(zip(f,val))])
    d=credentials()
    conn = mydb.connect(host=d['host'], port=d['write_port'], user=d['user'], password=d['pass'], database=d['db'])
    query="""INSERT INTO gDNA_vcf ("""+ ','.join(f) +""")
             VALUES ("""+ ','.join(val) + ") ON DUPLICATE KEY UPDATE " + updateLine
    #print (query)
    mycursor = conn.cursor()
    mycursor.execute(query)
    conn.commit()    

def suppressor_choice(file,screening):
    from prompt_toolkit.shortcuts import radiolist_dialog
    import os
    dirpath = os.path.basename(os.getcwd())
    proceed = radiolist_dialog(
        values=[(True, "Yes"),(False, "No"),(None,"No to All")],
        title=dirpath,
        text="Is this a suppressor screen?")
    if proceed:
        field = radiolist_dialog(
            values=[("ValidatedSuppress","ValidatedSuppress"),('PutativeSuppress','PutativeSuppress')],
            title=dirpath,
            text="Validated or Putative Suppressors?")
        with open(file) as fh:
            for line in fh:
                riga=line.rstrip().split("\t")
                if len(riga) > 2:
                    sample=riga[0]
                    aka=riga[1]
                    genotype=riga[2:]
                    gen=list(map(lambda x:(x,x), genotype))
                    gen.append(('unknown','unknown'))
                    result = radiolist_dialog(
                                 values=gen,
                                 title='Screening: '+dirpath,
                                 text='Please select the suppressor mutation for sample '+aka+' ('+sample+'):')
                    if result != None:
                        if result!='unknown':
                            if field == "ValidatedSuppress":
                                update_suppressor_db(sample,valid_sup=result, all_suppressors=" ".join(genotype))
                            elif field == "PutativeSuppress":
                                update_suppressor_db(sample,proposed_sup=result, all_suppressors=" ".join(genotype))
                        else:
                            update_suppressor_db(sample, all_suppressors=" ".join(genotype))
                    else:
                        break
    elif proceed == None:
        return('NoToAll')
