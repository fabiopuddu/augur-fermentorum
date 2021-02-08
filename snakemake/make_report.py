#!/usr/bin/env python
from reportlab.lib import colors,utils
from reportlab.lib.pagesizes import A4, landscape
from reportlab.platypus import Table, TableStyle, Paragraph, Image, SimpleDocTemplate, PageBreak, PageTemplate

from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.styles import getSampleStyleSheet
import PIL
from reportlab.lib.units import cm, inch
import sys

delfolder=str(sys.argv[1])

leftborder=45
rightborder=45
lowerborder=50

def get_data(file):
    with open (file) as f:
        d=list()
        for line in f:
            line=line.rstrip()
            row=line.split("\t")
            t=list()
            for tab in row:
                t.append(tab)
            d.append(t)
    return d

def AllPageSetup(canvas, canv):

    canvas.saveState()
    #header
    #canvas.drawString(0.5*inch, 8*inch, 'test')
    canvas.setFont('Helvetica-Bold', 10)
    canvas.drawRightString(10.5*inch, 7.5* inch, delfolder)

    #footers
    canvas.setFont('Helvetica-Bold', 8)
    canvas.drawString(0.5 * inch, 0.5 * inch, "augur-fermentorum; Fabio Puddu")
    canvas.setFont('Helvetica', 8)
    canvas.drawRightString(10.5 * inch, 0.5 * inch, 'Page %d' % (canv.page))

    canvas.setFont("Helvetica", 240)
    canvas.setStrokeGray(0.90)
    canvas.setFillGray(0.90)
    #canvas.drawCentredString(5.5 * inch, 3.25 * inch, canv.watermark)

    canvas.restoreState()

def add_header(file,d):
    styles = getSampleStyleSheet()
    style = styles["BodyText"]
    if file == delfolder+'/repDNA/results.txt':
        d.insert(0,list(map(lambda x:Paragraph("<bold><font size=12>"+x+"</font></bold>",style),['BarcodeID','rDNA','CUP1','mtDNA','2µ','Ty1', 'Ty2', 'Ty3', 'Ty4', 'Ty5', 'Telomeres', 'MAT', 'GWM', 'SampleName'])))
    elif file == delfolder+'/deletion_report1':
        d.insert(0,list(map(lambda x:Paragraph("<bold><font size=12>"+x+"</font></bold>",style),['Filename','ORF checked','Gene checked','Median coverage','Percent of CWM', 'Deleted?'])))

def get_image(path, height=1*cm):
    img = utils.ImageReader(path)
    iw, ih = img.getSize()
    aspect = ih / float(iw)
    return Image(path, width=(height/aspect), height=height)

def repDNA_page():
    data=get_data(delfolder+'/repDNA/results.txt')
    controls=list()
    style = styles["Normal"]
    for i in range (len(data)):
        if '*' in data[i][0]:
            controls=data[i].copy()
    for i in range (len(data)):
        if '*' not in data[i][0]:
            for j in range(len(data[i])):
                if len(controls)>0:
                    if j>0 and j < 11 and float(controls[j]) != 0:
                        if float(data[i][j])/float(controls[j]) > 1.2 or float(data[i][j])/float(controls[j]) < 0.8:
                            data[i][j] = Paragraph("<font size=8 color=red><strong>{}</strong></font>".format(data[i][j]),style)
                        else:
                            data[i][j] = Paragraph("<font size=8>{}</font>".format(data[i][j]),style)
                    else:
                        data[i][j] =Paragraph("<font size=8>"+str(data[i][j])+"</font>",style)
                else:
                    if j>0 and j < 11:
                        data[i][j] = Paragraph("<font size=8>{}</font>".format(data[i][j]),style)
                    else:
                        data[i][j] =Paragraph("<font size=8>"+str(data[i][j])+"</font>",style)   

        else:
            for j in range(len(data[i])):
                data[i][j] = Paragraph("<font size=8 color=green><strong>{}</strong></font>".format(data[i][j]),style)
    add_header(delfolder+'/repDNA/results.txt',data)
    t = Table(data)
    t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.25, colors.black),
                           ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    data_len = len(data)

    for each in range(data_len):
        if each % 2 == 0:
            bg_color = colors.whitesmoke
        else:
            bg_color = colors.lightgrey
        t.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), bg_color)]))
    return t

def repDNA_plots():
    out=list()
    headers=("","","RibosomalDNA", "CUP1", "mtDNA", "2 micron", "Transposon Ty1",
            "Transposon Ty2", "Transposon Ty3", "Transposon Ty4", "Transposon Ty5",
            "Telomeres")
    styles = getSampleStyleSheet()
    style = styles["BodyText"]
    out.append(Paragraph("<bold><font size=18>Repetitive DNA summary</font></bold>", style))
    for x in (2,3,4,5,6,7,8,9,10,11):
        out.append(Paragraph("<bold><font size=12>"+headers[x]+"</font></bold>", style))
        out.append(get_image(delfolder+'/reports/'+str(x)+'.png',height=6*cm))
    return out

def deletion_report():
    data=get_data(delfolder+'/reports/deletion_report.txt')
    data=list(filter(lambda row: "===" not in row[0], data))  #remove headers
    data=list(filter(lambda row: "BAM" in row[0], data))  #remove headers
    add_header('deletion_report1',data)
    t = Table(data)
    data_len = len(data)
    t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.25, colors.black),
                       ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    for each in range(data_len):
        if each % 2 == 0:
            bg_color = colors.whitesmoke
        else:
            bg_color = colors.lightgrey
        if  "Deleted:N" in data[each]:
            t.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), colors.mistyrose)]))
        else:
           t.setStyle(TableStyle([('BACKGROUND', (0, each), (-1, each), bg_color)]))
    return t

def genotype_report():
    data=get_data(delfolder+'/analysis/predicted_genotypes.txt')
    controls=list()
    style = styles["Normal"]
    outdata=list()
    for i in range (len(data)):
            outdata.append(list()) #the rows
            l=list()
            m=list()
            c=1
            if len(data[i]) < 40:
                for j in range(len(data[i])):   #the columns
                    if j<2:         #first two columns
                        outdata[i].append(Paragraph("<font size=8><strong>{}</strong></font>".format(data[i][j]),style))
                    elif j>=2:      #remianing columns (genotype)
                        if "Δ" in data[i][j] or "FS" in data[i][j]:
                            l.append ( Paragraph("<font size=6 color=red><strong>{}</strong></font>".format(data[i][j]),style))
                        elif "x" in data[i][j]:
                            l.append(Paragraph("<font size=6>{}</font>".format(data[i][j]),style))
                        elif ">" in data[i][j]:
                            l.append(Paragraph("<font size=6 color=royalblue>{}</font>".format(data[i][j]),style))
                        else:
                            l.append(Paragraph("<font size=6>{}</font>".format(data[i][j]),style))

                        if c>=8 or j == len(data[i])-1:
                            m.append(l)
                            c=1
                            l=list()
                        c+=1
                if len(m) > 0:
                    outdata[i].append(Table(m))
            else:
                for j in range(2):
                    outdata[i].append(Paragraph("<font size=8><strong>{}</strong></font>".format(data[i][j]),style))
                outdata[i].append("Too many mutations")
    #add_header('repDNA/results.txt',outdata)
    t = Table(outdata)
    t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.25, colors.black),
                   ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
                   ('BACKGROUND', (0, 0), (-1, -1), colors.whitesmoke)]))
    return t

flatten = lambda l: [item for item in l if type(l) == list]
# repDNA_plots()
# aneuploidy_report_page()

styles = getSampleStyleSheet()
#print ("Styles:",styles)
style = styles["Title"]

canv = SimpleDocTemplate(delfolder+"/reports/WGS_report.pdf", pagesize=landscape(A4),rightMargin=20, leftMargin=20, topMargin=40, bottomMargin=40)

canv.report_info = "%s %s" % ('a','b')
canv.watermark = 'DRAFT'

#page_template = PageTemplate(id="fund_notes", onPage=AllPageSetup, pagesize=landscape(A4))

header1 = Paragraph("<bold><font size=18>Deletion Report</font></bold>\n", style)
header2 = Paragraph("<bold><font size=18>Repetitive DNA Report</font></bold>\n", style)
header3 = Paragraph("<bold><font size=18>Aneuploidy Report</font></bold>\n", style)
header4 = Paragraph("<bold><font size=18>Genotype Report</font></bold>\n", style)
header5 = Paragraph("<bold><font size=18>Mutation Strand Bias</font></bold>\n", style)

elements=[header1, deletion_report(),PageBreak(),\
        header4, genotype_report(),PageBreak(),\
        header2,repDNA_page(),PageBreak(),\
        *repDNA_plots(), PageBreak(),\
        header3,get_image(delfolder+'/reports/aneuploidy_report.png',height=16*cm),PageBreak(),\
        header5,get_image(delfolder+'/reports/strand_bias.png',height=16*cm),PageBreak()]


canv.build(elements, onFirstPage=AllPageSetup, onLaterPages=AllPageSetup)
