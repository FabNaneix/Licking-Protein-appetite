# -*- coding: utf-8 -*-
"""
Data extraction and analyses for Quinine Protein Preference QPP

Created by Fabien Naneix
"""

import numpy as np
import GF_Naneix as GF

## Metadata folder (check pathway for each exp)
metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile.txt"
medfolder="C:\\Users\\fabie\\Documents\\QPP data\\Behaviour\\"

## Extraction of data and creation of lists for each column from the metafile
def Extract (metafile):
    f = open(metafile, 'r')
    f.seek(0)
    Rows = f.readlines()[1:]
    #open the metafile and go on column 0 row 1
    
    tablerows = []
    for i in Rows:
        items = i.split('\t')
        tablerows.append(items)
    
    Medfile, Rat, Session, Box, Diet, Date, Left, Right \
    = [],[],[],[],[],[],[],[]
    
    for i, lst in enumerate(tablerows):
         Medfile = Medfile + [lst[1]]
         Rat = Rat + [lst[2]]
         Session = Session + [lst[3]]
         Box = Box + [lst[4]]
         Diet = Diet + [lst[5]]
         Date = Date + [lst[6]]
         Left = Left + [lst[7]] ## content of the left bottle
         Right = Right + [lst[8]]  ## content of the right bottle
                       
    return ({'Medfile':Medfile, 'Rat':Rat, 'Session':Session, 'Box':Box, \
             'Diet':Diet, 'Date':Date, 'Left':Left, 'Right':Right})
    #create a dict for each column

Data = Extract(metafile)

## Extraction licks ('b' = licks left; 'e'= licks right) from medfile
Licksallrats = []
Licks = []
for Medfile in Data['Medfile']:
    Licks = []
    Licks.append(GF.medfilereader(medfolder+Medfile,varsToExtract = ['b','e'], remove_var_header = True))
    Licksallrats.append(Licks)

TotalLicks = []
for l in Licks:
    TotalLicks.append(len(l[0]) + len(l[1]))