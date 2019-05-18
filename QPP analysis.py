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

#05/11/2018 to 09/11/2018 
#12/11/2018 to 16/11/2018
#19/11/2018 to 21/11/2018
    
    
    
def ListsMaker(date, LLicks, RLicks):
    PrefT_C_LicksNR = []
    PrefT_M_LicksNR = []
    PrefT_C_LicksPR = []
    PrefT_M_LicksPR = []
    for rat in np.unique(Data['Rat']): #create a list of the rats in numerical order
        for index, session in enumerate(Data['Date']): #index the position in the Date list
            if Data['Date'][index] == date: #select the quinine concentration 0 (no quinine)
                    if Data['Rat'][index] == rat:
                        if 'NR' in Data['Diet'][index]:
                            if 'Cas' in Data['Left'][index]:
                                PrefT_C_LicksNR.append(LLicks[index])
                            if 'Cas' in Data['Right'][index]:
                                PrefT_C_LicksNR.append(RLicks[index])
                            if 'Malto' in Data['Left'][index]:
                                PrefT_M_LicksNR.append(LLicks[index])
                            if 'Malto' in Data['Right'][index]:
                                PrefT_M_LicksNR.append(RLicks[index])
                        if 'PR' in Data['Diet'][index]:
                           if 'Cas' in Data['Left'][index]:
                                PrefT_C_LicksPR.append(LLicks[index])
                           if 'Cas' in Data['Right'][index]:
                                PrefT_C_LicksPR.append(RLicks[index])
                           if 'Malto' in Data['Left'][index]:
                                PrefT_M_LicksPR.append(LLicks[index])
                           if 'Malto' in Data['Right'][index]:
                                PrefT_M_LicksPR.append(RLicks[index])
                                
    return([PrefT_C_LicksNR, PrefT_M_LicksNR, PrefT_C_LicksPR, PrefT_M_LicksPR])

# If length of list is zero then go away     
                                 
                                
Data = Extract(metafile)

## Extraction licks ('b' = licks left; 'e'= licks right) from medfile
Licksallrats = []
for Medfile in Data['Medfile']:
    Licksallrats.append(GF.medfilereader(medfolder+Medfile,varsToExtract = ['b','e'], remove_var_header = True))

TotalLicks = []
for l in Licksallrats:
    TotalLicks.append(len(l[0]) + len(l[1]))

LLicks = []
RLicks = []

for session in Licksallrats:
    LLicks.append(session[0])
    RLicks.append(session[1])

quininezero_t1 = ListsMaker('05/11/2018', LLicks, RLicks)