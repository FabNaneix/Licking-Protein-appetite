# -*- coding: utf-8 -*-
"""
Data extraction and analyses for Quinine Protein Preference QPP

Created by Fabien Naneix
"""

import numpy as np
import GF_Naneix as GF

## Metadata folder (check pathway for each exp)
#metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile.txt"
metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile_pref.txt"
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
    
    
def ListsMaker(date, LLicks, RLicks):
    PrefT_C_LicksNR = []
    PrefT_M_LicksNR = []
    PrefT_C_LicksPR = []
    PrefT_M_LicksPR = []
    for index, session in enumerate(Data['Date']): #index the position in the Date list
        if Data['Date'][index] == date: #select the quinine concentration 0 (no quinine)
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



## Extraction licks ('b' = licks left; 'e'= licks right) from medfile
Data = Extract(metafile)
Licksallrats = []
for Medfile in Data['Medfile']:
    Licksallrats.append(GF.medfilereader(medfolder+Medfile, varsToExtract = ['b', 'e', 'j'], remove_var_header = True))

TotalLicks = []
for l in Licksallrats:
    TotalLicks.append(len(l[0]) + len(l[1]))

LLicks = [] #left licks
RLicks = [] #right licks
Trials = [] #time stamp of 65 trials (45 forced + 20 free) Licksallrats list 3

for session in Licksallrats:
    LLicks.append(session[0])
    RLicks.append(session[1])
    Trials.append(session[2])


#Create different list for lefyt and right licks for forced choices (trials 1 to 45) and free choices (45 to 65)    
LLicks_forced = []
LLicks_free = []
RLicks_forced = []
RLicks_free = []

for index, session in enumerate(LLicks): #separate forced and free trials on left licks
    if len(Trials[index]) > 45:
        End_forced = Trials[index][45]
        Temp_licks = []
        Temp_licks2 = []
        for lick in session:
            if lick < End_forced:
                Temp_licks.append(lick)
            if lick > End_forced:
                Temp_licks2.append(lick)
        LLicks_forced.append(Temp_licks)
        LLicks_free.append(Temp_licks2)
        
for index, session in enumerate(RLicks): #separate forced and free trials on right licks
    if len(Trials[index]) > 45:
        End_forced = Trials[index][45]
        Temp_licks = []
        Temp_licks2 = []
        for lick in session:
            if lick < End_forced:
                Temp_licks.append(lick)
            if lick > End_forced:
                Temp_licks2.append(lick)
        RLicks_forced.append(Temp_licks)
        RLicks_free.append(Temp_licks2)
   
"""
Date of Preference tests with quinine.
2 tests per concentration. Date in format dd/mm/yyyy
Quin 0 mM 05/11 and 06/11/2018
Quin 0.03 mM 07/11 and 09/11/2018
Quin 0.06 mM 12/11 and 13/11/2018
Quin 0.1 mM 14/11 and 15/11/2018
Quin 0.5 mM 16/11 and 19/11/2018
Quin 1.0 mM 20/11 and 21/11/2018
"""

#Casein and Malto total licks (forced & free choice trials) for each preference tests
Quin0_t1 = ListsMaker('05/11/2018', LLicks, RLicks)
Quin0_t2 = ListsMaker('06/11/2018', LLicks, RLicks)
Quin0_03_t1 = ListsMaker('07/11/2018', LLicks, RLicks)
Quin0_03_t2 = ListsMaker('09/11/2018', LLicks, RLicks)
Quin0_06_t1 = ListsMaker('12/11/2018', LLicks, RLicks)
Quin0_06_t2 = ListsMaker('13/11/2018', LLicks, RLicks)
Quin0_1_t1 = ListsMaker('14/11/2018', LLicks, RLicks)
Quin0_1_t2 = ListsMaker('15/11/2018', LLicks, RLicks)
Quin0_5_t1 = ListsMaker('16/11/2018', LLicks, RLicks)
Quin0_5_t2 = ListsMaker('19/11/2018', LLicks, RLicks)
Quin1_0_t1 = ListsMaker('20/11/2018', LLicks, RLicks)
Quin1_0_t2 = ListsMaker('21/11/2018', LLicks, RLicks)

#Casein and Malto forced licks for each preference tests
Quin0_t1_forced = ListsMaker('05/11/2018', LLicks_forced, RLicks_forced)
Quin0_t2_forced = ListsMaker('06/11/2018', LLicks_forced, RLicks_forced)
Quin0_03_t1_forced = ListsMaker('07/11/2018', LLicks_forced, RLicks_forced)
Quin0_03_t2_forced = ListsMaker('09/11/2018', LLicks_forced, RLicks_forced)
Quin0_06_t1_forced = ListsMaker('12/11/2018', LLicks_forced, RLicks_forced)
Quin0_06_t2_forced = ListsMaker('13/11/2018', LLicks_forced, RLicks_forced)
Quin0_1_t1_forced = ListsMaker('14/11/2018', LLicks_forced, RLicks_forced)
Quin0_1_t2_forced = ListsMaker('15/11/2018', LLicks_forced, RLicks_forced)
Quin0_5_t1_forced = ListsMaker('16/11/2018', LLicks_forced, RLicks_forced)
Quin0_5_t2_forced = ListsMaker('19/11/2018', LLicks_forced, RLicks_forced)
Quin1_0_t1_forced = ListsMaker('20/11/2018', LLicks_forced, RLicks_forced)
Quin1_0_t2_forced = ListsMaker('21/11/2018', LLicks_forced, RLicks_forced)

#Casein and Malto free licks for each preference tests
Quin0_t1_free = ListsMaker('05/11/2018', LLicks_free, RLicks_free)
Quin0_t2_free = ListsMaker('06/11/2018', LLicks_free, RLicks_free)
Quin0_03_t1_free = ListsMaker('07/11/2018', LLicks_free, RLicks_free)
Quin0_03_t2_free = ListsMaker('09/11/2018', LLicks_free, RLicks_free)
Quin0_06_t1_free = ListsMaker('12/11/2018', LLicks_free, RLicks_free)
Quin0_06_t2_free = ListsMaker('13/11/2018', LLicks_free, RLicks_free)
Quin0_1_t1_free = ListsMaker('14/11/2018', LLicks_free, RLicks_free)
Quin0_1_t2_free = ListsMaker('15/11/2018', LLicks_free, RLicks_free)
Quin0_5_t1_free = ListsMaker('16/11/2018', LLicks_free, RLicks_free)
Quin0_5_t2_free = ListsMaker('19/11/2018', LLicks_free, RLicks_free)
Quin1_0_t1_free = ListsMaker('20/11/2018', LLicks_free, RLicks_free)
Quin1_0_t2_free = ListsMaker('21/11/2018', LLicks_free, RLicks_free)
















