# -*- coding: utf-8 -*-
"""
Data extraction and analyses for Quinine Protein Preference QPP

Created by Fabien Naneix
"""
import numpy as np
import GF_Naneix as GF

## Metadata folder (check pathway for each exp)
#metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile.txt" #metafile with all data
metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile_pref.txt" #metafile with only preference tests
medfolder="C:\\Users\\fabie\\Documents\\QPP data\\Behaviour\\"


"""
Extract function (from JEM_general_function)
Goes through the metafile (as .txt)
Extract each column of interest and return a list for each
"""
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
    
"""
ListMaker function
Use the following arguments
date: date of the session of interest
LLicks: list of left licks with the same index than 'date'
RLicks: list of Right licks with the same index than 'date'
Return for each session ('date') 4 lists of n animals as following:
NR_casein / NR_malto / PR_casein / NR_malto
"""    
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


#Create different list for left and right licks for forced choices (trials 1 to 45) and free choices (45 to 65)    
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


"""
Micro_anal function:
- Go throught a list of list (Quin_cond) corresponding of lick data for/
1 specific quinine day and condition (total, forced or free)
- Run lickCalc on each list to analyse licking microstructure
- return a list of dictionaries of 23 parameters
- each list is 1 rat/1 bottle in the following order:
NR_casein (4) / NR_Malto (4) / PR_casein (4) / PR_Malto (see ListsMaker function)
   
"""
def Micro_anal(Quin_cond):
    Lick_patterns_dicts = []
    for l in Quin_cond:
        for lst in l:
            x = GF.lickCalc(lst)
            Lick_patterns_dicts.append(x)
    return Lick_patterns_dicts
    
                      
#Licking microstructure for each quinine concentration for FORCED TRIALS (1 bottle available)
Micro_Quin0_t1_forced = Micro_anal(Quin0_t1_forced)
Micro_Quin0_t2_forced = Micro_anal(Quin0_t2_forced)
Micro_Quin0_03_t1_forced = Micro_anal(Quin0_03_t1_forced)
Micro_Quin0_03_t2_forced = Micro_anal(Quin0_03_t2_forced)
Micro_Quin0_06_t1_forced = Micro_anal(Quin0_06_t1_forced)
Micro_Quin0_06_t2_forced = Micro_anal(Quin0_06_t2_forced)
Micro_Quin0_1_t1_forced = Micro_anal(Quin0_1_t1_forced)
Micro_Quin0_1_t2_forced = Micro_anal(Quin0_1_t2_forced)
Micro_Quin0_5_t1_forced = Micro_anal(Quin0_5_t1_forced)
Micro_Quin0_5_t2_forced = Micro_anal(Quin0_5_t2_forced)
Micro_Quin1_0_t1_forced = Micro_anal(Quin1_0_t1_forced)
Micro_Quin1_0_t2_forced = Micro_anal(Quin1_0_t2_forced)


#Licking microstructure for each quinine concentration for FREE CHOICE TRIALS (2 bottles available)
Micro_Quin0_t1_free = Micro_anal(Quin0_t1_free)
Micro_Quin0_t2_free = Micro_anal(Quin0_t2_free)
Micro_Quin0_03_t1_free = Micro_anal(Quin0_03_t1_free)
Micro_Quin0_03_t2_free = Micro_anal(Quin0_03_t2_free)
Micro_Quin0_06_t1_free = Micro_anal(Quin0_06_t1_free)
Micro_Quin0_06_t2_free = Micro_anal(Quin0_06_t2_free)
Micro_Quin0_1_t1_free = Micro_anal(Quin0_1_t1_free)
Micro_Quin0_1_t2_free = Micro_anal(Quin0_1_t2_free)
Micro_Quin0_5_t1_free = Micro_anal(Quin0_5_t1_free)
Micro_Quin0_5_t2_free = Micro_anal(Quin0_5_t2_free)
Micro_Quin1_0_t1_free = Micro_anal(Quin1_0_t1_free)
Micro_Quin1_0_t2_free = Micro_anal(Quin1_0_t2_free)





#Transform lists of licks in Quinine conditions in array to calculate the pref with licks
NR_C_free = []
NR_M_free = []
for rat in Quin0_t1_free[0]:
    a = len(rat)
    NR_C_free.append(a)
for rat in Quin0_t1_free[1]:
    b = len(rat)
    NR_M_free.append(b)


Pref_NR = []
for index, len_licks in enumerate(NR_C_free):  
    cas = len_licks 
    malto = NR_M_free[index]
    all_free_licks = cas + malto
    Pref_NR.append(cas /all_free_licks * 100)

    
PR_C_free = []
PR_M_free = []
for rat in Quin0_t1_free[2]:
    a = len(rat)
    PR_C_free.append(a)
for rat in Quin0_t1_free[3]:
    b = len(rat)
    PR_M_free.append(b)

Pref_PR = []
for index, len_licks in enumerate(PR_C_free):  
    cas = len_licks 
    malto = PR_M_free[index]
    all_free_licks = cas + malto
    Pref_PR.append(cas /all_free_licks * 100)







