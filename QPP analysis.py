# -*- coding: utf-8 -*-
"""
Data extraction and analyses for Quinine Protein Preference QPP

@author: Fabien Naneix
"""
#import numpy as np
import GF_Naneix as GF


#Definition of functions used lated for the analyses
"""
Extracts function (from JEM_general_function)
Goes through the metafile (as .txt)
Extracts each column of interest and return a list for each
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
Uses the following arguments
date: date of the session of interest
LLicks: list of left licks with the same index than 'date'
RLicks: list of Right licks with the same index than 'date'
Returns for each session ('date') 4 lists of n animals as following:
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


"""
Micro_anal function:
- Goes throught a list of list (Quin_cond) corresponding of lick data for\
1 specific quinine day and condition (total, forced or free)
- Runs lickCalc on each list to analyse licking microstructure
- returns a list of dictionaries of 23 parameters
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

"""    
Sort_micro function:
- uses list of licks (here list of licks during forced or free trials)
- run the function Micro_anal (see above) on each list and return dictionnaries\
 of lick microstructure
- then sort these dicts in NR casein, NR maltodextrin, PR casein and\
 PR maltodextrin. 
- return a dict containing the 4 sorted analyses. You have a dict of 4 entries (diet_bottle)\
and each entry contains 12 lists (1 per quin day) containing 4 dicts (one per animal)\
with all licking parameters
"""
def Sort_micro(Quin_cond):
    Micro_lst = []
    for licks in Quin_cond:
        Micro_lst.append(Micro_anal(licks))
        
    Micro_NR_c = []
    Micro_NR_m = []
    Micro_PR_c =[]    
    Micro_PR_m = []
    
    for quin_test in Micro_lst: #for each day of quinine test
        for condition in quin_test[0:4]:
            Micro_NR_c.append(condition)        
        for condition in quin_test[4:8]:
            Micro_NR_m.append(condition)
        for condition in quin_test[8:12]:
            Micro_PR_c.append(condition)
        for condition in quin_test[12:]:
            Micro_PR_m.append(condition)

    Micro_NR_casein = []
    for rat in range(0, len(Micro_NR_c), 4):
        split = Micro_NR_c[rat:rat+4]
        Micro_NR_casein.append(split)
    Micro_NR_malto = []
    for rat in range(0, len(Micro_NR_m), 4):
        split = Micro_NR_m[rat:rat+4]
        Micro_NR_malto.append(split)
    Micro_PR_casein = []
    for rat in range(0, len(Micro_PR_c), 4):
        split = Micro_PR_c[rat:rat+4]
        Micro_PR_casein.append(split)    
    Micro_PR_malto = []
    for rat in range(0, len(Micro_PR_m), 4):
        split = Micro_PR_m[rat:rat+4]
        Micro_PR_malto.append(split)
    

    return {'Micro_NR_casein': Micro_NR_casein, 'Micro_NR_malto': Micro_NR_malto,\
        'Micro_PR_casein': Micro_PR_casein, 'Micro_PR_malto': Micro_PR_malto}


"""
Sort_micro function:
- uses dict created by Sort_micro
- extracts 6 parameters of licking patterns and split htem in different lists
- each list contain n=(animal*diet)session*bottle (here 192 in the following order\
NR-casein, NR-malto, PR-casein, PR-malto)
"""
def allratlistmaker(licking_patterns):
    bMean_all = []
    bNum_all = []
    clustMean_all = []
    clustNum_all = []
    freq_all = []
    total_all = []    

    for condition_name in licking_patterns:
        for quinine_conc in licking_patterns[condition_name]:
            for rat_dict in quinine_conc:
                bMean_all.append(rat_dict['bMean'])
                bNum_all.append(rat_dict['bNum'])
                clustMean_all.append(rat_dict['clustMean'])
                clustNum_all.append(rat_dict['clustNum'])
                freq_all.append(rat_dict['freq'])
                total_all.append(rat_dict['total'])
                
    return bMean_all, bNum_all, clustMean_all, clustNum_all, freq_all, total_all


"""
PrefCalc function
Calculate the % casein preference based on licks during free choice trials
% casein = licks casein / (licks casein + licks malto) *100
- Need the argument Quin_cond (sees ListMaker function) corresponding to \
a specific session of lick data for 1 specific quinine day and trial (total, forced or free)
- return 1 list of 2 lists: one for NR and one for PR diet groups
- each list contains the pref score for each rat of the diet group (n=4 per group)
"""
def PrefCalc(Quin_free):
    NR_C_free = []
    NR_M_free = []
    Pref_NR = []
    PR_C_free = []
    PR_M_free = []
    Pref_PR = []

    for rat in Quin_free[0]:
        a = len(rat)
        NR_C_free.append(a)
    for rat in Quin_free[1]:
        b = len(rat)
        NR_M_free.append(b)
    for index, len_licks in enumerate(NR_C_free):  
        cas = len_licks 
        malto = NR_M_free[index]
        all_free_licks = cas + malto
        Pref_NR.append(cas /all_free_licks * 100)

    for rat in Quin_free[2]:
        a = len(rat)
        PR_C_free.append(a)
    for rat in Quin_free[3]:
        b = len(rat)
        PR_M_free.append(b)
    
    for index, len_licks in enumerate(PR_C_free):  
        cas = len_licks 
        malto = PR_M_free[index]
        all_free_licks = cas + malto
        Pref_PR.append(cas /all_free_licks * 100)
        
    return ([Pref_NR, Pref_PR])


"-------------------------------------------------------------------------------"
## Metadata folder (check pathway for each exp)
#metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile.txt" #metafile with all data
metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile_pref.txt" #metafile with only preference tests
medfolder="C:\\Users\\fabie\\Documents\\QPP data\\Behaviour\\"

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


##Create different list for left and right licks for forced choices (trials 1 to 45) and free choices (45 to 65)    
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

##Casein and Malto total licks (forced & free choice trials) for each preference tests
Quin_all = []
Quin_all.append(ListsMaker('05/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('06/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('07/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('09/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('12/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('13/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('14/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('15/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('16/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('19/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('20/11/2018', LLicks, RLicks))
Quin_all.append(ListsMaker('21/11/2018', LLicks, RLicks))

##Casein and Malto forced licks for each preference tests
Quin_forced = []
Quin_forced.append(ListsMaker('05/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('06/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('07/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('09/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('12/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('13/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('14/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('15/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('16/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('19/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('20/11/2018', LLicks_forced, RLicks_forced))
Quin_forced.append(ListsMaker('21/11/2018', LLicks_forced, RLicks_forced))

##Casein and Malto free licks for each preference tests
Quin_free = []
Quin_free.append(ListsMaker('05/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('06/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('07/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('09/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('12/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('13/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('14/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('15/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('16/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('19/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('20/11/2018', LLicks_free, RLicks_free))
Quin_free.append(ListsMaker('21/11/2018', LLicks_free, RLicks_free))



##Calculate casein preference for each rat for a each quinine concentration test\
#based on licks during free choice trials
#Split the results in separate lists for NR and PR
Cas_Pref = []
for Quin in Quin_free:
    Cas_Pref.append(PrefCalc(Quin)) #calculation of casein pref

Cas_Pref_NR = []
Cas_Pref_PR = []
for pref in Cas_Pref:
    Cas_Pref_NR.append(pref[0]) #isolate casein preference for NR group
    Cas_Pref_PR.append(pref[1]) #isolate casein preference for PR group

#Add quin concentration as a key for both NR and PR groups casein preferences    
Cas_Pref_NR_quinine = {'Quin 0': Cas_Pref_NR[0:2], 'Quin 0.03':Cas_Pref_NR[2:4],\
                      'Quin 0.06': Cas_Pref_NR[4:6], 'Quin 0.1': Cas_Pref_NR[6:8],\
                      'Quin 0.5': Cas_Pref_NR[8:10], 'Quin 1.0': Cas_Pref_NR[10:]}
    
Cas_Pref_PR_quinine = {'Quin 0': Cas_Pref_PR[0:2], 'Quin 0.03':Cas_Pref_PR[2:4],
                      'Quin 0.06': Cas_Pref_PR[4:6], 'Quin 0.1': Cas_Pref_PR[6:8],
                      'Quin 0.5': Cas_Pref_PR[8:10], 'Quin 1.0': Cas_Pref_PR[10:]}
     
##Licking microstructure during forced and free choice trials
Micro_licks_forced = Sort_micro(Quin_forced)
Micro_licks_free = Sort_micro(Quin_free)

#Sort and export final results (1 list for each parameters with 192 values:\
#0-47: NR-casein; 48-95: NR-malto; 96-143: PR-casein; 144-191: PR-malto)
bMean_all_forced, bNum_all_forced, clustMean_all_forced, \
clustNum_all_forced, freq_all_forced, total_all_forced = allratlistmaker(Micro_licks_forced) #Forced trials

bMean_all_free, bNum_all_free, clustMean_all_free, \
clustNum_all_free, freq_all_free, total_all_free = allratlistmaker(Micro_licks_free) #Free trials