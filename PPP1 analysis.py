# -*- coding: utf-8 -*-
"""
Data extraction and analyses for Protein Preference PPP1
Preference tests 1, 2 and 3. Licking microstructure.

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
def ListsMaker(test, LLicks, RLicks):
    PrefT_C_LicksNR = []
    PrefT_M_LicksNR = []
    PrefT_C_LicksPR = []
    PrefT_M_LicksPR = []
    for index, session in enumerate(Data['Session']): #index the position in the Date list
        if Data['Session'][index] == test: 
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
- Goes throught a list of list corresponding of lick data for\
1 specific preference test and condition (total, forced or free)
- Runs lickCalc on each list to analyse licking microstructure
- returns a list of dictionaries of 23 parameters
- each list is 1 rat/1 bottle in the following order:
NR_casein (4) / NR_Malto (4) / PR_casein (4) / PR_Malto (see ListsMaker function)
   
"""
def Micro_anal(Pref_test):
    Lick_patterns_dicts = []
    for l in Pref_test:
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
and each entry contains 3 lists (1 per preference test) containing 4 dicts (one per animal)\
with all licking parameters
"""
def Sort_micro(Pref_test):
    Micro_lst = []
    for licks in Pref_test:
        Micro_lst.append(Micro_anal(licks))
        
    Micro_NR_c = []
    Micro_NR_m = []
    Micro_PR_c =[]    
    Micro_PR_m = []
    
    for pref_test in Micro_lst: #for each day of quinine test
        for condition in pref_test[0:9]:
            Micro_NR_c.append(condition)        
        for condition in pref_test[9:18]:
            Micro_NR_m.append(condition)
        for condition in pref_test[18:25]:
            Micro_PR_c.append(condition)
        for condition in pref_test[25:]:
            Micro_PR_m.append(condition)

    Micro_NR_casein = []
    for rat in range(0, len(Micro_NR_c), 9):
        split = Micro_NR_c[rat:rat+9]
        Micro_NR_casein.append(split)
    Micro_NR_malto = []
    for rat in range(0, len(Micro_NR_m), 9):
        split = Micro_NR_m[rat:rat+9]
        Micro_NR_malto.append(split)
    Micro_PR_casein = []
    for rat in range(0, len(Micro_PR_c), 7):
        split = Micro_PR_c[rat:rat+7]
        Micro_PR_casein.append(split)    
    Micro_PR_malto = []
    for rat in range(0, len(Micro_PR_m), 7):
        split = Micro_PR_m[rat:rat+7]
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
def PrefCalc(Pref_free):
    NR_C_free = []
    NR_M_free = []
    Pref_NR = []
    PR_C_free = []
    PR_M_free = []
    Pref_PR = []

    for rat in Pref_free[0]:
        a = len(rat)
        NR_C_free.append(a)
    for rat in Pref_free[1]:
        b = len(rat)
        NR_M_free.append(b)
    for index, len_licks in enumerate(NR_C_free):  
        cas = len_licks 
        malto = NR_M_free[index]
        all_free_licks = cas + malto
        Pref_NR.append(cas /all_free_licks * 100)

    for rat in Pref_free[2]:
        a = len(rat)
        PR_C_free.append(a)
    for rat in Pref_free[3]:
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
metafile="C:\\Users\\fabie\\Documents\\PPP\\PPP_metafile_preference.txt" #metafile with only preference tests
medfolder="C:\\Users\\fabie\\Documents\\PPP\\"


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
Preference tests.
Pref 1
Pref 2
Pref 3
"""
##Casein and Malto total licks (forced & free choice trials) for each preference tests
tests = ['pref1', 'pref2', 'pref3']

Pref_Licks_all = []
Pref_Licks_forced = [] #Casein and Malto forced licks for each preference tests
Pref_Licks_free = [] #Casein and Malto free licks for each preference tests

for test in tests:
    Pref_Licks_all.append(ListsMaker(test, LLicks, RLicks))
    Pref_Licks_forced.append(ListsMaker(test, LLicks_forced, RLicks_forced))
    Pref_Licks_free.append(ListsMaker(test, LLicks_free, RLicks_free))


##Calculate casein preference for each rat for a each preference test\
#based on licks during free choice trials
#Split the results in separate lists for NR and PR
Cas_Pref = []
for Test in Pref_Licks_free:
    Cas_Pref.append(PrefCalc(Test)) #calculation of casein pref

Cas_Pref_NR = []
Cas_Pref_PR = []
for pref in Cas_Pref:
    Cas_Pref_NR.append(pref[0]) #isolate casein preference for NR group
    Cas_Pref_PR.append(pref[1]) #isolate casein preference for PR group
    
##Add Preference test number as a key for both NR and PR groups casein preferences    
Cas_Pref_NR_tests = {'PrefTest1': Cas_Pref_NR[0], 'PrefTest2': Cas_Pref_NR[1], \
                     'PrefTest3': Cas_Pref_NR[2]}
Cas_Pref_PR_tests = {'PrefTest1': Cas_Pref_PR[0], 'PrefTest2': Cas_Pref_PR[1], \
                     'PrefTest3': Cas_Pref_PR[2]}
     
##Licking microstructure during forced and free choice trials
Micro_licks_forced = Sort_micro(Pref_Licks_forced)
Micro_licks_free = Sort_micro(Pref_Licks_free)

#Sort and export final results (1 list for each parameters with 192 values:\
#0-47: NR-casein; 48-95: NR-malto; 96-143: PR-casein; 144-191: PR-malto)
bMean_all_forced, bNum_all_forced, clustMean_all_forced, \
clustNum_all_forced, freq_all_forced, total_all_forced = allratlistmaker(Micro_licks_forced) #Forced trials

bMean_all_free, bNum_all_free, clustMean_all_free, \
clustNum_all_free, freq_all_free, total_all_free = allratlistmaker(Micro_licks_free) #Free trials
