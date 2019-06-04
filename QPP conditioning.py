# -*- coding: utf-8 -*-
"""
Data extraction and analyses for Quinine Protein Preference QPP
Conditioning data

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
Uses 2 list of licks (Left and Right)
Returns 4 lists of n animals as following:
NR_casein / NR_malto / PR_casein / NR_malto
Each list contains n values*2 (1 value for each conditioning day)
here 4 rats -> 8 values [0:3] day 1, [4:7] day 2
"""    
def ListsMaker(Llicks, RLicks):
    Cond_C_LicksNR = []
    Cond_M_LicksNR = []
    Cond_C_LicksPR = []
    Cond_M_LicksPR = []
    
    for index, session in enumerate(Data['Date']): #index the position in the Date list
        if 'Cas' in Data['Left'][index]:
               if 'NR' in Data['Diet'][index]:
                      Cond_C_LicksNR.append(LLicks[index])
               if 'PR' in Data['Diet'][index]:
                      Cond_C_LicksPR.append(LLicks[index])
        if 'Cas' in Data['Right'][index]:
               if 'NR' in Data['Diet'][index]:
                      Cond_C_LicksNR.append(RLicks[index])
               if 'PR' in Data['Diet'][index]:
                      Cond_C_LicksPR.append(RLicks[index]) 
                      
        if 'Malt' in Data['Left'][index]:
               if 'NR' in Data['Diet'][index]:
                      Cond_M_LicksNR.append(LLicks[index])
               if 'PR' in Data['Diet'][index]:
                      Cond_M_LicksPR.append(LLicks[index])
        if 'Malt' in Data['Right'][index]:
               if 'NR' in Data['Diet'][index]:
                      Cond_M_LicksNR.append(RLicks[index])
               if 'PR' in Data['Diet'][index]:
                      Cond_M_LicksPR.append(RLicks[index])                      
           
    return([Cond_C_LicksNR, Cond_M_LicksNR, Cond_C_LicksPR, Cond_M_LicksPR])
    
"""
Micro_anal function:
- Goes throught a list of list (Quin_cond) corresponding of lick data for\
1 specific quinine day and condition (total, forced or free)
- Runs lickCalc on each list to analyse licking microstructure
- returns a list of dictionaries of 23 parameters
- each list is 1 rat/1 bottle
   
"""
def Micro_anal(lickslist):
    Lick_patterns_dicts = []
    for l in lickslist:
        for lst in l:
            x = GF.lickCalc(lst)
            Lick_patterns_dicts.append(x)
    return Lick_patterns_dicts

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


    for rat_dict in licking_patterns:
         bMean_all.append(rat_dict['bMean'])
         bNum_all.append(rat_dict['bNum'])
         clustMean_all.append(rat_dict['clustMean'])
         clustNum_all.append(rat_dict['clustNum'])
         freq_all.append(rat_dict['freq'])
         total_all.append(rat_dict['total'])
                
    return bMean_all, bNum_all, clustMean_all, clustNum_all, freq_all, total_all
    
"-------------------------------------------------------------------------------"
## Metadata folder (check pathway for each exp)
#metafile="C:\\Users\\fabie\\Documents\\QPP data\\QPPh1_metafile.txt" #metafile with all data
metafile="D:\\QPPh1\\QPPh1_metafile_conditioning.txt" #metafile with only preference tests
medfolder="D:\\QPPh1\\Behaviour\\"

## Extraction licks ('b' = licks left; 'e'= licks right) from medfile
Data = Extract(metafile)
Licksallrats = []
for Medfile in Data['Medfile']:
    Licksallrats.append(GF.medfilereader(medfolder+Medfile, varsToExtract = ['b', 'e'], remove_var_header = True))

TotalLicks = []
for l in Licksallrats:
    TotalLicks.append(len(l[0]) + len(l[1]))
    

LLicks = [] #left licks
RLicks = [] #right licks

for session in Licksallrats:
    LLicks.append(session[0])
    RLicks.append(session[1])

#Extract licks for each bottle and conditioning day    
Conditioning = ListsMaker(LLicks, RLicks)

##Licking microstructure during forced and free choice trials
Micro_licks_cond = Micro_anal(Conditioning)

#Sort and export final results (1 list for each parameters with 32 values:\
#0-7: NR-casein d1 and d2; 8-15: NR-malto d1 and d2; 16-23: PR-casein d1 and d2; \
# 24-31: PR-malto d1 and 2)
bMean_all_cond, bNum_all_cond, clustMean_all_cond, clustNum_all_cond, \
freq_all_cond, total_all_cond = allratlistmaker(Micro_licks_cond)

