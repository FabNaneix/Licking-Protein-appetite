# -*- coding: utf-8 -*-
"""
Created on Thu May 30 15:45:19 2019

@author: fabie
"""
bMean_all = []
@ -374,6 +372,53 @@ for condition_name in Micro_licks_forced:
            clustNum_all.append(rat_dict['clustNum'])
            freq_all.append(rat_dict['freq'])
            total_all.append(rat_dict['total'])
  count = 0 
    bMean_all = []
    bNum_all = []
    clustMean_all = []
    clustNum_all = []
    freq_all = []
    total_all = [] 
    
    for index, condition in enumerate(Micro_licks_dict):
        for quinine_conc in condition:
            #print(quinine_conc)
            for rat_dict in quinine_conc:
                #print(type(rat_dict))
               # print(rat_dict)

                count += 1 
            
#                bMean_all.append(rat['bMean'])
#                bNum_all.append(rat['bNum'])
#                clustMean_all.append(rat['clustMean'])
#                clustNum_all.append(rat['clustNum'])
#                freq_all.append(rat['freq'])
#                total_all.append(rat['total'])
                
 #   return(bMean_all, bNum_all, clustMean_all, clustNum_all, freq_all, total_all)

    return count

bMean_all_forced, bNum_all_forced, clustMean_all_forced, \
clustNum_all_forced, freq_all_forced, total_all_forced = allratlistmaker(Micro_licks_forced)


bMean_all = []
bNum_all = []
clustMean_all = []
clustNum_all = []
freq_all = []
total_all = []

for condition_name in Micro_licks_forced:
    for quinine_conc in Micro_licks_forced[condition_name]:
        for rat_dict in quinine_conc:
            bMean_all.append(rat_dict['bMean'])
            bNum_all.append(rat_dict['bNum'])
            clustMean_all.append(rat_dict['clustMean'])
            clustNum_all.append(rat_dict['clustNum'])
            freq_all.append(rat_dict['freq'])
            total_all.append(rat_dict['total'])