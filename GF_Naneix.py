# -*- coding: utf-8 -*-
"""
General functions
Created on Wed May  8 16:21:11 2019

@author: Fabien Naneix (based on Jaime McCutcheon)
"""

import numpy as np

### medfilereader is used to extract data from the medfile datasets
"""
medfilereader takes the following arguments:
    filename - including path
    sessionToExtract - 1 is default, for situations in which more than one session is included in a single file
    varsToExtract - as strings, number of outputs must match
    verbose - False is default
    remove_var_header - False is default, removes first value in array, useful when negative numbers are used to signal array start
"""        
def medfilereader(filename, varsToExtract = 'all',
                  sessionToExtract = 1,
                  verbose = False,
                  remove_var_header = False):
    if varsToExtract == 'all':
        numVarsToExtract = np.arange(0,26)
    else:
        numVarsToExtract = [ord(x)-97 for x in varsToExtract]
    
    f = open(filename, 'r')
    f.seek(0)
    filerows = f.readlines()[8:]
    datarows = [isnumeric(x) for x in filerows]
    matches = [i for i,x in enumerate(datarows) if x == 0.3]
    if sessionToExtract > len(matches):
        print('Session ' + str(sessionToExtract) + ' does not exist.')
    if verbose == True:
        print('There are ' + str(len(matches)) + ' sessions in ' + filename)
        print('Analyzing session ' + str(sessionToExtract))
    
    varstart = matches[sessionToExtract - 1]
    medvars = [[] for n in range(26)]
    
    k = int(varstart + 27)
    for i in range(26):
        medvarsN = int(datarows[varstart + i + 1])
        
        medvars[i] = datarows[k:k + int(medvarsN)]
        k = k + medvarsN
        
    if remove_var_header == True:
        varsToReturn = [medvars[i][1:] for i in numVarsToExtract]
    else:
        varsToReturn = [medvars[i] for i in numVarsToExtract]

    if np.shape(varsToReturn)[0] == 1:
        varsToReturn = varsToReturn[0]
    return varsToReturn

def isnumeric(s):
    try:
        x = float(s)
        return x
    except ValueError:
        return float('nan')

### metafilereader is used to read the data from the metafile before to combine with medile data
"""
metafilereader takes the following argument:
    filename - including path
"""        
def metafilereader(filename):
    
    f = open(filename, 'r')
    f.seek(0)
    header = f.readlines()[0]
    f.seek(0)
    filerows = f.readlines()[1:]
    
    tablerows = []
    
    for i in filerows:
        tablerows.append(i.split('\t'))
        
    header = header.split('\t')
    # need to find a way to strip end of line \n from last column - work-around is to add extra dummy column at end of metafile
    return tablerows, header
    
### lickCalc is used to extract lick parameters
"""
This function will calculate data for bursts from a train of licks. The threshold
for bursts and clusters can be set. It returns all data as a dictionary.

lickCalc uses the following arguments:
    licks as a list []
    offset is used to detect and isoltate the longlicks
    burstThreshold defines the ILI to separate bursts (usually 0.25sec)
    runThreshold defines the ILI to separate clusters (usually 0.5sec)
    adjustforlonglicks adjusts lick data base on longlicks detection with the offset
"""
def lickCalc(licks, offset = [], burstThreshold = 0.25, clustThreshold = 0.5, 
             adjustforlonglicks='none'): ### makes dictionary of data relating to licks, bursts and clusters. threshold is determined in the arguments
    if type(licks) != np.ndarray or type(offset) != np.ndarray:
        try:
            licks = np.array(licks)
            offset = np.array(offset)
        except:
            print('Licks and offsets need to be arrays and unable to easily convert.')
            return

    lickData = {}
    
    if len(offset) > 0:      ##detect and isolate long licks  
        lickData['licklength'] = offset - licks[:len(offset)]
        lickData['longlicks'] = [x for x in lickData['licklength'] if x > 0.3]
    else:
        lickData['licklength'] = []
        lickData['longlicks'] = []
    
    if adjustforlonglicks != 'none':
        if len(lickData['longlicks']) == 0:
            print('No long licks to adjust for.')
        else:
            lickData['median_ll'] = np.median(lickData['licklength'])
            lickData['licks_adj'] = int(np.sum(lickData['licklength'])/lickData['median_ll'])
            if adjustforlonglicks == 'interpolate':
                licks_new = []
                for l, off in zip(licks, offset):
                    x = l
                    while x < off - lickData['median_ll']:
                        licks_new.append(x)
                        x = x + lickData['median_ll']
                licks = licks_new
    try:    
        lickData['licks'] = licks #find timing licks
        lickData['ilis'] = np.diff(np.concatenate([[0], licks])) #find inter-licks intervals ILIs
        lickData['shilis'] = [x for x in lickData['ilis'] if x < burstThreshold] #find short(not supposed to be) ILIs
        lickData['freq'] = 1/np.mean([x for x in lickData['ilis'] if x < burstThreshold]) #calculate licks frequency
        lickData['total'] = len(licks) #total number of licks
        
    ### Calculates start, end, number of licks and time for each BURST 
        lickData['bStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)] #start of bursts  
        lickData['bInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > burstThreshold)] #find indexation of burst
        lickData['bEnd'] = [lickData['licks'][i-1] for i in lickData['bInd'][1:]] #find end of bursts
        lickData['bEnd'].append(lickData['licks'][-1])
    
        lickData['bLicks'] = np.diff(lickData['bInd'] + [len(lickData['licks'])]) #nbr licks per bursts   
        lickData['bTime'] = np.subtract(lickData['bEnd'], lickData['bStart']) #bursts duration
        lickData['bNum'] = len(lickData['bStart']) #bursts number
        if lickData['bNum'] > 0:
            lickData['bMean'] = np.nanmean(lickData['bLicks']) #mean licks/bursts
        else:
            lickData['bMean'] = 0
        
        lickData['bILIs'] = [x for x in lickData['ilis'] if x > burstThreshold] #ILIs within bursts
    
        ### Calculates start, end, number of licks and time for each CLUSTER (parameters are the same than for bursts)
        lickData['clustStart'] = [val for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > clustThreshold)]  
        lickData['clustInd'] = [i for i, val in enumerate(lickData['licks']) if (val - lickData['licks'][i-1] > clustThreshold)]
        lickData['clustEnd'] = [lickData['licks'][i-1] for i in lickData['clustInd'][1:]]
        lickData['clustEnd'].append(lickData['licks'][-1])
    
        lickData['clustLicks'] = np.diff(lickData['clustInd'] + [len(lickData['licks'])])    
        lickData['clustTime'] = np.subtract(lickData['clustEnd'], lickData['clustStart'])
        lickData['clustNum'] = len(lickData['clustStart'])
        if lickData['clustNum'] > 0:
            lickData['clustMean'] = np.nanmean(lickData['clustLicks'])
        else:
            lickData['clustMean'] = 0
            
        lickData['clustILIs'] = [x for x in lickData['ilis'] if x > clustThreshold]        
    
    except:
        lickData['freq'] = 0
        lickData['total'] = 0
        lickData['bNum'] = 0
        lickData['bMean'] = 0
        lickData['clustNum'] = 0
        lickData['clustMean'] = 0
        
        
    return lickData


