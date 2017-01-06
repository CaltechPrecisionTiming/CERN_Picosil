import numpy as np
import os
import sys
from glob import glob

def macro(datdir, nEvents, runNums):
    '''Analyzes .dat files, and outputs them as .root.
    Run this from directory containing dat2rootCP.
    dat dir: directory containing .dat files, as a string.
    nEvents: maximum number of events to be processed in any individual file.
             i.e. nEvents = max(2000,5000,2000) = 5000
    
    e.g. macro('raw/', 5000, [0]) OR macro('raw/', 5000, [1,2,3,5,7])
    Then the output .root files would be in the raw/ subdirectory.
    '''
    
    analyzer = './dat2rootCP '

    
    datdir_array = [] # list of dat files with their dir's.
    
    if runNums != [0]:
        datfile_array = runNums
        for i in range(len(datfile_array)):
            datfile_array[i] = 't1065-jun-2016-' + str(datfile_array[i]) + '.dat'
        for i in datfile_array:
            datdir_array.append(datdir + '/' + i)
    else:
        datdir_array = glob(datdir + '/t1065-jun-2016-*.dat')
            
    for i in datdir_array:
        final = analyzer + i + ' ' + str(nEvents)
        print 'Entering command: ' + final
        os.system(final)
    

if __name__ == '__main__':
    
    if len(sys.argv) < 4:
        print "usage:"
        print " to analyze runs 1-3 of <= 5000 events each in raw subdir:"
        print " python python/batch_analyze.py raw/ 5000 1-3,5,7"
        print " enter '0' instead of last argument for all .dat files"
        sys.exit()
        
    datdir = sys.argv[1] # 'raw/'
    nEvents = int(sys.argv[2]) # 5000
    
    strNumbersToConvert = sys.argv[3] # '1-3,5,7' ## Enter 0 for ALL .dat files
    strNumbersToConvertList = strNumbersToConvert.split(",") #['1-3','5','7']
    
    intNumersToConvertList = [] # Will be [1,2,3,5,7]
    
    for strNum in strNumbersToConvertList:
        if '-' in strNum:
            num1,num2 = strNum.split("-")
            num1 = int(num1)
            num2 = int(num2)
            fill = np.arange(num2 + 1)[num1:] # [1,2,3]
            for i in fill:
                intNumersToConvertList.append(i)
        else:
            intNumersToConvertList.append(int(strNum))

    macro(datdir, nEvents, intNumersToConvertList)
