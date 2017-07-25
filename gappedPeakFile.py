"""
Date - 7/20/17

Writen for the Gamble Lab @
Albert Einstein College of Medicine

Simple program for reading gappedPeak Files
into 'universe' and 'significant' gappedPeak objects.

NOTE : pyinterval package will need to be downloaded.


Output obj Structure:

{ Chr :
    {
    'TotalCoverage' : int(bp's for the chr) ,
    'Peaks' : intervalObj,
    }
}



Future Revisions:
- Converted gappedPeakFile into a class rather than a function
    - The class will include all of the parts of the gappedPeak
    - The current version trashes the columns not required for overlap
- Profile the code to speed it up
- Convert into a D program

"""

__author__ =  'Gregory A. Hamilton'
__version__ = '0.0.1'
__license__ = ''
__email__ = 'ghamilto@mail.einstein.yu.edu'



#####################
# Required Modules
#####################

import interval
import argparse
from csv import writer


#####################
# Functions
#####################

def gappedPeakReader(chrList, gappedPeakFile, pValFilter=.05, signalFilter=0):
    """
    Reads a gapped Peak file and outputs the 'interval'
    universe and significant 'intervals' as objects described above.

    :param gappedPeakFile:
    :param pValFilter: The pValue cut off you'd like to
    :param signalFilter: The ISOR signal value you'd like to filter by for
        significant 'intervals'
    :param chrList: List of chromosomes you would like to analyze ( must
        be just the autosome number, x, y or m ) i.e. 1 not chr1
    :return:
    """
    sigInt = interval.interval()
    uniInt = interval.interval()
    universeIntervals = dict.fromkeys(chrList)
    sigIntervals = dict.fromkeys(chrList)

    for i in chrList:
        universeIntervals[i] = {
            'Total Coverage': 0,
            'Peaks': uniInt
        }
        sigIntervals[i] = {
            'Total Coverage': 0,
            'Peaks': sigInt
        }

    with open(gappedPeakFile) as infile:
        for line in infile:
            line = line.strip().lower().split('\t')
            line[0] = line[0].replace('chr', '')

            if line[0] not in chrList:
                continue

            else:
                chrom = line[0]
                chrStart = int(line[1])
                chrStop = int(line[2])

                if chrStart > chrStop:
                    continue

                blockSizes = line[10]
                blockStart = line[11]
                signalVal = float(line[12])
                chrInt = interval.interval()

                if 'inf' in line[13]:
                    pVal = float('inf')

                else:
                    pVal = float(line[13])

                #if signalVal > 0:  # Delete this line if we want both negative and positive signal in the universe
                blockStart = blockStart.replace(' ', '').split(',')
                blockSizes = blockSizes.replace(' ', '').split(',')
                count = 0

                while count < len(blockStart):
                    blockStart[count] = int(blockStart[count])
                    blockSizes[count] = int(blockSizes[count])
                    gapStop = (blockStart[count] + blockSizes[count]) + chrStart
                    gapStart = (blockStart[count] + chrStart)
                    tmpInt = interval.interval([gapStart, gapStop])
                    chrInt = tmpInt | chrInt
                    count = count + 1
                coverage = 0
                for i in chrInt:
                    coverage = coverage + (i[1] - i[0])

                universeIntervals[chrom]['Total Coverage'] = universeIntervals[chrom]['Total Coverage'] + coverage
                universeIntervals[chrom]['Peaks'] = universeIntervals[chrom]['Peaks'] | chrInt

                #if signalVal > signalFilter:
                if pVal < pValFilter and signalVal > signalFilter:
                    sigIntervals[chrom]['Total Coverage'] = sigIntervals[chrom]['Total Coverage'] + coverage
                    sigIntervals[chrom]['Peaks'] = sigIntervals[chrom]['Peaks'] | chrInt

    return universeIntervals, sigIntervals

def testOutput(universeIntervals,sigIntervals,outFile,chrList):
    sigFile = outFile + '.sig'
    uniFile = outFile + '.uni'
    with open(sigFile, "w",-1) as f:
        write = writer(f,delimiter='\t')
        for i in chrList:
            cover = sigIntervals[i]['Total Coverage']
            write.writerow(("Chr :", i))
            write.writerow(("Total Coverage :", cover))
            intervals = sigIntervals[i]['Peaks']
            write.writerow(("Peaks", ""))
            for j in intervals:
                write.writerow(("Start:", j[0], "stop:", j[1]))
        f.close()

    with open(uniFile, "w",-1) as f:
        write = writer(f,delimiter='\t')
        for i in chrList:
            cover = universeIntervals[i]['Total Coverage']
            write.writerow(("Chr :", i))
            write.writerow(("Total Coverage :", cover))
            intervals = universeIntervals[i]['Peaks']
            write.writerow(("Peaks", ""))
            for j in intervals:
                write.writerow(("Start:", j[0], "stop:", j[1]))
        f.close()
    return


#####################
# Main
#####################

if __name__ == '__main__':
    ### This is meant for testing the code
    parser = argparse.ArgumentParser()
    parser.add_argument("-f","--inFile")
    parser.add_argument("-o","--oFile")
    arg = parser.parse_args()
    chrList = ['1','2','3','4','5']
    universeIntervals, sigIntervals = gappedPeakReader(chrList, arg.inFile)
    testOutput(universeIntervals, sigIntervals, arg.oFile, chrList)

