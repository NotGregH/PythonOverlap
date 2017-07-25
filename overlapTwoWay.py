"""
Date - 7/25/17

Writen for the Gamble Lab @
Albert Einstein College of Medicine

Simple program for overlapping gappedPeak objects
created by gappedPeakReader(). The program will output
to a fisher exact test result text file.

NOTE : pyinterval package will need to be downloaded.

Output Fisher Exact Test :

FileName A vs FileName B
Chromosomes overlapped: Autosomes, All or All + m
Significant Overlap Base Pairs : #
Significant A not in Sig Overlap : #
Significant B not in Sig Overlap : #
Universe not in any Sig : #
pValue : float
Odds Ratio : float


Future Revisions:
 - Convert gappedPeakFile into a class rather than a function
    - The current version trashes the columns not required for overlap
- Profile the code to speed it up
- Convert into a D program

"""

__author__ =  'Gregory A. Hamilton'
__version__ = '0.0.2'
__license__ = ''
__email__ = 'ghamilto@mail.einstein.yu.edu'

###################
# Required Modules
###################

import interval
import argparse
from csv import writer
from gappedPeakFile import gappedPeakReader
import scipy.stats as spst

######################
# Functions
######################


def pickChromosomes(selection='A'):
    """
    Selects the chromosomes you want to overlap.
    :param selection: A(autosomes), AS (Auto+Sex), or ASM (Auto+Sex+M)
    :return: chrList
    """
    chrList = [
        '1', '2', '3', '4', '5',
        '6', '7', '8', '9', '10',
        '11', '12', '13', '14', '15',
        '16', '17', '18', '19', '20',
        '21', '22'
               ]
    if selection == "AS":
        chrList.append("x")
        chrList.append("y")
    if selection == "ASM":
        chrList.append("x")
        chrList.append("y")
        chrList.append("m")
    return chrList


def intervalCoverage(interval):
    """
    Calculates the base pairs covered by an interval.
    :param interval:
    :return:
    """
    coverage = 0
    for i in interval:
        coverage = coverage + (i[1] - i[0])
    return coverage


def overLap2Way(uniA,uniB,sigA,sigB):
    """
    Overlaps the interval objects and
    outputs the variables required for a fisher exact test.
    Also outputs the intervals associated with that group.

    :param uniA: All intervals in A
    :param uniB: All intervals in B
    :param sigA: Intervals that meet the significance cut off for A
    :param sigB: Intervals that meet the significance cur off for B
    :return: A, B, C, and D for a fishers exact test
    """
    yesYes = 0
    yesNo = 0
    noYes = 0
    noNo = 0

    universe = uniA & uniB
    noNo = intervalCoverage(universe)
    yesYesInterval = sigA & sigB
    yesYes = intervalCoverage(yesYesInterval)
    yesNoInterval = sigA & universe
    yesNo = intervalCoverage(yesNoInterval) - yesYes
    noYesInterval = sigB & universe
    noYes = intervalCoverage(noYesInterval) - yesYes

    noNo = (((noNo - yesNo) - yesYes) - noYes)
    return yesYes, yesNo, noYes, noNo


def compileOverlapData(uniA,uniB,sigA,sigB,chrList):
    """
    Piles-up the overlap coverage data for the fisher exact test.

    :param uniA: Interval Universe A
    :param uniB: Interval Universe B
    :param sigA: Interval Significant A
    :param sigB: Interval Significant B
    :param chrList: List of chromosomes to analyze
    :return:
    """
    yesYes = 0
    yesNo = 0
    noYes = 0
    noNo = 0
    for i in chrList:
        a, b, c, d = overLap2Way(uniA[i]["Peaks"], uniB[i]["Peaks"], sigA[i]["Peaks"], sigB[i]["Peaks"])
        yesYes = yesYes + a
        yesNo = yesNo + b
        noYes = noYes + c
        noNo = noNo + d

    return yesYes, yesNo, noYes, noNo


def fisherTest(yesYes,yesNo,noYes,noNo):
    """
    Performs a fisher exact test on the overlap results
    returning pValue and Odds Ratio
    :param yesYes: The Significant Overlap Basepairs for A and B
    :param yesNo: The number of sig basepairs for only A
    :param noYes: The number of sig basepairs for only B
    :param noNo: The Total Universe of non significant basepairs.
    :return: pValue and Odds Ratio
    """
    oddsRatio, pValue = spst.fisher_exact([[yesYes, yesNo], [noYes, noNo]])

    return pValue,oddsRatio


def writeOutput(A, B, C, D, fileA, fileB,
                outFile, pValue, oddsRatio, chrList):
    """
    Writes the fisher's exact test results to the
    designated output file.
    """
    fileA = fileA.replace("segGappedPeak", "")
    fileB = fileB.replace("segGappedPeak", "")
    with open(outFile,'w') as f:
        write = writer(f,delimiter="\t",lineterminator='\n')
        write.writerow((fileA,"vs",fileB))
        write.writerow(("Chromosomes Analyzed", ":" ))
        write.writerow(chrList)
        write.writerow(("Significant Overlap :",A))
        write.writerow(("Significant", fileA,"only :", B))
        write.writerow(("Significant", fileB,"only :", C))
        write.writerow(("Universe :", D))
        write.writerow(("p-Value :", pValue))
        write.writerow(("Odds Ratio :", oddsRatio))

    return

###########################################
# Main function
###########################################

def main(fileA, fileB, pValFilterA, signalFilterA,
         pValFilterB, signalFilterB, outFile, selection):
    chrList = pickChromosomes(selection)
    uniA, sigA = gappedPeakReader(chrList, fileA, pValFilterA, signalFilterA)
    uniB, sigB = gappedPeakReader(chrList, fileB, pValFilterB, signalFilterB)
    yesYes, yesNo, noYes, noNo = compileOverlapData(uniA, uniB, sigA, sigB, chrList)
    pValue, oddsRatio = fisherTest(yesYes, yesNo, noYes, noNo)
    writeOutput(yesYes, yesNo, noYes, noNo, fileA, fileB, outFile, pValue, oddsRatio,chrList)
    return


##########################################
##########################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-fA','--fileA')
    parser.add_argument('-fB', '--fileB')
    parser.add_argument('-pA','--pFilterA', default=0.05, type = float)
    parser.add_argument('-pB','--pFilterB' , default=0.05, type = float)
    parser.add_argument('-sA','--sFilterA', default = 1, type= float)
    parser.add_argument('-sB','--sFilterB', default = 1, type = float)
    parser.add_argument('-o','--outFile')
    parser.add_argument('-chr','--chromo',default='A')
    arg = parser.parse_args()
    main(arg.fileA, arg.fileB, arg.pFilterA, arg.sFilterA,
         arg.pFilterB, arg.sFilterB, arg.outFile, arg.chromo)

