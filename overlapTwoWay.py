"""
Date - 7/20/17

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
__version__ = '0.0.1'
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


def overLap2Way(uniA,uniB,sigA,sigB,chrList):
    """
    Overlaps the interval objects from gappedPeakReader()
    and outputs the variables required for a fisher exact test.
    Also outputs the intervals associated with that group.

    :param uniA: All intervals in A
    :param uniB: All intervals in B
    :param sigA: Intervals that meet the significance cut off for A
    :param sigB: Intervals that meet the significance cur off for B
    :param chrList: List of Chromosomes to analyze
    :return: The returned dictionaries are organized by chromosomes.
    """
    sigOverlap = dict.fromkeys(chrList)
    sigOnlyA = dict.fromkeys(chrList)
    sigOnlyB = dict.fromkeys(chrList)
    universe = dict.fromkeys(chrList)
    for i in chrList:
        sigOverlap[i] = {
            'Total Coverage': 0
        }
        sigOnlyA[i] = {
            'Total Coverage': 0
        }
        sigOnlyB[i] = {
            'Total Coverage': 0
        }
        universe[i] = {
            'Total Coverage': 0,
        }

    for i in chrList:
        ### Significant Overlap
        sigIntA = sigA[i]["Peaks"]
        sigIntB = sigB[i]["Peaks"]
        overlap = sigIntA & sigIntB
        sigCover = 0
        for j in overlap:
            sigCover = sigCover + (j[1] - j[0])
        sigOverlap[i]["Total Coverage"] = sigCover
        ## Universe
        uniOverlap = uniA[i]["Peaks"] & uniB[i]["Peaks"]
        uniCover = 0
        for k in uniOverlap:
            uniCover = uniCover + (k[1] - k[0])
        ### Sig A Only
        overlapA = uniOverlap & sigIntA
        coverA = 0
        for l in overlapA:
            coverA = coverA + (l[1] - l[0])
        sigOnlyA[i]["Total Coverage"] = coverA - sigCover
        ### Sig B Only
        overlapB = sigIntB & uniOverlap
        coverB = 0
        for m in overlapB:
            coverB = coverB + (m[1] - m[0])
        sigOnlyB[i]["Total Coverage"] = coverB - sigCover
        ## Universe
        uniCover = (((uniCover - sigCover) - sigOnlyA[i]["Total Coverage"]) - sigOnlyB[i]["Total Coverage"])
        universe[i]["Total Coverage"] = uniCover

    return sigOverlap, sigOnlyA, sigOnlyB, universe


def fisherTest(sigOverlap,sigOnlyA,sigOnlyB,universe,chrList):
    """
    Performs a fisher exact test on the overlap results
    returning pValue and Odds Ratio
    :param sigOverlap: The Significant Overlaps for A and B
    :param sigOnlyA: The number of sig basepairs for only A
    :param sigOnlyB: The number of sig basepairs for only B
    :param universe: The Universe of non significant basepairs with peaks.
    :return: pValue and Odds Ratio
    """
    A = 0
    B = 0
    C = 0
    D = 0
    for i in chrList:
        A = A + sigOverlap[i]["Total Coverage"]
        B = B + sigOnlyA[i]["Total Coverage"]
        C = C + sigOnlyB[i]["Total Coverage"]
        D = D + universe[i]["Total Coverage"]
    A = int(A)
    B = int(B)
    C = int(C)
    D = int(D)
    oddsRatio, pValue = spst.fisher_exact([[A, B], [C, D]])

    return pValue,oddsRatio, A, B, C, D


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
    sigOverlap, sigOnlyA, sigOnlyB, universe = overLap2Way(uniA, uniB, sigA, sigB, chrList)
    pValue, oddsRatio, A, B, C, D = fisherTest(sigOverlap, sigOnlyA, sigOnlyB, universe, chrList)
    writeOutput(A, B, C, D, fileA, fileB, outFile, pValue, oddsRatio,chrList)
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

