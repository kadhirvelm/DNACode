from __future__ import division
import sys
import re
import random

from Bio import SeqIO


def main():
    size = 21
    mySeq = loadFasta()  # load the fasta file
    myDict, indices = parseFasta(mySeq, size)  # parse the fasta file for .{21}GG (regex)
    myDict = callOnParsers(myDict, mySeq)
    writeFile(myDict, indices)


def loadFasta():
    # Opening the input file, with exception handling.
    try:
        fullSeq = ""
        inFile = open(sys.argv[1])
        for seq in SeqIO.parse(inFile, "fasta"):
            fullSeq += seq.seq
    except IndexError:
        print("You must specify a fasta filename to run.")
        sys.exit()
    return str(fullSeq)


def parseFasta(mySeq, size):
    indices = {}
    matchDict = {}
    newMatchDict = {}
    currIndex = []
    pattern = ".{" + str(size) + "}GG"
    m = re.findall(pattern, mySeq)
    for kmer in re.finditer(pattern, mySeq):
        matchDict[kmer.group(0)] = 0
        indices[kmer.group(0)] = str(kmer.start()) + " " + str(kmer.end())
    for x in range(0, 100):
        index = random.randint(x, len(m) - 1)
        if index in currIndex:
            index = x
            currIndex.append(x)
        else:
            currIndex.append(index)
        newMatchDict[m[index]] = 0
    return newMatchDict, indices


def parseOffTarget(mySeq, myDict,kmer):
    total = (len(re.findall(kmer, mySeq))-1)*10  # arbitrary
    myDict[kmer] += total
    return myDict


def parseOffTarget1(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    total = 0
    for bp in range(start, end):
        temp1 = temp + helperResEnd(start, bp, kmer)
        total += len(re.findall(temp1, mySeq))/numOffTargets
    myDict[kmer] += total
    return myDict


def parseOffTarget2(mySeq, myDict, start, end, temp, kmer,numOffTargets):
    for bp in range(start, end - 1):
        temp2 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget1(mySeq, myDict, bp + 1, end, temp2, kmer,numOffTargets)
    return myDict


def parseOffTarget3(mySeq, myDict, start, end, temp, kmer,numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget2(mySeq, myDict, bp + 1, end, temp1, kmer,numOffTargets)
    return myDict


def parseOffTarget4(mySeq, myDict, start, end, temp, kmer,numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget3(mySeq, myDict, bp + 1, end, temp1, kmer,numOffTargets)
    return myDict


def parseOffTarget5(mySeq, myDict, start, end, temp, kmer,numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget4(mySeq, myDict, bp + 1, end, temp1, kmer,numOffTargets)
    return myDict


def parseOffTarget6(mySeq, myDict, start, end, temp, kmer,numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget5(mySeq, myDict, bp + 1, end, temp1, kmer,numOffTargets)
    return myDict


def parseOffTarget7(mySeq, myDict, start, end, kmer,numOffTargets):
    for bp in range(start, end - 2):
        temp = helperRes(start, bp, kmer)
        myDict = parseOffTarget6(mySeq, myDict, bp + 1, end, temp, kmer,numOffTargets)
    return myDict


def callOnParsers(myDict, mySeq):
    counter = 0
    print("Status:")
    for kmer in myDict:
        totalOff = 5
        totalOffStandard = 5
        while totalOff > -1:
            if totalOff == 7:
                parseOffTarget7(mySeq, myDict, 0, len(kmer) - 2, kmer, totalOff)
            elif totalOff == 6:
                parseOffTarget6(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 5:
                parseOffTarget5(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 4:
                parseOffTarget4(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 3:
                parseOffTarget3(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 2:
                myDict = parseOffTarget2(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 1:
                myDict = parseOffTarget1(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 0:
                myDict = parseOffTarget(mySeq, myDict, kmer)
            totalOff -= 1
            percent = (1/(len(myDict)+1))*(totalOffStandard-totalOff)/totalOffStandard
            percent = ((counter/len(myDict)) + percent)*100
            percent = round(percent, 2)
            print(str(percent)+"%")
        counter += 1
        percent = (counter/len(myDict))*100
        percent = round(percent, 2)
        print(str(percent)+"%")
    return myDict


def helperRes(start, bp, kmer):
    return kmer[start:bp] + regExDefiner(kmer[bp])


def helperResEnd(start, bp, kmer):
    return helperRes(start, bp, kmer) + kmer[bp + 1:]


def regExDefiner(char):
    if char == 'A':
        return '[T,C,G]'
    if char == 'T':
        return '[A,C,G]'
    if char == 'C':
        return '[A,T,G]'
    if char == 'G':
        return '[A,T,C]'


def writeFile(myDict, indicies):
    counter = 0
    sortDict = sorted(myDict, key=myDict.get)
    fileTemp = open('output.txt', 'w+')
    fileTemp.write("Crispi -- Version 66.0123 Beta \n")
    fileTemp.write("Note: \n-->Indicies are exclusive. "
                       + "\n-->The lower the score, the better the site, aka the less off-target probability\n\n")
    fileTemp.write("Sequence (start, end) Score: --\n")
    for kmer in sortDict:
        counter += 1
        fileTemp.write(str(kmer) + " ( " + indicies[kmer] + " )" + " Score: " + str(myDict[kmer])+"\n")
        if counter == 10:
            break


main()
