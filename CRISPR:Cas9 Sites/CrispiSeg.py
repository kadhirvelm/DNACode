##! /usr/bin/python
from __future__ import division
import sys
import random

import regex
from Bio import SeqIO


def main():
    size = 21
    print("Loading fasta file...")
    mySeq, segmentSeq, name = loadFasta()  # load the fasta file
    print("Fasta file loaded.")
    print("Total sequence length: " + str(len(mySeq)))
    maxMismatches, numberOutputted, cas9Site = loadParameters()
    print("Locating and parsing NGG sites...")
    myDict, indices, totalCas9 = parseFasta(segmentSeq, size, cas9Site)  # parse the fasta file for .{21}GG
    print("NGG Sites found.")
    print("Calculating off-targets...")
    myDict = callOnParsers(myDict, mySeq, maxMismatches)
    print("Calculated off-targets.")
    print("Writing output file...")
    writeFile(myDict, indices, name, mySeq, segmentSeq, str(maxMismatches), str(numberOutputted), cas9Site,
              totalCas9)
    print("Finished. Check " + name + ".txt")


# Will open and save the fasta file along with the name
def loadFasta():
    # Opening the input file, with exception handling.
    try:
        fullSeq = ""
        segSeq = ""
        inFile = open(sys.argv[1])
        inFile2 = open(sys.argv[2])
        name = sys.argv[3]
        for seq in SeqIO.parse(inFile, "fasta"):
            fullSeq += seq.seq
        for seq in SeqIO.parse(inFile2, "fasta"):
            segSeq += seq.seq
    except IndexError:
        print("Usage: python3 CrispiSeg.py ReferenceGenome.fasta Segment.fasta" +
              " OutputName MaxMismatches NumberOutputted [NGG/NAG]")
        sys.exit()
    # except FileNotFoundError:
    #     print("Can't find " + sys.argv[1])
    #     sys.exit()
    return str(fullSeq), str(segSeq), name


# Loads the parameters from terminal for the max mismatches, the selection size,
# the number outputted and the Cas9 sites
def loadParameters():
    try:
        maxMismatches = sys.argv[4]
        if maxMismatches == '-':
            maxMismatches = 3
        numberOutputted = sys.argv[5]
        if numberOutputted == '-':
            numberOutputted = 100
        cas9Site = sys.argv[6]
        if cas9Site == '-':
            cas9Site = 'NGG'
    except IndexError:
        print("Usage: python3 CrispiSeg.py ReferenceGenome.fasta Segment.fasta" +
              " OutputName MaxMismatches NumberOutputted [NGG/NAG]")
        sys.exit()
    return int(maxMismatches), int(numberOutputted), cas9Site


# Given the sequence and total size, this will parse through the given file, returning 50 random
# Cas9 NGG sites, along with their indices
def parseFasta(segmentSeq, size, cas9Site):
    indices = {}
    matchDict = {}
    if cas9Site == 'NGG':
        pattern = ".{" + str(size) + "}GG"
    elif cas9Site == 'NAG':
        pattern = ".{" + str(size) + "}AG"
    else:
        print("Neither NGG or NAG selected")
        sys.exit()
    for kmer in regex.finditer(pattern, segmentSeq):
        matchDict[kmer.group(0)] = 0
        indices[kmer.group(0)] = str(kmer.start()) + " " + str(kmer.end())
    return matchDict, indices, str(len(matchDict))


# Given the sequence, dictionary with NGG sites, and a particular NGG site, will check for exact matches
# If it finds any, it will delete the NGG site out of the dictionary, since it's not worth checking anymore
# Then it returns the dictionary
def parseOffTarget(mySeq, myDict, kmer):
    total = (len(regex.findall(kmer, mySeq)) - 1) * 2  # arbitrary
    if total > 10:
        del myDict[kmer]
    elif kmer in myDict:
        myDict[kmer] += total
    return myDict


# Will check for 1 off target mismatch
def parseOffTarget1(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    total = 0
    for bp in range(start, end):
        temp1 = temp + helperResEnd(start, bp, kmer)
        total += len(regex.findall(temp1, mySeq)) / numOffTargets
    if total > 2:
        del myDict[kmer]
    elif kmer in myDict:
        myDict[kmer] += total
    return myDict


# Will check for 2 off target mismatches
def parseOffTarget2(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    for bp in range(start, end - 1):
        temp2 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget1(mySeq, myDict, bp + 1, end, temp2, kmer, numOffTargets)
    return myDict


# Will check for 3 off target mismatches
def parseOffTarget3(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget2(mySeq, myDict, bp + 1, end, temp1, kmer, numOffTargets)
    return myDict


# Will check for 4 off target mismatches
def parseOffTarget4(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget3(mySeq, myDict, bp + 1, end, temp1, kmer, numOffTargets)
    return myDict


# Will check for 5 off target mismatches
def parseOffTarget5(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget4(mySeq, myDict, bp + 1, end, temp1, kmer, numOffTargets)
    return myDict


# Will check for 6 off target mismatches
def parseOffTarget6(mySeq, myDict, start, end, temp, kmer, numOffTargets):
    for bp in range(start, end - 2):
        temp1 = temp + helperRes(start, bp, kmer)
        myDict = parseOffTarget5(mySeq, myDict, bp + 1, end, temp1, kmer, numOffTargets)
    return myDict


# Will check for 7 off target mismatches
def parseOffTarget7(mySeq, myDict, start, end, kmer, numOffTargets):
    for bp in range(start, end - 2):
        temp = helperRes(start, bp, kmer)
        myDict = parseOffTarget6(mySeq, myDict, bp + 1, end, temp, kmer, numOffTargets)
    return myDict


# Given the dictionary full of NGG sites and the total sequence, this will go through and check for up to
# totalOffStandard mismatches in each of the NGG sites, updating the percentage done as it goes
def callOnParsers(myDict, mySeq, maxMismatches):
    counter = 0
    print("Number of Cas9 Sites:" + str(len(myDict)))
    print("Status: ")
    myDictCopy = myDict.copy()
    if maxMismatches > 7:
        print("Max Mismatches can't be greater than 7.")
        sys.exit()
    for kmer in myDictCopy:
        totalOff = 0
        totalOffStandard = maxMismatches
        while totalOff < totalOffStandard + 1:
            if totalOff == 7:
                myDict = parseOffTarget7(mySeq, myDict, 0, len(kmer) - 2, kmer, totalOff)
            elif totalOff == 6:
                myDict = parseOffTarget6(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 5:
                myDict = parseOffTarget5(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 4:
                myDict = parseOffTarget4(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 3:
                myDict = parseOffTarget3(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 2:
                myDict = parseOffTarget2(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 1:
                myDict = parseOffTarget1(mySeq, myDict, 0, len(kmer) - 2, "", kmer, totalOff)
            elif totalOff == 0:
                myDict = parseOffTarget(mySeq, myDict, kmer)
            totalOff += 1
            percent = (1 / (len(myDictCopy)*2)) * totalOff / totalOffStandard
            percent = ((counter / len(myDictCopy)) + percent) * 100
            percent = round(percent, 2)
            print(str(percent) + "% ")
            if kmer not in myDict:
                break
        counter += 1
    return myDict


# Given a start, an end, and a particular kmer in question, will return a regex expression looking for
# the site, and replacing the last bp with a mismatch
def helperRes(start, bp, kmer):
    return kmer[start:bp] + regExDefiner(kmer[bp])


# The end of the regular expression specified by helperRes
def helperResEnd(start, bp, kmer):
    return helperRes(start, bp, kmer) + kmer[bp + 1:]


# Given a base pair, will return a regular expression for a mismatch
def regExDefiner(char):
    if char == 'A':
        return '[T,C,G]'
    if char == 'T':
        return '[A,C,G]'
    if char == 'C':
        return '[A,T,G]'
    if char == 'G':
        return '[A,T,C]'


# Writes the output file, given the dictionary, the indices, the name of the file and the sequence
def writeFile(myDict, indicies, name, mySeq, segSeq, maxMismatches, numberOutputted, cas9Site, totalCas9):
    counter = 0
    sortDict = sorted(myDict, key = myDict.get)
    fileTemp = open(name + ".txt", 'w+')
    fileTemp.write("Crispi -- Version 66.0123 Beta \n")
    fileTemp.write("Note: \n-->Indices are exclusive. "
                   + "\n-->The lower the score, the better the site, aka the less off-target probability\n\n")
    fileTemp.write("-- Parameters --\n")
    fileTemp.write("Total Genome Length: " + str(len(mySeq)) + "\n")
    fileTemp.write("Total Target Sequence Length: " + str(len(segSeq)) + "\n")
    fileTemp.write("Max Mismatches: " + maxMismatches + ", Cas9 Site: " + cas9Site + ", Total " +
                   cas9Site + " Sites: " + totalCas9 + "\n\n")
    fileTemp.write("Sequence (start, end) Score: --\n")
    for kmer in sortDict:
        counter += 1
        fileTemp.write(str(kmer) + " ( " + indicies[kmer] + " )" + " Score: " + str(myDict[kmer]) + "\n")
        if counter == int(numberOutputted):
            break
    if counter < int(numberOutputted):
        print("Not enough cas9 sites, outputted: " + str(counter))


main()
