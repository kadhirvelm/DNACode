#! /usr/bin/python
from __future__ import division
import operator
import sys
import re

from Bio import SeqIO

import math
import math

__author__ = 'Kadhir'


def reading(readfile, sequenceFiles):  # outputs a parsed Stockholm file that has the [identifier, sequence]
    for line in readfile:
        if line.startswith("//"):
            break
        if not line.startswith("#"):
            line = line.replace(" ", "")
            tempname = re.search('[A-Z]+[0-9]+.+[0-9]+.[0-9]+', line)
            tempseq = re.search(".+", line[len(tempname.group())::])
            sequenceFiles.append(tempname.group())
            sequenceFiles.append(tempseq.group())
    return sequenceFiles


def counter(sequenceFiles):  # outputs a dictionary filled with the Shannon entropy of each column
    colDictionary = {}
    colnums = len(sequenceFiles[1])
    for col in range(colnums):
        tempentropy = probcalculator(sequenceFiles, col)
        colDictionary[col] = tempentropy
    return colDictionary


def counter2(sequenceFiles):  # outputs a dictionary filled with the Shannon entropy of two columns
    colDictionary = {}
    colnums = len(sequenceFiles[1])
    for col1 in range(colnums - 1):
        prob1 = probcalculator(sequenceFiles, col1)
        for col2 in range(col1 + 1, colnums):
            prob2 = probcalculator(sequenceFiles, col2)
            prob12 = prob2calculator(sequenceFiles, col1, col2)
            tempmutualinfo = prob1 + prob2 - prob12
            colDictionary[col1, col2] = tempmutualinfo
    return colDictionary


def probcalculator(sequenceFiles, col):  # outputs the Shannon entropy at column col
    tempprob = 0
    numrows = len(sequenceFiles[1::2])
    letters = countAll(sequenceFiles, col)
    for temp in letters:
        if temp > 0:
            tempprob += temp / numrows * math.log(temp / numrows, 2)  # calculates probability and entropy at once
    return -tempprob


def prob2calculator(sequenceFiles, col1, col2):  # outputs Shannon entropy of col1 and col2
    tempprob = 0
    numrows = len(sequenceFiles[1::2])
    probtemp = countAll2(sequenceFiles, col1, col2)
    for temp in probtemp.values():
        tempprob += temp / numrows * math.log(temp / numrows, 2)  # calculates probability and entropy at once
    return -tempprob


def countAll(sequenceFiles, col):  # send in the sequencefile, and it'll return a list with nucleotide numbers
    letters = [0, 0, 0, 0, 0]
    for segment in sequenceFiles[1::2]:
        templetter = segment[col]
        if templetter == 'A':
            letters[0] += 1
        elif templetter == 'C':
            letters[1] += 1
        elif templetter == 'U':
            letters[2] += 1
        elif templetter == 'G':
            letters[3] += 1
        else:
            letters[4] += 1
    return letters


def countAll2(sequenceFiles, col1, col2):  # outputs a dictionary filled with the different variations of nucleotides
    dict_letters = {}  # + blanks and how frequently they pop up
    for segment in sequenceFiles[1::2]:
        templetter1 = segment[col1]
        templetter2 = segment[col2]
        try:
            dict_letters[templetter1, templetter2] += 1
        except KeyError:
            dict_letters[templetter1, templetter2] = 1
    return dict_letters


def sorter(colDictionary, rang):  # prints the lowest rang values in colDictionary
    temp = sorted(colDictionary.items(), key = operator.itemgetter(1))
    for tempP in range(rang):
        print(temp[tempP][0])


def sorter2(colDictionary, rang):  # prints the highest rang values in colDictionary
    temp = sorted(colDictionary.items(), key = operator.itemgetter(1))
    total = len(temp) - 1
    for tempP in range(rang):
        print(temp[total - tempP][0])


def main():
    try:
        readfile = open(sys.argv[1], "r")  # read the input file
        sequenceFiles = []
        sequenceFiles = reading(readfile, sequenceFiles)  # outputs a list with [identifier, sequence]+
        colDictionary = counter(sequenceFiles)  # outputs a dictionary with each column identifier, entropy value
        sorter(colDictionary, 10)  # prints the lowest 10 entropy values
        colDictionary = counter2(sequenceFiles)  # outputs a dictionary with each column pair identifier and entropy val
        sorter2(colDictionary, 50)  # prints the highest 50 entropy values
    except IOError:  # If there's an invalid input file
        print("Error: Can't find " + sys.argv[1])
    except IndexError:  # If there's no input file
        print("Error: Please enter a file")


main()
