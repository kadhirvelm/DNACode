#! /usr/bin/python
__author__ = 'Kadhir'

import sys


def reverseCompliment(a, rnapresent):  # method that when given a basepair 'a' and boolean rnaPresent
    if a == 'A':
        # it will spit out the complimentary strand
        if rnapresent:
            return 'U'
        else:
            return 'T'
    elif a == 'G':
        return 'C'
    elif a == 'C':
        return 'G'
    elif a == 'T':
        return 'A'
    else:
        return a


def printSequenceToTerminal(finallist, count, description):
    # finalList = list of base pairs, count = number of base pairs, description = description line
    print(description.rstrip() + ", %d bp" % count)  # prints description
    loopindex = len(finallist)
    while loopindex > 0:  # keep looping until nothing left to print
        if loopindex - 80 >= 0:  # if there are more than 80 characters left
            print(''.join(finallist[loopindex:loopindex - 80:-1]))
        else:  # less than 80 characters left
            print(''.join(finallist[loopindex::-1]))  # print the rest of it, backwards
        loopindex -= 80  # move to the next 80 indices

def main():
    global infile, description
    try:
        infile = open(sys.argv[1])  # open fasta file
    except IOError:
        print("Invalid File Name - Now Exiting")  # if invalid
        quit()
    try:
        rnapresent = sys.argv[2]  # look for rna
        if rnapresent.upper() == "RNA":
            rnapresent = True
    except IndexError:
        rnapresent = False  # set false if dna desired
    count = 0
    finalist = []  # base pairs array
    printdescription = False  # initializing that the description has not been printed yet
    for basepair in infile:  # looping through each line in the fasta file
        if not basepair.startswith('>'):
            for loop in basepair:  # loop through each base pair in each line in the fasta file
                if loop != '\n':  # don't want new lines, it messes up formatting
                    finalist.append(reverseCompliment(loop, rnapresent))
                    # go to reverseCompliment to get the compliment, don't know why I called it reverse
                    count += 1 #add one to the count of this particular sequence
        else:
            if printdescription:
                printSequenceToTerminal(finalist, count, description)
                del finalist[:]
                count = 0
            description = basepair  # if the line starts with '>', set it equal to the description
            printdescription = True
    printSequenceToTerminal(finalist, count, description)  # last sequence in the file printed
    infile.close()  # close the file
main()  # start the program by going to main
