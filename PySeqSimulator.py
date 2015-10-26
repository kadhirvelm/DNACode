#! /usr/local/bin/python

import random
import sys

from Bio import SeqIO

__author__ = 'Kadhir'


def parsing(my_seq, compList):  # will calculate composition of my_seq and store it in compList
    totalSeq = 0
    for sequence in SeqIO.parse(my_seq, "fasta"):  # begins parsing through the fasta file my_seq
        totalSeq += 1  # adds 1 to the total sequence count
        for letter in sequence.seq:  # parses through and appends compList accordingly
            if letter == 'A':
                compList[0] += 1
            elif letter == 'C':
                compList[1] += 1
            elif letter == 'T':
                compList[2] += 1
            elif letter == 'G':
                compList[3] += 1
    for index in range(4):  # change compLists to averages, since it's currently just totals
        compList[index] = compList[index] / totalSeq
    return compList, sum(compList), totalSeq


def createSequence(length, output, comp, userMotifs):  # creates a random sequence base on length and composition
    lenTemp = 3
    userMotifs = userMotifs[1:]  # takes out the extraneous userMotif
    usermindex = len(userMotifs) - 1  # index of which userMotif to add defined
    random.shuffle(userMotifs)  # shuffles them around to add randomness
    temp = redefineLength(length, comp)
    length = int(temp[0])
    comp = temp[1]
    compListTemp = list(comp)
    hasToBeFixed = False
    while lenTemp <= length - 3:  # leave space for stop codon
        if (compListTemp[4] >= 1) & (
                    (random.randint(0, int(length / (usermindex + 2))) == 0) | (lenTemp > length - 10)):
            nextletter = 4  # condition for adding a motif to the sequence
        else:
            nextletter, hasToBeFixed = chooseRant(lenTemp, output)  # if not adding motif, then add a random letter
        if (compListTemp[nextletter] > -2) | (
                hasToBeFixed):  # if the particular letter needs more letters printed, or has to be printed
            compListTemp[nextletter] -= 1  # to prevent a premature stop codon, then print it
            lenTemp += 1
            if hasToBeFixed:
                hasToBeFixed = False  # boolean for needing to fix the premature stop codons
            if nextletter == 4:
                output += userMotifs[usermindex]
                usermindex -= 1
            else:
                output = addOutput(nextletter, output)  # go find the appropriate letter to print
    return output


def redefineLength(length, comp):  # redefines the length according to a Gaussian distribution
    for index in range(4):  # don't want other method, because we need float, not rounded numbers here
        comp[index] /= float(length)  # convert comp into percentages
    length = random.gauss(length, length / 4)  # random.gauss distribution for the length
    for index in range(4):
        comp[index] *= length  # to ensure that the percentages still remain the same, but works with the earlier code
    return length, comp


def addOutput(nextletter, output):  # Given what the next letter needs to be, will add the appropriate letter to output
    if nextletter == 0:
        output += "A"
    elif nextletter == 1:
        output += "C"
    elif nextletter == 2:
        output += "T"
    elif nextletter == 3:
        output += "G"
    return output


def chooseRant(lenTemp, output):  # given the current output and current length, will check and prevent premature
    if output[lenTemp - 2:lenTemp] == "TA":  # stop codons
        return random.randint(1, 2), True
    elif output[lenTemp - 2:lenTemp] == "TG":
        return random.randint(1, 3), True
    else:
        return random.randint(0, 3), False


def startOutput(compList, usermotifs):  # given numbers of nucleotides and motifs, will create the sequence output
    compList.append(len(usermotifs) - 1)  # first append the number of motifs to the end
    output = createSequence(sum(compList), "ATG", compList, usermotifs)  # create the sequence
    for leftover in range(len(output) % 3):  # make sure the output is in the proper
        nextletter = chooseRant(len(output), output)  # choose random numbers to fill in the places and make it 3mers
        addOutput(nextletter, output)
    output += "TAG"  # add the stop codon to the end
    return output


def extractCompList(readfile):  # given a file, this method will parse out the following:
    compList = [0, 0, 0, 0]  # the composition of nucleotides
    numNucleotides = 0  # total number of nucleotides
    numSequences = 0  # number of sequences
    userspecified = []  # the user specific motifs
    userindex = 0
    for line in readfile:
        if compList == [0, 0, 0, 0]:  # get the composition
            compList = list(extractList(compList, line))
        elif numNucleotides == 0:  # get the nucleotide number
            tempIndex = findSpace(line)
            numNucleotides = float(line[0:tempIndex])
        elif numSequences == 0:  # get the sequences
            tempIndex = findSpace(line)
            numSequences = int(line[0:tempIndex])
        else:  # get the user motifs
            userspecified.append(line[:-1])
            userindex += 1

    for index in range(4):  # convert the percentage of nucleotides into numbers
        compList[index] = compList[index] / 100 * numNucleotides

    return compList, numNucleotides, numSequences, userspecified


def extractList(compList, line):  # pulls out the numbers out of a line in readfile
    prevIndex = 0
    currIndex = 0
    counter = 0
    findInts = False
    for letters in line:  # parses through each letter in the line
        if letters == '[':
            findInts = True
            prevIndex = counter + 1
        elif (letters == " ") & findInts:
            compList[currIndex] = float(line[prevIndex:counter - 1])  # convert the percents to floats and add them
            prevIndex = counter + 1  # change the previous index to the current one
            currIndex += 1  # change the current index
            if currIndex == 4:
                break
        counter += 1  # move to the next step, increasing counter by one
    return compList


def findSpace(line):  # given a line, it'll return the index priot to the first space
    templine = line
    counter = 0
    tempindex = 0
    for num in templine:
        counter += 1
        if num == " ":
            tempindex = counter - 1
            break
    return tempindex


def adjustCompToPercent(compList):  # changes composition to percent instead of number
    total = float(sum(compList))
    newlist = map(float, compList)  # experimenting with map function, but change compList to floats
    for index in range(4):
        newlist[index] = round(compList[index] / total * 100, 1)  # only one decimal place for the percent
    return newlist


def main():  # main function, does all of the parsing of the terminal commands
    global complist
    try:
        if sys.argv[1] == "--calc":  # if there's a calc present, do the following:
            fastafile = sys.argv[2]  # get the file out
            complist = parsing(fastafile, [0, 0, 0, 0])  # parse it, getting the composition out
            output = startOutput(complist[0], "")  # output using the composition data
            print(output)  # print to the terminal
        elif sys.argv[1] == "--load":  # if there's a load instead of a --calc
            readfile = open(sys.argv[2], "r")  # read the param file
            complist = extractCompList(readfile)  # extract the needed parameters out of it
            output = startOutput(complist[0], complist[3])  # send the composition to get created
            print(output)  # print the whole thing out
    except IndexError:  # if there's nothing proper there
        if len(sys.argv) <= 1:  # and there are no arguments
            complist = [249, 250, 249, 249]  # A,C,T,G totalling up to 1 kilobase - gets overwritten due to gauss
            usermotifs = ["Ignore"]  # trash usermotifs, to prevent bugs mostly
            output = startOutput(complist, usermotifs)
            print(output)  # print it to terminal
        else:
            print("Error, wrong input!")  # if the input is just wrong, print an error message
    except IOError:
        print("Not a valid file!")  # if the file is incorrect, print an error message
    try:
        if (sys.argv[1] == "--h") | (sys.argv[1] == "--help"):  # if the user wants help, print instructions
            print("Use '--calc fasta.fasta' to create a random sequence based on fasta.fasta")
            print("Or Use '--load params' to create a random sequence based on the parameters in params")
            print("Follow either command above with '--save file' to save the logistics data to file")
        if sys.argv[3] == "--save":  # if the user wants to save the composition and other data
            file = open(sys.argv[4], "w")  # open/create the necessary file and set it to write
            file.write(str(adjustCompToPercent(complist[0])))  # change the composition to percent to write to file
            file.write(" # Percent of: A,C,T,G, # of user defined motifs (<=0 for no motifs) \n")
            file.write(str(int(sum(complist[0]))) + " #Number of nucleotides " + "\n")
            file.write(str(complist[2]) + " #Number of sequences: " + "\n")
            file.write("User specified motifs, each one on its own line:\n")
            file.close()  # make sure to close the file
            # EXAMPLE OF OUTPUT OF WRITTEN FILE:
            # [27.1, 30.0, 17.7, 25.2, -1.0] # Percent of: A,C,T,G, # of user defined motifs (<=0 for no motifs)
            # 2648 #Number of nucleotides
            # 2 #Number of sequences:
            # User specified motifs, each one on its own line:
    except IndexError:
        pass  # if there's no --h or --save, that's okay, just pass it along


main()
