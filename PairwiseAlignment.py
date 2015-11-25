#! /usr/bin/python
import sys
import re

from Bio import SeqIO


# The letters being looked at
topLetters = []
# The letter scores
compiledDictionary = {}
# Gap score for the matrix
gapScore = 0
# The first sequence
firstSequence = ""
# The first name
firstName = ""
# The second sequence
secondSequence = ""
# The second name
secondName = ""
# final scored matrix
finalMatrix = []
# The first alignment
alignment1 = []
# The second alignment
alignment2 = []


# Taking in the string file, this will parse and retreive the two sequences
# spitting out errors if one of the sequences is blank
def parseStringFile(stringFile):
    counter = 0
    for sequence in SeqIO.parse(stringFile, "fasta"):
        if counter == 0:
            global firstSequence
            global firstName
            firstSequence = sequence.seq
            firstName = sequence.description
            counter = 1
        elif counter == 1:
            global secondSequence
            global secondName
            secondName = sequence.description
            secondSequence = sequence.seq
    if secondSequence == "":
        try:
            raise ValueError
        except ValueError:
            print("Error: Two sequences not present.")


# This will parse the matrixFile, creating the scoring dictionary and defining
# the gap score
def parseMatrixFile(matrixFile):
    global topLetters
    global compiledDictionary
    global gapScore
    counter1 = 0
    for line in matrixFile:
        temp = line.strip().split(" ")
        counter2 = 0
        for letters in temp:
            if (counter1 == 0):
                topLetters.append(letters.upper())
            elif (counter1 <= len(topLetters)):
                compiledDictionary[topLetters[counter1 - 1] + topLetters[counter2]] = int(letters)
            else:
                gapScore = int(letters)
            counter2 += 1
        counter1 += 1


# This will in take a matrix and check if the format is correct
def checkSimilarityMatrix(matrixFile):
    counter = 0
    step1 = 0
    step2 = 0
    step3 = 0
    for line in matrixFile:
        if counter == 0:
            re.match("[A-Z] *[A-Z]", line, re.IGNORECASE)
            length = len(line.split(" "))
            step1 += 1
        elif len(line.split(" ")) == length:
            re.match("-?[0-9]+ *-?[0-9]+", line)
            step2 += 1
        elif counter == length + 1:
            re.match("-?[0-9]+", line)
            step3 += 1
        else:
            raise NameError
        counter += 1
    checkSteps([step1, step2, step3])


# Given a step, will check if all the steps have taken place
def checkSteps(steps):
    for bools in range(0, 3):
        if steps[bools] == 0:
            raise NameError


# Check if the sequences have random letters
def checkStringFile():
    for letters in firstSequence:
        if checkInTop(letters):
            raise EnvironmentError
    for letters in secondSequence:
        if checkInTop(letters):
            raise EnvironmentError


# Loops through the top letters, checking if char is within it
def checkInTop(char):
    for letters in topLetters:
        if char == letters:
            return False
    return True


# Calculates the score matrix with the gap penalties
def calculateMatrix():
    array = [' ', '-']
    temparray = []
    global finalMatrix
    for letters in firstSequence:
        array.append(letters)
    finalMatrix.append(array)
    for letters in secondSequence:
        temparray = [letters]
        finalMatrix.append(temparray)
    for currCounter in range(len(array) - 1):
        finalMatrix[1].append(currCounter * gapScore)
    for currCounter in range(2, len(finalMatrix)):
        finalMatrix[currCounter].append((currCounter - 1) * gapScore)
    for currRow in range(2, len(finalMatrix)):
        for currCol in range(2, len(array)):
            valueTemp = calcValue(currRow, currCol)
            value = max(valueTemp)
            finalMatrix[currRow].append(value)
    return finalMatrix


# Given the current row and current col, will calculate the diagonal score
# And the horizontal and vertical scores
def calcValue(currRow, currCol):
    value = [0, 0, 0]
    char1 = finalMatrix[0][currCol]
    char2 = finalMatrix[currRow][0]
    value[0] = finalMatrix[currRow - 1][currCol - 1] + compiledDictionary.get(char1 + char2)
    value[1] = finalMatrix[currRow - 1][currCol] + gapScore
    value[2] = finalMatrix[currRow][currCol - 1] + gapScore
    return value


# Prints the finalMatrix
def printMatrix(finalMatrix):
    for line in finalMatrix:
        print(line)


# Writes the final alignment into a file called output
def printAlignments():
    firstSequence = ""
    secondSequence = ""
    for letters in reversed(alignment1):
        firstSequence += letters
    for letters in reversed(alignment2):
        secondSequence += letters
    fileTemp = open('output', 'w+')
    fileTemp.write(firstName + " aligned to " + secondName)
    fileTemp.write("\n")
    for letter in range(1, len(firstSequence) + 1):
        fileTemp.write(firstSequence[letter - 1])
        if letter % 80 == 0:
            fileTemp.write("\n")
            for letter in range(letter - 79, letter + 1):
                fileTemp.write(secondSequence[letter - 1])
            fileTemp.write("\n\n")
        elif letter == len(firstSequence):
            fileTemp.write("\n")
            if letter - 80 > 0:
                for letter in range((letter - letter % 80) + 1, len(firstSequence) + 1):
                    fileTemp.write(secondSequence[letter - 1])
            else:
                for letter in range(1, len(firstSequence) + 1):
                    fileTemp.write(secondSequence[letter - 1])


# Generates the alignment sequence
def generateSequence():
    determineSequence(len(finalMatrix) - 1, len(finalMatrix[0]) - 1)


# Determines the alignment sequence
def determineSequence(row, col):
    if (row == 1) & (col == 1):
        return
    elif row == 1:
        temp = horizontalAlign(row, col)
    elif col == 1:
        temp = verticalAlign(row, col)
    else:
        valueTemps = calcValue(row, col)
        tempInt = finalMatrix[row][col]
        if valueTemps[0] == tempInt:
            temp = diagonalAlign(row, col)
        elif valueTemps[1] == tempInt:
            temp = verticalAlign(row, col)
        elif valueTemps[2] == tempInt:
            temp = horizontalAlign(row, col)
        else:
            return
    determineSequence(temp[0], temp[1])


# Changes to the diagonal alignment given row and col
def diagonalAlign(row, col):
    alignment1.append(finalMatrix[0][col])
    alignment2.append(finalMatrix[row][0])
    return row - 1, col - 1


# Changes to the vertical alignment given row and col
def verticalAlign(row, col):
    alignment1.append(finalMatrix[0][1])
    alignment2.append(finalMatrix[row][0])
    return row - 1, col


# Changes to the horizontal alignment given row and col
def horizontalAlign(row, col):
    alignment1.append(finalMatrix[0][col])
    alignment2.append(finalMatrix[1][0])
    return row, col - 1


# The main function
def main():
    try:
        matrixFile = open(sys.argv[1], "r")  # read the input file
        checkSimilarityMatrix(matrixFile)
        matrixFile.seek(0)
        parseMatrixFile(matrixFile)
        stringFile = open(sys.argv[2], "r")
        parseStringFile(stringFile)
        checkStringFile()
        global secondSequence
        secondSequence = '-' + secondSequence
        calculateMatrix()
        generateSequence()
        printAlignments()
        fileTemp = open("output", "r")
        for line in fileTemp:
            print(line.strip())
    except EnvironmentError:
        print("Error: The sequences have characters not present in the similarity matrix.")
    except IOError:  # If there's an invalid input file
        print("Error: Can't find " + sys.argv[1] + " and/or " + sys.argv[2])
    except IndexError:  # If there's no input file
        print("Error: Please enter two files, a matrix file and a string file.")
    except NameError:
        print("Error: The Similarity Matrix is incorrectly written.")


main()
