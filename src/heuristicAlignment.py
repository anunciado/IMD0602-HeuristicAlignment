# -*- coding: utf-8 -*-

from __future__ import division
import sys
import gzip
import time
import operator

"""
Heuristic Aligner 1.0 by Luís Eduardo.

This script consists in the implementation of a heuristic aligner. The script
read a query.fa.gz file and align the sequences with those included in the
database.fa.gz file. At the end of the routine return the 10 database sequences
that best align to each query.fa.gz sequence.

Example:
    To run the script given the files query.fa.gz and database.fa.gz just:

        $ python heuristicAlignment.py query.fa.gz database.fa.gz
"""

def show(hashIdentity, time):
    """
    This is a method that a dict with DNA bases alignments, show the alignments.

    Args:
        hashIdentity (dict): A dict with DNA bases alignments.
        time (float): A float that represents the time of execution of the alignment.

    """
    print ("Heuristic Aligner 1.0 by Luís Eduardo.")
    print ("Time since submission: %.2f" %time + "s" + "\n")
    for x, (key, value) in enumerate(hashIdentity.items()):
        print ("Sequence Name Query: " + key.split(" ", 1)[1].rstrip('\n'))
        print ("Sequence ID: " + key.split(' ')[0])
        print ("Query Length: " + str(hashIdentity[key][0][3]) + "\n")
        print ("Sequences producing significant alignments:" + "\n")
        for y in value:
            print ("Sequence Name: " + y[0].split(" ", 1)[1].rstrip('\n'))
            print ("Sequence ID: " + y[0].split(' ')[0])
            print ("Length: " + str(y[4]))
            print ("Score: " + str(y[1]))
            print ("Identities: " + str(len(y[2])-y[2].count('-')) + "/" + str(len(y[2])) + "(" + "{:.0%}".format((len(y[2])-y[2].count('-'))/len(y[2])) + ")")
            print ("Sequence Alignment: " + "\n")
            print (y[2] + "\n")

def score(value1, value2):
    """
    This is a method that given two DNA bases parameters returns a values for
    match, missmatch and indel based on some alignment table of the parameters
    for alignments.

    Args:
        value1 (str): A DNA base of the first sequence.
        value2 (str): A DNA base of the second sequence.

    Returns:
        int: A values for match, missmatch or indel.

    """
    dictionary = {'AA': [5], 'CC': [5], 'GG': [5], 'TT': [5], 'AC': [4], 'CA': [4], 'AT': [3], 'TA': [3], 'AG' : [2], 'TG' : [2], 'TC' : [2], 'GC': [1], 'GA' : [2], 'GT' : [2], 'CT' : [2], 'CG': [1], 'A-' : [-2], 'C-' : [-2], 'G-' : [-2], 'T-' : [-2], '-A' : [-2], '-C' : [-2], '-G' : [-2], '-T' : [-2], '--' : [-2]}
    if (value1 + value2) in dictionary:
        return dictionary[value1 + value2][0]
    else:
        return -2

def identity(hashIdentity, labelsI, labelsQ, score, sequence, lengthQuery, lengthData):
    """
    This is a method that given a hash table with the best alignments in tuples,
    add a new tuple (where the first value is the sequence name aligned the
    table key and the second value is the alignment score) and return the best
    10 alignments.

    Args:
        hashIdentity (dict): A hash table of the best alignments.
        labelsI (str): The name of the reference sequence.
        labelsQ (str): The name of the database sequence aligned to the
        reference sequence.
        score (int): A value relative to the sequence alignment.
        sequence (str): A sequence alignmented with a database.
        lengthQuery (int): A value relative to the length of the query sequence.
        lengthData (int): A value relative to the length of the database sequence.

    Returns:
        dict: A dictionary where the keys are the names of the reference
        sequences and their values are the 10 best alignments in tuples
        format, where the first value is the sequence name aligned the table
        key and the second value is the alignment score

    """
    if labelsQ in hashIdentity.keys():
        hashIdentity[labelsQ].append((labelsI, score, sequence, lengthQuery, lengthData))
    else:
        hashIdentity[labelsQ] = [(labelsI, score, sequence, lengthQuery, lengthData)]
    for key, value in hashIdentity.items():
        hashIdentity[key] = sorted(value, key=lambda x: x[1], reverse=True)[:10]
    return hashIdentity

def localAlignment(sequenceDNA1, sequenceDNA2):
    """
    This is a method that does the local alignment of two given sequences and
    returns an arbitrary value of the alignment of the two sequences.

    Args:
        sequenceDNA1 (str): The reference sequence.
        sequenceDNA2 (str): The database sequence that will be aligned with
        the reference sequence.

    Returns:
        int: An arbitrary value of the local alignment of the two sequences.
        str: The best sequence of the local alignment of the two sequences.

    """
    sequenceDNA3 = ''
    matrix = [[0 for x in range(len(sequenceDNA2))] for y in range(len(sequenceDNA1))]
    for x, value1 in enumerate(sequenceDNA1):
        for y, value2 in enumerate(sequenceDNA2):
            lista = []
            if x == 0 and y > 0 and y == 0 and x > 0:
                matrix[x][y] = 0
            elif y > 0 and x > 0:
                lista.append(matrix[x-1][y-1] + score(value1, value2))
                lista.append(matrix[x-1][y] + score('-', value2))
                lista.append(matrix[x][y-1] + score(value1, '-'))
                matrix[x][y] = max(lista)
            del lista[:]

    aux1 = len(sequenceDNA1) - 1
    aux2 = len(sequenceDNA2) - 1

    for x, row in reversed(list(enumerate(matrix))):
        for y, element in reversed(list(enumerate(row))):
            if aux1 == x and aux2 == y and aux1 != 0 and aux2 != 0:
                if matrix[x][y] == matrix[x-1][y-1] + score(sequenceDNA1[x], sequenceDNA2[y]):
                    aux1 = x-1
                    aux2 = y-1
                    sequenceDNA3 += sequenceDNA1[x]
                elif matrix[x][y] == matrix[x-1][y] + score(sequenceDNA1[x-1], '-'):
                    aux1 = x - 1
                    sequenceDNA3 += "-"
                elif matrix[x][y] == matrix[x][y-1] + score('-', sequenceDNA2[y-1]):
                    aux2 = y - 1
                    sequenceDNA3 += "-"

    return max([sublist[-1] for sublist in matrix]), sequenceDNA3[::-1]

def readQuery():
    """
    This is a method that open and read the file 'query.fa.gz', after put in
    two list.

    Returns:
        list: A list of reference sequences.
        list: A list of the reference sequences names.
        float: A float that represents the time of execution of the function.

    """
    start = time.time()
    queries = []
    labels = []
    with open(sys.argv[1], 'r') as infile:
        novo = ''
        for line in infile:
            if(line[0]=='>'):
                labels.append(line[1:])
                if(novo):
                    queries.append(novo)
                    novo = ''
            else:
                novo += line.strip('\n')
    end = time.time()
    return queries, labels, (end - start)

def readDatabase(queries, labelsQ):
    """
    This is the main script method that open and read the file 'database.fa.gz'
    and makes alignments of all sequences with the reference sequences and
    return the best 10 alignments with every reference sequenceDNA1.

    Args:
        queries (list): A list of reference sequences.
        labelsQ (list): A list of the reference sequences names.

    Returns:
        dict: A dictionary where the keys are the names of the reference
        sequences and their values are the 10 best alignments in tuples
        format, where the first value is the sequence name aligned the table
        key and the second value is the alignment score
        float: A float that represents the time of execution of the function.

    """
    start = time.time()
    hashQuery = {}
    hashIdentity = {}

    with open(sys.argv[2], 'r') as infile:
        labelsI = ''
        novo = ''
        dictionaryI = {}
        for line in infile:
            if(line[0]=='>'):
                labelsI = line[1:]
                if(novo):
                    for w, query in enumerate(queries):
                        longestSequence = [0, 0, 0, 0, 0]
                        secondlongestSequence = [0, 0, 0, 0, 0]
                        for x in range(0, len(query) - 10):
                            key = query[x : x + 11]
                            if key in hashQuery:
                                hashQuery[key].append(x + 1)
                            else:
                                hashQuery[key] = [x + 1]
                        for y in range(0, len(novo) - 10):
                            key = novo[y : y + 11]
                            if key in hashQuery:
                                if (y - 1) in dictionaryI:
                                    dictionaryI[y] = [dictionaryI[y - 1][0], dictionaryI[y-1][1] + 1, dictionaryI[y - 1][2]]
                                    del dictionaryI[y-1]
                                else:
                                    dictionaryI[y] = [y, 1, hashQuery[key][0]]
                        # dictionaryI[last position in novo] = [first position in novo, number of sets de matchs, first position in query]
                        for key, item in dictionaryI.items():
                            if dictionaryI[key][1] > longestSequence[2]:
                                secondlongestSequence = longestSequence
                                longestSequence = [key, dictionaryI[key][0], dictionaryI[key][1], dictionaryI[key][2]]
                            elif dictionaryI[key][1] > secondlongestSequence[2]:
                                secondlongestSequence = [key, dictionaryI[key][0], dictionaryI[key][1], dictionaryI[key][2]]
                        if (longestSequence[2] > 90) or (secondlongestSequence[2] > 3):
                            if secondlongestSequence[2] == 0:
                                hashIdentity = identity(hashIdentity, labelsI, labelsQ[w], (longestSequence[2]+10)*5, query[longestSequence[3]:longestSequence[3]+longestSequence[2]+10], len(query), len(novo))
                            else:
                                if longestSequence[0] < secondlongestSequence[1] and len(query[longestSequence[3]:secondlongestSequence[3]+secondlongestSequence[2]+10]) > 1 and len(novo[longestSequence[1]:secondlongestSequence[0]+10]) > 1:
                                    score, sequence = localAlignment(query[longestSequence[3]:secondlongestSequence[3]+secondlongestSequence[2]+10],novo[longestSequence[1]:secondlongestSequence[0]+10])
                                    hashIdentity = identity(hashIdentity, labelsI, labelsQ[w], score, sequence, len(query), len(novo))
                                else:
                                    if len(query[secondlongestSequence[3]:longestSequence[3]+secondlongestSequence[2]+10]) > 1 and len(novo[longestSequence[1]:secondlongestSequence[0]+10]) > 1:
                                        score, sequence = localAlignment(query[secondlongestSequence[3]:longestSequence[3]+secondlongestSequence[2]+10],novo[longestSequence[1]:secondlongestSequence[0]+10])
                                        hashIdentity = identity(hashIdentity, labelsI, labelsQ[w], score, sequence, len(query), len(novo))
                            hashQuery.clear()
                        dictionaryI.clear()
                novo = ''
            else:
                novo += line.strip('\n')
    end = time.time()
    return hashIdentity, (end - start)

def main():
    """
    This is a method main only used to run all the method for this script.

    """
    queries, labelsQ, time1 = readQuery()
    hashIdentity, time2 = readDatabase(queries, labelsQ)
    show(hashIdentity, time1 + time2)

if __name__ == "__main__":
    main()
