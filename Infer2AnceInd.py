import numpy
import math
import random
import time
from TransitionsInd import *
from emissions2Ind import *

def InferAnceInd(Nc, Ng, Recomb, PositionFile, ChildrenFile, Ge, AFAFile, AFBFile, dataz):
    #Read PositionFile:
    datag = Readfile(PositionFile)

    position = [0.999999999]

    for h in range(1,len(datag)):
        position.append(float(datag[h][0])-float(datag[h-1][0]))

    #Read Alleles Files:
    A = Readfilefloat(AFAFile)
    B = Readfilefloat(AFBFile)


    #Generate States:
    states = []
    number = 2
    for i in range(2**number):
        states.append(bin(i)[2:].zfill(number))
    sn = 2**number

    #Generate Recombinations:
    probability = []

    for i in range(len(position)):
        RF = position[i] * Ng * Recomb
        if(RF == 0.0):
            RF = 0.000000001
        elif(RF >= 1.0):
            RF = 0.999999999
        probability.append(RF)

    #Define model:
    model = 1
    Comodel = 1

    #Read Phasing and Recombinations:
    dataR = []
    for i in range(len(dataz)):
        dataR.append(dataz[i][:Nc])

    #Read Observations:
    datac = Readfile(ChildrenFile)
    observations = []

    for i in range(len(datac[0])):
        datax = datac[0][i]
        for j in range(1, Nc):
            datax = datax + datac[j][i]
        observations.append(datax)

    #Generate Alphabet:
    alphabet = []
    GenoC = Nc
    for i in range(2**(GenoC)):
        alphabet.append(bin(i)[2:].zfill(GenoC))
    alphabet_index = {}
    for i, alphabet in enumerate(alphabet):
        alphabet_index[alphabet] = i
    alphabet = []
    for i in range(2**GenoC):
        alphabet.append(bin(i)[2:].zfill(GenoC))

    EC = numpy.zeros((len(states), len(alphabet)), dtype=object)
    for i in range(len(states)):
        for j in range(len(alphabet)):
            EC[i][j] = states[i]+alphabet[j]

    #Forward Algorithm:
    forward = numpy.zeros((len(observations),sn), dtype=object)

    emission = emI2(EC,len(states),alphabet_index[observations[0]],A[0],B[0], Nc, dataR[0], Ge)

    for x in range(sn):
        forward[0][x] = math.log(emission[x]/sn)
    for num in range(len(observations)-1):

        emission = emI2(EC,len(states),alphabet_index[observations[num+1]],A[num+1],B[num+1], Nc, dataR[num+1], Ge)

        forwardline = forward[num]

        T2P = forwardline.reshape(-1,1)

        TP = forwardbackward(probability[num+1], GenoC, Nc, number, model, T2P, Comodel, 0,0,0,0)

        for x in range(sn):
            emit = math.log(emission[x])
            forward[num+1][x] = TP[0][x] + emit


    #Backward Algorithm:
    backward = numpy.zeros((len(observations),sn), dtype=object)

    weight = numpy.zeros((1,sn), dtype=object)
    for x in range(sn):
        backward[-1][x] = math.log(1.0)
    for num in range(-1,-len(observations),-1):
        emission = emI2(EC,len(states),alphabet_index[observations[num]],A[num],B[num], Nc, dataR[num], Ge)

        for y in range(sn):
            weight[0][y] = math.log(emission[y])+backward[num][y]
        backwardline = weight[0]
    
        T2P = backwardline.reshape(-1,1)

        TP = forwardbackward(probability[num], GenoC, Nc, number, model, T2P, Comodel, 0,0,0,0)

        for x in range(sn):
            backward[num-1][x] = TP[0][x]

    #Sequence:
    sequence = []
    for i in range(len(observations)):
        s = forward[i][0]+backward[i][0]
        for j in range(1,len(states)):
            s = logA(forward[i][j]+backward[i][j],s)
        sequence.append(s)

    #Compute Conditional Probabilities:
    condit_prob = [ [0.0 for x in range(len(states))] for x in range(len(observations))]
    for i in range(len(observations)):
        for j in range(len(states)):
            condit_prob[i][j] = forward[i][j]+backward[i][j]-sequence[i]

    state0 = []
    for prob in condit_prob:
        state0.append(states[numpy.argmax(prob)])

    return state0