import numpy
import math
import random
from emissions1 import *
from Transitions import *


def InferRecom(Nc, Ng, Recomb, PositionFile, ChildrenFile, PPE, Ge):

    #Read PositionFile:
    datag = Readfile(PositionFile)

    position = [0.999999999]

    for h in range(1,len(datag)):
        position.append(float(datag[h][0])-float(datag[h-1][0]))

    #Generate States:
    states = []
    number = Nc + 2*Nc + 2*2
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


    #Generate Phasing Errors:
    phasingerror = []
    for i in range(len(position)):
        PE = position[i] * Ng * PPE
        if(PE == 0.0):
            PE = 0.000000001
        elif(PE >= 1.0):
            PE = 0.999999999
        phasingerror.append(PE)

    #Define model and Comodel:
    model = 0
    Comodel = 0

    #Read Observations:
    datac = Readfile(ChildrenFile)
    observations = []

    for i in range(len(datac[0])):
        datax = datac[0][i] + datac[1][i]
        for j in range(2, Nc*2):
            datax = datax + datac[j][i]
        observations.append(datax)

    #Generate Alphabet:
    alphabet = []
    GenoC = 2*Nc
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

    #Generate Emission Matrix:
    emission = em1(EC,Ge,len(states),len(alphabet), Nc)


    #Forward Algorithm:
    forward = numpy.zeros((len(observations),sn), dtype=object)
    for x in range(sn):
        forward[0][x] = math.log(emission[x][alphabet_index[observations[0]]]/sn)
    for num in range(len(observations)-1):
        forwardline = forward[num]

        T2P = forwardline.reshape(-1,1)

        TP = forwardbackward(probability[num+1], phasingerror[num+1], GenoC, Nc, number, model, T2P, Comodel, 0, 0,0,0,0,0,0,0)

        for x in range(sn):
            emit = math.log(emission[x][alphabet_index[observations[num+1]]])
            forward[num+1][x] = TP[0][x] + emit

    #Backward Algorithm:
    backward = numpy.zeros((len(observations),sn), dtype=object)
    weight = numpy.zeros((1,sn), dtype=object)
    for x in range(sn):
        backward[-1][x] = math.log(1.0)
    for num in range(-1,-len(observations),-1):
        for y in range(sn):
            weight[0][y] = math.log(emission[y][alphabet_index[observations[num]]])+backward[num][y]
        backwardline = weight[0]
    
        T2P = backwardline.reshape(-1,1)

        TP = forwardbackward(probability[num], phasingerror[num], GenoC, Nc, number, model, T2P, Comodel, 0, 0,0,0,0,0,0,0)

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


    max = []
    for prob in condit_prob:
        max.append(numpy.where(prob==numpy.max(prob)))

    #Find Orders:
    state0 = ['0.0' for x in range(len(observations))]
    fixed = []
    for x in range(len(observations)):
        if(len(max[x][0]) == 1):
            state0[x] = states[max[x][0][0]]
            fixed.append(x)

    if len(fixed) != 0:
        if(fixed[0] != 0):
            for i in range(fixed[0]-1, -1, -1):
                tube =[]
                for j in range(len(max[i][0])):
                    pair = 0
                    for x in range(GenoC+Nc):
                        if(states[max[i][0][j]][x] == state0[i+1][x]):
                            pair += 1
                    tube.append(pair)
                maxtube = numpy.argmax(tube)
                state0[i] = states[max[i][0][maxtube]]

        for i in range(len(fixed)-1):
            for j in range(fixed[i]+1,fixed[i+1]):
                tube =[]
                for x in range(len(max[j][0])):
                    pair = 0
                    for y in range(GenoC+Nc):
                        if(states[max[j][0][x]][y] == state0[j-1][y]):
                            pair += 1
                    tube.append(pair)
                maxtube = numpy.argmax(tube)
                state0[j] = states[max[j][0][maxtube]]

        for i in range(fixed[-1],len(observations)):
            tube = []
            for j in range(len(max[i][0])):
                pair = 0
                for x in range(GenoC+Nc):
                    if(states[max[i][0][j]][x] == state0[i-1][x]):
                        pair += 1
                tube.append(pair)
            maxtube = numpy.argmax(tube)
            state0[i] = states[max[i][0][maxtube]]

    else:
        state0[0] = states[max[0][0][0]]
        for i in range(1,len(observations)):
            tube = []
            for j in range(len(max[i][0])):
                pair = 0
                for x in range(GenoC+Nc):
                    if(states[max[i][0][j]][x] == state0[i-1][x]):
                        pair += 1
                tube.append(pair)
            maxtube = numpy.argmax(tube)
            state0[i] = states[max[i][0][maxtube]]
    return state0, len(observations)