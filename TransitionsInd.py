import numpy
import math

#Read Files
def Readfile(File):
	data=[]
	with open(File, 'r') as f:
		for x in f.readlines():
			x=x.strip('\n')
			x=x.split(' ')
			data.append(x)
	return data

def Readfilefloat(File):
	data=[]
	with open(File, 'r') as f:
		for x in f.readlines():
			x=x.strip('\n')
			x=x.split(' ')
			data.append(float(x[0]))
	return data

#Log Add between Matrix
def logadd(X,Y):
    Z = numpy.ones((len(X),len(X[0])), dtype=object)
    for i in range(len(X)):
        for j in range(len(X[0])):
            if(X[i][j] > Y[i][j]):
                Z[i][j] = X[i][j] + math.log(1+math.exp(Y[i][j]-X[i][j]))
            else:
                Z[i][j] = Y[i][j] + math.log(1+math.exp(X[i][j]-Y[i][j]))
    return Z

#Log Add between Numbers
def logA(a,b):
    if(a>b):
        c = a + math.log(1+math.exp(b-a))
    else:
        c = b + math.log(1+math.exp(a-b))
    return c

def forwardbackward(probability,GenoC, Nc,number, model, t2p, Comodel, LDF1, LDF2, LD1, LD2):
    T2P = t2p
    if(model == 0):
        T1 = numpy.ones((2,2), dtype = object)
        T1 = T1*math.log(0.5)
        T2 = numpy.ones((2,2), dtype = object)
        T2 = T2*math.log(0.5)
        T = numpy.concatenate((T1,T2))
    elif(model == 1):
        Tp = (1-probability)**10
        T1 = numpy.ones((2,2), dtype = object)
        T1[0][0] = math.log(0.5+0.5*Tp)
        T1[0][1] = math.log(0.5-0.5*Tp)
        T1[1][0] = math.log(0.5-0.5*Tp)
        T1[1][1] = math.log(0.5+0.5*Tp)
        T2 = numpy.ones((2,2), dtype = object)
        T2[0][0] = math.log(0.5+0.5*Tp)
        T2[0][1] = math.log(0.5-0.5*Tp)
        T2[1][0] = math.log(0.5-0.5*Tp)
        T2[1][1] = math.log(0.5+0.5*Tp)
        T = numpy.concatenate((T1,T2))
    elif(model == 2):
        T1 = numpy.ones((2,2), dtype = object)
        if(LDF2 == 0):
            T1[0][0] = math.log(LD2)
            T1[0][1] = math.log(1-LD2)
            T1[1][0] = math.log(LD2)
            T1[1][1] = math.log(1-LD2)
        elif(LDF2 == 1):
            T1[0][0] = math.log(float(LD2[0]))
            T1[0][1] = math.log(float(LD2[1]))
            T1[1][0] = math.log(float(LD2[2]))
            T1[1][1] = math.log(float(LD2[3]))

        T2 = numpy.ones((2,2), dtype = object)
        if(LDF1 == 0):
            T2[0][0] = math.log(LD1)
            T2[0][1] = math.log(1-LD1)
            T2[1][0] = math.log(LD1)
            T2[1][1] = math.log(1-LD1)
        elif(LDF1 == 1):
            T2[0][0] = math.log(float(LD1[0]))
            T2[0][1] = math.log(float(LD1[1]))
            T2[1][0] = math.log(float(LD1[2]))
            T2[1][1] = math.log(float(LD1[3]))

        T = numpy.concatenate((T1,T2))

    if(Comodel == 0):
        for i in range(GenoC):
            T3 = numpy.ones((2,2), dtype = object)
            T3[0][0] = math.log(1-probability)
            T3[0][1] = math.log(probability)
            T3[1][0] = math.log(probability)
            T3[1][1] = math.log(1-probability)
            T = numpy.concatenate((T,T3))

    sn = 2**number
    for m in range(0,number,2):
        n = m*2
        j = 2**m
        T1P = logadd(T[n:n+2][0][0]+T2P[:j],T[n:n+2][0][1]+T2P[j:j*2])
        T1P = numpy.concatenate((T1P,logadd(T[n:n+2][1][0]+T2P[:j],T[n:n+2][1][1]+T2P[j:j*2])))
        for i in range(j*2,sn,j*2):
            T1P = numpy.concatenate((T1P,logadd(T[n:n+2][0][0]+T2P[i:i+j],T[n:n+2][0][1]+T2P[i+j:i+j*2])))
            T1P = numpy.concatenate((T1P,logadd(T[n:n+2][1][0]+T2P[i:i+j],T[n:n+2][1][1]+T2P[i+j:i+j*2])))
        if(m == number - 1):
            break
        T2P = logadd(T[n+2:n+4][0][0]+T1P[:j*2],T[n+2:n+4][0][1]+T1P[j*2:j*4])
        T2P = numpy.concatenate((T2P,logadd(T[n+2:n+4][1][0]+T1P[:j*2],T[n+2:n+4][1][1]+T1P[j*2:j*4])))
        for k in range(j*4,sn,j*4):
            T2P = numpy.concatenate((T2P,logadd(T[n+2:n+4][0][0]+T1P[k:k+j*2],T[n+2:n+4][0][1]+T1P[k+j*2:k+j*4])))
            T2P = numpy.concatenate((T2P,logadd(T[n+2:n+4][1][0]+T1P[k:k+j*2],T[n+2:n+4][1][1]+T1P[k+j*2:k+j*4])))

    if(number%2 == 1):
        TP = T1P.T
    else:
        TP = T2P.T
    return TP
