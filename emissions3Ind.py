import numpy

def emI3(EC,ls,ab,A,B,Nc,dr,Ge):
    emission = numpy.ones((ls), dtype=object) 

    for i in range(ls):
        for x in range(Nc):
            if(dr[x] == '0'):
                if(EC[i][ab][0] == '0'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * (1-Ge)
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * Ge
                elif(EC[i][ab][0] == '1'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * Ge
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * (1-Ge)

            elif(dr[x] == '1'):
                if(EC[i][ab][1] == '0'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * (1-Ge)
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * Ge
                elif(EC[i][ab][1] == '1'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * Ge
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * (1-Ge)
    return emission

