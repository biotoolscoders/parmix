import numpy

def emI1(EC,Ge,ls,la,Nc):
    emission = numpy.ones((ls, la), dtype=object)

    for i in range(ls):
        for j in range(la):
            for x in range(Nc):
                if(EC[i][j][x] == '0'):
                    if(EC[i][j][Nc] == '0'):
                        if(EC[i][j][Nc+2+x] == '0'):
                            emission[i][j] = emission[i][j] * (1-Ge)
                        elif(EC[i][j][Nc+2+x] == '1'):
                            emission[i][j] = emission[i][j] * (Ge)
                    elif(EC[i][j][Nc] == '1'):
                        if(EC[i][j][Nc+2+x] == '0'):
                            emission[i][j] = emission[i][j] * (Ge)
                        elif(EC[i][j][Nc+2+x] == '1'):
                            emission[i][j] = emission[i][j] * (1-Ge)      

                elif(EC[i][j][x] == '1'):
                    if(EC[i][j][Nc+1] == '0'):
                        if(EC[i][j][Nc+2+x] == '0'):
                            emission[i][j] = emission[i][j] * (1-Ge)
                        elif(EC[i][j][Nc+2+x] == '1'):
                            emission[i][j] = emission[i][j] * (Ge)
                    elif(EC[i][j][Nc+1] == '1'):
                        if(EC[i][j][Nc+2+x] == '0'):
                            emission[i][j] = emission[i][j] * (Ge)
                        elif(EC[i][j][Nc+2+x] == '1'):
                            emission[i][j] = emission[i][j] * (1-Ge) 
    return emission

