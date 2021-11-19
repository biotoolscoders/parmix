import numpy

def emI2(EC,ls,ab,A, B, Nc, dr, Ge):
    emission = numpy.ones((ls), dtype=object)

    for i in range(ls):
        for x in range(Nc):
            if(dr[x] == '0'):
                if(EC[i][ab][0] == '0'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * (A * (1-Ge) + (1-A) * Ge)
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * ((1-A) *(1-Ge) + A * Ge)  
                elif(EC[i][ab][0] == '1'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * (B * (1-Ge) + (1-B) * Ge)
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * ((1-B) *(1-Ge) + B * Ge)

            elif(dr[x] == '1'):
                if(EC[i][ab][1] == '0'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * (A * (1-Ge) + (1-A) * Ge)
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * ((1-A) *(1-Ge) + A * Ge)  
                elif(EC[i][ab][1] == '1'):
                    if(EC[i][ab][2+x] == '0'):
                        emission[i] = emission[i] * (B * (1-Ge) + (1-B) * Ge)
                    elif(EC[i][ab][2+x] == '1'):
                        emission[i] = emission[i] * ((1-B) *(1-Ge) + B * Ge) 
    return emission
