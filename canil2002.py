def canil2002(d):
    #Reference: Canil (2002)

    #d is the partition coefficient for V in olv

    param = [0.31, 1.53, 0.9]

    delNNO = (np.log10(param[2]*1/d-1)-param[1])/param[0]

    return  delNNO
