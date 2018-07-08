def fo2buffer(T, P, delta, buff):
    #Reference: Reviews in Mineralogy Volume 25
    #T is in C
    #P is in MPa
    #delta is the delta value from NNO or QFM
    #buff - text indicating NNO/FMQ/QFM

    T += 273.15
    P *= 10.0

    if buff in ['FMQ','QFM']:
        FO2 = 10**((-25096.3/T)+8.735+(0.110*(P-1)/T)+delta)
    else:
        FO2 = 10**((-24930/T)+9.36+(0.046*(P-1)/T)+delta)

    return FO2
