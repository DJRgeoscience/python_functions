def fecalc(d, columns, FO2, T):
    #Reference: Kress and Carmichael (1991)
    #d - input data consisting of major element concentrations in weight percent
    #columns- column headers for d (must be in the mm dictionary below)
    #FO2 - oxygen fugacity (e.g., fo2buffer output)
    #T is in C

    #Molar masses
    mm = {'SIO2': 60.08, 'TIO2': 79.866, 'AL2O3': 101.96, 'FE2O3': 159.69, 'FEO': 71.844, 'MNO': 70.9374, 'MGO': 40.3044, 'CAO': 56.0774, 'NA2O': 61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528, 'CO2': 44.01, 'S': 32.065, 'CL': 35.453, 'F': 18.9984}

    #param - Kress & Carmichael parameters
    param = [0.196, 11492.0, -6.675, -2.243, -1.828, 3.201, 5.854, 6.215]
    #param = [a, b, c, dAl2O3, dFeO*, dCaO, dNa2O, dK2O]

    mole_pro = []
    mole_headers = []
    count = 0
    for i in d:
        if columns[count] in mm:
            mole_headers.append(columns[count])
            mole_pro.append(i/mm[columns[count]])
        elif columns[count] == 'FEOT':
            mole_headers.append(columns[count][:-1])
            mole_pro.append(i/mm[columns[count][:-1]])
        count += 1

    mole_frac = [x/sum(mole_pro) for x in mole_pro]
    T += 273.15

    Fe2O3FeO = np.exp((param[0]*np.log(FO2))+(param[1]/T)+param[2]+((param[3]*mole_frac[mole_headers.index('AL2O3')])+(param[4]*mole_frac[mole_headers.index('FEO')])+(param[5]*mole_frac[mole_headers.index('CAO')])+(param[6]*mole_frac[mole_headers.index('NA2O')])+(param[7]*mole_frac[mole_headers.index('K2O')])))

    Fe3FeT = 2.0/(2.0+(1.0/Fe2O3FeO))

    return Fe3FeT
