def canil2002(d):
    #Reference: Canil (2002)

    #d is the partition coefficient for V in olv

    param = [0.31, 1.53, 0.9]

    delNNO = (np.log10(param[2]*1/d-1)-param[1])/param[0]

    return  delNNO

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

def putirka2007(H2O, pressure, columns, olv_columns, cat_frac, olv_cat_frac):
    #Reference: Putirka et al. (2007) - eq. 4
    #H2O is in weight percent
    #pressure is in MPa
    #column is list of the headers for cat_frac
    #olv_columns is a list of the headers for olv_cat_frac
    #cat_frac is list of the cation fractions in the melt
    #olv_cat_frac is a list of the cation fractions in the olivine

    #convert pressure from MPa to GPa
    pressure = pressure/1000.0

    #calculate magnesium partitioning
    Dmg = olv_cat_frac[olv_columns.index('MGO')]/cat_frac[columns.index('MGO')]

    #calculate C (L NM) factor #THE LAST PLACE ON THE CAT_FRAC LIST IS FEOT
    Cnm = cat_frac[-1]+cat_frac[columns.index('MNO')]+cat_frac[columns.index('MGO')]+cat_frac[columns.index('CAO')]

    #calculate C (L SiO2) factor
    Csio2 = cat_frac[columns.index('SIO2')]

    #calculate NF factor
    NF = (7./2)*np.log(1.-cat_frac[columns.index('AL2O3')])+7.*np.log(1.-cat_frac[columns.index('TIO2')])

    #solv eqn
    temperature = (15294.6+1318.8*pressure+2.4834*pressure**2.)/(8.048+2.8352*np.log(Dmg)+2.097*np.log(1.5*Cnm)+2.575*np.log(3.0*Csio2)-1.41*NF+0.222*H2O+0.5*pressure)

    return temperature

def sugawara2000(mole_frac, columns, pressure):
    #Reference: Sugawara (2000) - eqn. 6a
    #mole_frac is a list of the molar fractions in the melt
    #columns is a list of the headers for mole_frac
    #Pressure is in MPa

    pressure *= 10.

    temperature = (1446.0-144*(mole_frac[columns.index('SIO2')])-50.0*(mole_frac[-1])+1232.0*(mole_frac[columns.index('MGO')])-389.9*(mole_frac[columns.index('CAO')])+(0.0043*pressure)-540.3*(mole_frac[columns.index('H2O')]))-273.15

    return temperature
