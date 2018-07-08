def sugawara2000(mole_frac, columns, pressure):
    #Reference: Sugawara (2000) - eqn. 6a
    #mole_frac is a list of the molar fractions in the melt
    #columns is a list of the headers for mole_frac
    #Pressure is in MPa

    pressure *= 10.

    temperature = (1446.0-144*(mole_frac[columns.index('SIO2')])-50.0*(mole_frac[-1])+1232.0*(mole_frac[columns.index('MGO')])-389.9*(mole_frac[columns.index('CAO')])+(0.0043*pressure)-540.3*(mole_frac[columns.index('H2O')]))-273.15

    return temperature
