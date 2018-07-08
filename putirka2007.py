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
