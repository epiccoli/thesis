import numpy as np

import cantera as ct
import pdb
import matplotlib.pyplot as plt



def list_spec_names(cantera_order,gas):
    species_names=[]
    for i in range(len(cantera_order)):
        if cantera_order[i] == -1:
            species_names.append('T')
        else:
            species_names.append(gas.species_name(cantera_order[i]))


    return species_names




def set_mixture_wagner(gas,Z,phi_j,press=1):

    # DEFINE ZO (oxydiser) side mixture
    Z0compo = get_propane_vitiated_gas()
    T_Z0 = 1500
    C2H4_Z0 = 0.0
    CO2_Z0 = Z0compo['CO2'] # 0.15988194019
    O2_Z0   = Z0compo['O2'] # 0.0289580664922
    N2_Z0   = Z0compo['N2'] # 0.723951662304
    H2O_Z0 = Z0compo['H2O'] # 0.087208331013
    phi_ZO = 0


    # DEFINE Z1 (fuel) side mixture
    Z1compo = getEthyleneJetCompo(phi_j)          # phi_j = 0.8, 1.0, 1.2
    T_Z1 = 300              
    C2H4_Z1 = Z1compo['C2H4']
    CO2_Z1 = 0.0
    O2_Z1   = Z1compo['O2'] 
    N2_Z1   = Z1compo['N2'] 
    H2O_Z1 = 0.0
    phi_Z1 = C2H4_Z1/O2_Z1*96/28

    print phi_Z1, 'phi of the jet'

    # CALCULATE gas state as function of Z (mixture fraction)
    Tempi   = ((T_Z0-T_Z1)/(0.0-1.0))*Z + T_Z0                  # used approximation of Cp = uniform?
    yc2h4   = ((C2H4_Z0-C2H4_Z1)/(0.0-1.0))*Z + C2H4_Z0
    yco2    = ((CO2_Z0-CO2_Z1)/(0.0-1.0))*Z + CO2_Z0
    yo2 = ((O2_Z0-O2_Z1)/(0.0-1.0))*Z + O2_Z0
    yn2 = ((N2_Z0-N2_Z1)/(0.0-1.0))*Z + N2_Z0
    yh2o    = ((H2O_Z0-H2O_Z1)/(0.0-1.0))*Z + H2O_Z0
    phi     = yc2h4/yo2*96/28

    # print "The equivalence ratio is", phi
    # phi is equivalence ratio: stoech ratio massic for methane was 1/4 fuel/oxygen
    # for C3H8 stoech ratio massic is 44/160
    # for C2H4 it is 28/96, corresponding to 1 mole of C2H4 to 3 of O2
    compo="C2H4:"+str(yc2h4)+" O2:"+str(yo2)+" N2:"+str(yn2)+" CO2:"+str(yco2)+" H2O:"+str(yh2o)
    # print "********** Initial state **********"
    print "  - C2H4 mass fraction: "+str(yc2h4)
    print "  - CO2 mass fraction: "+str(yco2)
    print "  - O2 mass fraction : "+str(yo2)
    print "  - N2 mass fraction : "+str(yn2)
    print "  - H2O mass fraction: "+str(yh2o)
    print "  - sum mass fraction: "+str(yc2h4+yo2+yn2+yco2+yh2o), "\n \n"

    print "Temperature of mixture:"
    print(Tempi)

    gas.TPY = Tempi, press*1.01325e5, compo

    return gas, phi

    
def getEthyleneJetCompo(phi):

    # Calculation of fresh gases composition
    # phi = 0.8 - 1.2 equivalence ratio

    # C2H4 combustion reaction:
    # C2H4 + 3 O2 => 2 CO2 + 2 H2O

    # molar mixing ratio of dry air:
    # oxygen: 0.21
    # nytrogen: 0.78
    # rest is neglected

    # everything done in moles, then converted to mass at the end using molar masses


    # Molar masses
    mm_C2H4 = 2*12.0 + 4*1.0
    mm_O2 = 2*16.0
    mm_CO2 = 12.0 + 2*16.0
    mm_H20 = 2*1.0 + 16.0
    mm_N2 = 2*14.0

    Fuel2Oxygen_mole=1/3.0      # stoechiometric combustion

    # moles of fresh gases
    n_C2H4 = 1          # keep moles of fuel as reference = 1 
    n_O2 = 1/Fuel2Oxygen_mole/phi
    n_N2 = 0.78/0.21*n_O2


    # masses of cross-jet components
    m_C2H4 = n_C2H4*mm_C2H4
    m_O2 = n_O2*mm_O2
    m_N2 = n_N2*mm_N2

    # total number of moles and total mass
    n_tot = n_N2 + n_O2 + n_C2H4
    m_tot = m_N2 + m_O2 + m_C2H4

    # mass fractions of vitiated atmosphere components
    f_O2_mass = m_O2/m_tot
    f_N2_mass = m_N2/m_tot
    f_C2H4_mass = m_C2H4/m_tot

    check_sum = f_N2_mass + f_O2_mass + f_C2H4_mass

    # print f_O2_mass, " O2"
    # print f_N2_mass, " N2"
    # print f_C2H4_mass, " C2H4"



    # print "Sum of mass fractions: "
    # print check_sum


    compo = {'O2':f_O2_mass, 'N2':f_N2_mass, 'C2H4':f_C2H4_mass}

    return compo


def get_propane_vitiated_gas():

    # Calculation of first combustion vitiated atmosphere
    # phi = 0.87 equivalence ratio

    # C3H8 combustion:
    # C3H8 + 5 O2 => 3 CO2 + 4 H2O

    # molar mixing ratio of dry air:
    # oxygen: 0.21
    # nytrogen: 0.78
    # rest is neglected --> 1 O2 + 3.76 N2

    # everything done in moles, then converted to mass at the end using molar masses

    # Molar masses
    mm_C3H8 = 3*12.0 + 8*1.0
    mm_O2 = 2*16.0
    mm_CO2 = 12.0 + 2*16.0
    mm_H20 = 2*1.0 + 16.0
    mm_N2 = 2*14.0

    phi = 0.87 
    FuelOxygen_mole=1/5.0       # stoechiometric combustion

    # moles of fresh gases
    n_C3H8_fresh = 1            # keep moles of fuel as reference = 1 
    n_O2_fresh = 5/phi
    n_N2_fresh = 3.76*n_O2_fresh

    # moles of vitiated atmosphere components
    n_N2_vitiated = n_N2_fresh
    n_O2_left = n_O2_fresh - n_C3H8_fresh/FuelOxygen_mole
    n_CO2 = 3
    n_H2O = 4
    # masses of vititated atmosphere components
    m_N2 = n_N2_vitiated*mm_N2
    m_O2_left = n_O2_left*mm_O2
    m_CO2 = n_CO2*mm_CO2
    m_H20 = n_H2O*mm_H20

    # total number of moles and total mass
    n_tot = n_N2_vitiated + n_O2_left + n_CO2 + n_H2O
    m_tot = m_N2 + m_O2_left + m_CO2 + m_H20

    # mass fractions of vitiated atmosphere components
    f_O2_mass = m_O2_left/m_tot
    f_N2_mass = m_N2/m_tot
    f_CO2_mass = m_CO2/m_tot
    f_H2O_mass = m_H20/m_tot

    check_sum = f_O2_mass + f_N2_mass + f_CO2_mass + f_H2O_mass

    print f_O2_mass, " O2"
    print f_N2_mass, " N2"
    print f_CO2_mass, " CO2"
    print f_H2O_mass, "H2O"

    print "Sum of mass fractions: "
    print check_sum

    print "Mass fraction of unreduced oxygen in vitiated atmosphere (% total)"
    print f_O2_mass*100
    print "Corresponds to 3\% value given"
    # compo = 'O2:' + str(f_O2_mass) + ' N2:' + str(f_N2_mass) + ' CO2:' + str(f_CO2_mass) + ' H2O:' + str(f_H2O_mass)
    compo = {'O2':f_O2_mass, 'N2':f_N2_mass, 'CO2':f_CO2_mass, 'H2O':f_H2O_mass}
    return compo




