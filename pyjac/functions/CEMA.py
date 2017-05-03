import numpy as np
import pyjacob
import scipy.linalg as LA
import cantera as ct
import pdb
import matplotlib.pyplot as plt

def solve_eig_gas(gas):


    # input arg: gas cantera object
    # output arg: eigenvalue vector, left and right eigenvector matrices 

    T = gas.T
    P = gas.P
    
    # #setup the state vector
    y = np.zeros(gas.n_species)
    y[0] = T
    y[1:] = gas.Y[:-1]
    
    # #create a dydt vector
    dydt = np.zeros_like(y)
    pyjacob.py_dydt(0, P, y, dydt)

    #create a jacobian vector
    jac = np.zeros(gas.n_species * gas.n_species)


    #evaluate the Jacobian
    pyjacob.py_eval_jacobian(0, P, y, jac)

    jac = jac.reshape(gas.n_species,gas.n_species)

    # Solve eigenvalue PB > D: eigenvalues
    D, vl, vr = LA.eig(jac, left = True)
    D=D.real

    return D, vl, vr


def EI(D,l,r,k):

    # D is 1D array of eigenvalues, corresponding to l(eft) eigenvector and r(ight) eigenvector matrices
    # returns the E(xplosive) I(ndex) calculated as 

    
    a = r[:,k]
    b = l[k,:]  
    
    jac_dim = len(a)

    ab=np.zeros(jac_dim)
    for j in range(jac_dim):

        ab[j] = abs(a[j]*b[j])

    S = np.sum(ab)

    ab = ab/S

    return ab

# def PI()


def find_max_eigenvalue(gas):

    # input arg: gas cantera object
    # output arg: max chemical eigenvalue

    T = gas.T
    P = gas.P
    
    # #setup the state vector
    y = np.zeros(gas.n_species)
    y[0] = T
    y[1:] = gas.Y[:-1]
    
    # #create a dydt vector
    dydt = np.zeros_like(y)
    pyjacob.py_dydt(0, P, y, dydt)

    #create a jacobian vector
    jac = np.zeros(gas.n_species * gas.n_species)


    #evaluate the Jacobian
    pyjacob.py_eval_jacobian(0, P, y, jac)

    jac = jac.reshape(gas.n_species,gas.n_species)

    # Solve eigenvalue PB > D: eigenvalues
    D, W, V = LA.eig(jac, left = True)

    # FINDS HIGHEST EIGENVALUE FOR EACH POINT
    # Take only real part into account
    D = D.real
    max_val = 0
    for i in range(len(D)):
        if D[i] > max_val:
            # print D[i], gas.species_name(i)
            max_val = D[i]
            species_idx = i
            print species_idx

    return max_val


def solve_eig_flame(f,gas):
    N_eig = 3 # selects number of eigenvalues to store for each flame location
    N_EI = 1 # number of EI species to track 

    T = f.T # 1D array with temperatures
    Y = f.Y # matrix array with lines corresponding to the 29 species, columns corresponding to the grid points
    P = f.P # single value

    n_species = gas.n_species 
    grid_pts = len(f.grid)

    eigenvalues = np.zeros([N_eig, grid_pts])
    # manual_eig = np.zeros(grid_pts)
    eig_position = np.zeros(grid_pts)
    
    track_specs=[]      # initialised as list, then converts to np.array when using np.union1d
    global_expl_indices = np.zeros([n_species, grid_pts])
    
    # iterate over x-span of 1D flame domain
    for loc in range(grid_pts):
        y=np.zeros(n_species)
        y[0] = T[loc]
        # find the position of N2 species
        N2_idx = gas.species_index('N2')   
        y_massfr = np.concatenate([Y[0:N2_idx,loc], Y[N2_idx+1:,loc]])
        y[1:] = y_massfr

        jac = create_jacobian(T,P,y)
        
        D, vl, vr = LA.eig(jac, left = True)
        D = D.real
        # print eigenvalues
        # if loc % 10 == 0:
        #     print D
        #     pdb.set_trace()
        # Store the N most positive eigenvalues
        eigenvalues_loc = D[np.argsort(D)[-N_eig:]]

        max_idx = np.argmax(D)
        max_eig = D[max_idx]    
        
        # Store explosive indices corresponding to chemical explosive mode at x loc(ation) in 1D domain
        expl_indices = EI(D,vl,vr,max_idx)
        # Find indices that would sort the EI (from lowest to highest)
        sorted_idx = np.argsort(expl_indices) 

        # select indices of 5 most important species (attention: pyjac position index!! need converting with reorder_species)
        main_species_local = sorted_idx[-N_EI:]
        # Union operation between indices of main species at x loc(ation) and previous important species
        track_specs = np.union1d(main_species_local,track_specs)
        
        global_expl_indices[:,loc] = expl_indices 
        eigenvalues[:,loc] = eigenvalues_loc
        
        # if loc == 0:
        #     idx_max_eig = np.argmax(D)

        # manual_eig[loc] = D[idx_max_eig]
        eig_position[loc] = np.argmax(D)

    plt.figure()
    plt.plot(f.grid,eig_position)
    plt.show()
        
    
    track_specs = track_specs.astype(np.int64)

    return eigenvalues, global_expl_indices, track_specs

def reorder_species(pyjac_indices, N2_idx):

    # subtract one position to pyjac indices up to index where N2 was taken out of species list:
    # example: cantera species [O H N2 CH3 CH4], pyjac y vector [T O H CH3 CH4]
    # --> to bring back to original cantera species indices, need to subtract 1 to positions up to index of N2 (2 here)
    
    cantera_idx = pyjac_indices
    
    for i in range(len(cantera_idx)):

        if cantera_idx[i] == 0:
            # in the case that temperature has one of the highest explosive indices
            print "Temperature is EI --> add treatment of this case"
            pdb.set_trace()

        elif cantera_idx[i] <= N2_idx:
            cantera_idx[i] -= 1  

    return cantera_idx

def list_spec_names(cantera_order,gas):
    species_names=[]
    for i in range(len(cantera_order)):

        species_names.append(gas.species_name(cantera_order[i]))


    return species_names

def get_species_names(tracked_species_idx, gas):

    # Returns dictionary with trackes species names as keys, indices in pyjac notation as values
    # tracked_species_idx is in pyjac notation (reordered N2 at the end).
    N2_idx = gas.species_index('N2')   
    # Revert pyjac ordering
    cantera_species_idx = reorder_species(tracked_species_idx, N2_idx)
    # get species names corresponding to cantera_species_idx 
    species_names = list_spec_names(cantera_species_idx, gas)
    # create dictionary with species names (string) as keys and pyjac indices as values
    dictionary = dict(zip(species_names,tracked_species_idx))

    return dictionary


def create_jacobian(T,P,y):

    n_species = len(y)
    dydt = np.zeros_like(y)
    pyjacob.py_dydt(0, P, y, dydt)
    
    #create a jacobian vector
    jac = np.zeros(n_species*n_species)

    #evaluate the Jacobian
    pyjacob.py_eval_jacobian(0, P, y, jac)
    jac = jac.reshape(n_species,n_species)

    return jac


def setSolutionProperties(gas,Z,press=1):

    # DEFINE ZO (oxydiser) side mixture
    T_Z0 = 1500
    C2H4_Z0 = 0.0
    CO2_Z0 = 0.15988194019
    O2_Z0   = 0.0289580664922
    N2_Z0   = 0.723951662304
    H2O_Z0 = 0.087208331013
    phi_ZO = 0
    # DEFINE Z1 (fuel) side mixture

    Z1compo = getEthyleneJetCompo(.8)          # phi_j = 0.8, 1.0, 1.2
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

    # print "Temperature of mixture:"
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