import pdb
import numpy as numpy

import cantera as ct 
import sys

sys.path.append('../../ethylene/')       # path to pyjacob.so
sys.path.append('../')               # path to read_slice.py
# from read_slice import *
# from CEMA import *
import pyjacob
import matplotlib.pyplot as plt
import os 
import pandas as pd
import scipy.spatial
import numpy as np
import numpy.ma as ma 
import time as t
import scipy.linalg as LA

def parallelism(v1,v2):
    # returns [-1,1]
    num = np.dot(v1,v2)
    den = LA.norm(v1)*LA.norm(v2)
    # high value of parallelism 
    parallelism = num/den

    return parallelism

def EI(D,l,r,k):

    # D is 1D array of eigenvalues, corresponding to l(eft) eigenvector and r(ight) eigenvector matrices
    # returns the E(xplosive) I(ndex) calculated as 

    
    a = r[:,k]
    b = l[:,k] # changed this according to documentation 
    
    jac_dim = len(a)

    ab=np.zeros(jac_dim)
    for j in range(jac_dim):

        ab[j] = abs(a[j]*b[j])

    S = np.sum(ab)

    ab = ab/S

    return ab


def get_state_vec_keys(gas):

    n_species = gas.n_species

    # phi is the vector of variable names to look for in the columns first row (ordered for pyjac)
    phi = []
    phi.append('temperature')
    for i in range(gas.n_species):
        if gas.species_name(i) != 'N2':
            phi.append(gas.species_name(i))
    return phi

def solve_jacobian(P,y):

    n_species = len(y)
    dydt = np.zeros_like(y)
    pyjacob.py_dydt(0, P, y, dydt)
    
    #create a jacobian vector
    jac = np.zeros(n_species*n_species)

    #evaluate the Jacobian
    pyjacob.py_eval_jacobian(0, P, y, jac)
    jac = jac.reshape(n_species,n_species)

    D, L, R = LA.eig(jac, left = True)
    
    D = D.real 
    L = L.real
    R = R.real

    return D, L, R

def build_state_df(df_row,state_vec_keys):
    # df_row is one row of dataframe, state_vec_keys contains the headers of interest for the state vector

    n_species = len(state_vec_keys)
    # fill the state vector y with values from 
    y = np.zeros(n_species)
    for i in range(n_species):
        try:
            y[i] = df_row[state_vec_keys[i]]

            if y[i] < 0:
                y[i] = 0
        
        except KeyError:
            print('Data column {} not found in csv'.format(state_keys[i]))   
            pass
        
    return y

def eig_line(df_row,state_vec_keys):

    y = build_state_df(df_row,state_vec_keys)
    P = df_row['pressure']

    # Solve eigenvalue problem 
    D, L, R = solve_jacobian(P,y)

    return D, L, R 

def do_kdtree(coordinates,query_point,nb_neigh):

    # [IN]:     coordinates is an array containing n rows, m columns, n being the number of coordinates
    #               m is the dimensionality of the space
    #           query_point is a 1 d array with the m coordinates of the point
    # 
    # [OUT]:    returns the indices of the closest nb_neigh points
    
    mytree = scipy.spatial.cKDTree(coordinates)
    dist, ind = mytree.query(query_point, k=nb_neigh)
    return ind

def create_tree(coordinates):

    mytree = scipy.spatial.cKDTree(coordinates)
    return mytree

def query_tree(mytree,query_point,nb_neigh):

    dist, ind = mytree.query(query_point, k=nb_neigh)
    return ind

def eig_line(df_row,state_vec_keys):

    y = build_state_df(df_row,state_vec_keys)
    P = df_row['pressure']

    # Solve eigenvalue problem 
    D, L, R = solve_jacobian(P,y)

    return D, L, R 


########### REAL START OF SCRIPT ############

# processed_csv = './processed-csv/panda_processed_J87_phi12_snapshot_withQSS.csv'
processed_csv = './processed-csv/panda_subdomain_J87_phi12_snapshot_withQSS.csv'

# READ FILE IN DATAFRAME

# read contents of file into dataframe
whole_df = pd.read_csv(processed_csv)
all_coord = whole_df[['x','y','z']].values

# Create tree of coordinates
start = t.time()
coord_tree = create_tree(all_coord)
print('Time to create tree: {:.4f} sec'.format(t.time()-start))

# select starting row with max eigenvalue in subdomain
start_row = whole_df.loc[whole_df['max_eig'].idxmax()]      # is row dataframe
start_row_idx = whole_df['max_eig'].argmax()
whole_df.set_value(start_row_idx,'CEM',start_row['max_eig']) # set CEM to maximum value (by def)

# COMPUTE EI
###### create state vector ######
gas = ct.Solution('../Skeletal29_N.cti')    
n_species = gas.n_species           
# create list with all names (str) of components of the state vector for the chemical jacobian
state_vec_keys = get_state_vec_keys(gas)

###### Solve eigenvalue problem at node i ######
D, L, R = eig_line(start_row,state_vec_keys)
max_eig_idx = np.argmax(D)

###### Compute EI of start point ######
ei_prev = EI(D,L,R,max_eig_idx)
coord_prev = start_row[['x','y','z']].values

############# DEBUG FOLLOW #############
Dj=9.53e-3
processed_idx=[]
############# DEBUG FOLLOW ############# 

processed = 0
stop_LOOP = 0

big_start = t.time()

while pd.isnull(whole_df['CEM']).any() == True:

    start_LOOP = t.time()        # time 
    ###### Follow CEM in closest neigh with EI ######
    # Find indices of nb_neigh closer points (neighbouring points)
    nb_neigh = 1000
    
    
    neigh_idx = query_tree(coord_tree, coord_prev, nb_neigh)
    # print('Takes {:.5f} seconds with {:d} neighbours'.format(t.time()-start,nb_neigh))    # time

    
    # Find index of closest non processed point (neigh_idx[count])
    count = 0
    
    while pd.isnull(whole_df.iloc[neigh_idx[count]]['CEM']) == False:
        count +=1

    D, L, R = eig_line(whole_df.iloc[neigh_idx[count]], state_vec_keys)
        
    # Process and append real CEM value
    alignment = np.zeros(n_species)

    for pos in xrange(n_species):
        ei_test = EI(D,L,R,pos)             # spare time
        alignment[pos] = parallelism(ei_test,ei_prev)

    best_fit_idx = np.argmax(alignment)
    # print('Sec eval align : {:.4f}'.format(t.time()-start))

    # Write best aligned eig
    CEM = D[best_fit_idx]
    whole_df.set_value(neigh_idx[count], 'CEM', CEM)


    
    ############# DEBUG FOLLOW #############
    # pdb.set_trace()
    if processed > 1000:
        break
    fig, ax = plt.subplots(figsize=(10,10))

    if processed == 0:
        ax.plot(all_coord[neigh_idx[count]][0]/Dj,all_coord[neigh_idx[count]][1]/Dj,'bo',markerfacecolor='white',label='next')
        ax.plot(coord_prev[0]/Dj,coord_prev[1]/Dj,'ro',label='curr')
        processed_idx.append(neigh_idx[count])
    else:
        ax.plot(np.array(all_coord[processed_idx])[:,0]/Dj,np.array(all_coord[processed_idx])[:,1]/Dj,'ko',markerfacecolor='white',label='past')
        ax.plot(coord_prev[0]/Dj,coord_prev[1]/Dj,'ro',label='curr')
        ax.plot(all_coord[neigh_idx[count]][0]/Dj,all_coord[neigh_idx[count]][1]/Dj,'bo',markerfacecolor='white',label='next')
        processed_idx.append(neigh_idx[count])
    plt.legend()
    plt.savefig('ei_debug{:d}.pdf'.format(processed))
    
    ############# DEBUG FOLLOW #############

    ## END TASKS OF LOOP ## 
    
    ei_prev = EI(D,L,R,best_fit_idx)
    coord_prev = whole_df.iloc[neigh_idx[count]][['x','y','z']].values

    processed += 1 



    if processed > whole_df.shape[0]:
        break

    stop_LOOP += t.time() - start_LOOP      # time
    if processed % 1000 == 0:
      print('Process time for 1000 points: {:.3f} sec'.format(stop_LOOP))
      stop_LOOP = 0


    # pdb.set_trace()

if pd.isnull(whole_df['CEM']).any() == True:
    print "unfinished"

else:
    print "finished"
    print('Total processing time (Ei follow): {:.4f}'.format(t.time()-big_start))

whole_df.to_csv('./processed-csv/ei_followed_J87_phi12_snapshot_withQSS.csv')



        # test how many alignment are below threshold 1e-2 : take a representative value #savetime
        




