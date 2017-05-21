import csv
import pdb
import numpy as np
import time as t
import cantera as ct
import scipy.linalg as LA

def load_row(file_in,row_nb):
    # returns only line number "row_nb" of "file_in"
    with open(file_in,'rb') as f:   # suggested to open in binary mode
        reader = csv.reader(f)
        
        for i, split_line in enumerate(reader):
            if i == row_nb:
                return split_line

def load_val(file_in,row_nb,col_nb):

    # row_nb corresponds to a certain point in the domain, 
    # col_nb corresponds to the index of the columnt of a certain variable

    split_line = load_row(file_in,row_nb)
    
    try:
        value = float(split_line[col_nb])
    except ValueError:
        print "Error: trying to extract non numeric value"
    
    return value

# def make_state_vector(file_in)

def count_lines(file_in):
    # returns number of lines in file_in
    num_lines = sum(1 for line in open(file_in,'rb'))   # suggested to open in binary mode

    return num_lines

def get_variable_dict(file_in):
    # Load header of csv file (contains variable names stored in the file)
    var_names = load_row(file_in,0)
    # create dictionary with variable names as keys, values are indices of columns
    variable_dict = dict(zip(var_names,range(0,len(var_names))))
    return variable_dict

def build_state_vector(file_in,gas,row_nb):
    # y is state vector for jacobian calculation
    # phi are the names identifying the columns in the csv data
    # row_nb is the row of interest 

    variable = get_variable_dict(file_in)
    
    n_species = gas.n_species

    # phi is the vector of variable names to look for in the columns first row (ordered for pyjac)
    phi = []
    phi.append('temperature')
    for i in range(gas.n_species):
        if gas.species_name(i) != 'N2':
            phi.append(gas.species_name(i))
    
    # fill the state vector y with values from 
    y = np.zeros(n_species)
    for i in range(n_species):
        try:
            y[i] = load_val(file_in, row_nb, variable[phi[i]])
        except KeyError:
            # print('still need to extract those values of QSS')   
            pass

    return y

def get_coordinates(file_in, row_nb):

    column_name = get_variable_dict(file_in)

    # Check if actually not mixed up order
    x = load_val(file_in,row_nb,column_name['Points:0'])
    y = load_val(file_in,row_nb,column_name['Points:1'])
    z = load_val(file_in,row_nb,column_name['Points:2'])

    return np.array([x, y, z])


def csv_append(line, path):

    with open(path, 'ab') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        
        writer.writerow(line)

            # line_write.tofile('foo.csv',sep=',',format='%1.5e')


# data = np.array([1.001012001230,2,3.0,4])
# stromg = 'Flat,second,thri'.split(',')



