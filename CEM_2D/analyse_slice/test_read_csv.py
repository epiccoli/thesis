import csv
import pdb
import numpy as np
import time as t
import cantera as ct
import scipy.linalg as LA
import sys
sys.path.append('../ethylene/')       # path to pyjacob.so
from CEMA import *
from read_slice import *
import pyjacob


# with open('slice_test.58.csv','rb') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=',', skipinitialspace=True) # quotechar = 'comment_char'
#     for row in spamreader:
#         pdb.set_trace()
#         for element in row:
#             print(element)
        # print(', '.join(row))


def csv_append(line, path):

    with open(path, 'ab') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        
        writer.writerow(data)

data = np.array([1.001012001230,2,3.0,4])
stromg = 'Flat,second,thri'.split(',')

def get_coordinates(file_in, row_nb):

    column_name = get_variable_dict(file_in)

    # Check if actually not mixed up order
    x = load_val(file_in,row_nb,column_name['Points:0'])
    y = load_val(file_in,row_nb,column_name['Points:1'])
    z = load_val(file_in,row_nb,column_name['Points:2'])

    return np.array([x, y, z])



def main():
    file_in = 'slice_test.58.csv'

    column_name = get_variable_dict(file_in)

    ###### create state vector ######
    gas = ct.Solution('Skeletal29_N.cti')

    point = 1 # point in space (idx)
    # Load coordinates of node
    coord = get_coordinates(file_in,point)
    # Create state vector for node
    y = build_state_vector(file_in,gas,point)
    # Load pressure at node
    P = load_val(file_in,point,column_name['pressure'])

    
    # Solve eigenvalue problem at node 
    D, L, R = solve_jacobian(P,y)

    max_eig = np.amax(D)
    pdb.set_trace()
    line_write = np.append(coord, max_eig)
    print line_write

    csv_append(line_write,'test_EIG_write.csv')
    pdb.set_trace()


    # for i in range(1,count_lines(file_in)):
    #     print i
    #     pdb.set_trace()

    print load_val(file_in,1,variable['zeta_y'])



if __name__ == '__main__':
    main()
pdb.set_trace()

with open(file_in,'rb') as f:

    for element in csv.reader(first_line,delimiter=','):
        print element


    pdb.set_trace()
    start = t.time()
    data=f.read()
    print t.time() - start

    
new_data = data.replace('"','')

header=[]
for i_line,row in enumerate(csv.reader(new_data.splitlines(),delimiter=',')):
    # pdb.set_trace()
    for element in row:
        try:
            float(element)

        except ValueError:
            if i_line == 0:
                header.append(element)
            elif i_line > 0:
                print "Non float encountered outside of header row -> check format of csv file"



            

    # print dictionary


            
# presidents = ["wash","all"]

# for num, name in enumerate(presidents, start=1):
#     print("President {}: {}".format(num,name))