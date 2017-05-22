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
import matplotlib.pyplot as plt


# with open('slice_test.58.csv','rb') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=',', skipinitialspace=True) # quotechar = 'comment_char'
#     for row in spamreader:
#         pdb.set_trace()
#         for element in row:
#             print(element)
        # print(', '.join(row))

def main():
    file_in = 'slice_test.58.csv'                       # argument 1
    file_in = 'phi1dot2_Tub320K_sL079ms_Ldomain0dot03m.csv'
    test_out = '0test_EIG_write.csv'
    test_out = 'new_debug_flame.csv'

    column_name = get_variable_dict(file_in)
    n_lines = count_lines(file_in)
    print n_lines
    ###### create state vector ######
    gas = ct.Solution('Skeletal29_N.cti')               # argument 2

    header_line = 'Points:0,Points:1,Points:2,max_eig'.split(',')
    csv_append(header_line,test_out)
    pdb.set_trace()
    for point in range(1,n_lines):#count_lines(file_in)):
        big_start = t.time()
        # Load coordinates of node i
        coord = get_coordinates(file_in,point)
        # Create state vector for node i
        y = build_state_vector(file_in,gas,point)
        # Load pressure value at node i
        P = load_val(file_in,point,column_name['pressure'])
        # Solve eigenvalue problem at node i
        D, L, R = solve_jacobian(P,y)
        max_eig = np.amax(D)

        line_write = np.append(coord, max_eig)
        csv_append(line_write,test_out)

        print t.time() - big_start, 'one loop'
    
    column_name = get_variable_dict(test_out)
    pdb.set_trace()

    n_lines = count_lines(test_out)
    x_coord = np.zeros(n_lines-1)
    CEM = np.zeros(n_lines-1)


    
    for point in range(1,count_lines(test_out)-1):
        coord = get_coordinates(test_out,point)
        x_coord[point] = coord[0]
        CEM[point] = load_val(test_out,point,column_name['max_eig'])


    plt.figure()
    plt.plot(x_coord,CEM/1e6,'.')
    plt.show()


    # point = 1 # point in space (idx)
    # # Load coordinates of node
    # coord = get_coordinates(file_in,point)
    # # Create state vector for node
    # y = build_state_vector(file_in,gas,point)
    # # Load pressure at node
    # P = load_val(file_in,point,column_name['pressure'])

    
    # # Solve eigenvalue problem at node 
    # D, L, R = solve_jacobian(P,y)

    # max_eig = np.amax(D)
    # pdb.set_trace()
    # line_write = np.append(coord, max_eig)
    # print line_write


    # csv_append(line_write,'0test_EIG_write.csv')
    # pdb.set_trace()


    # for i in range(1,count_lines(file_in)):
    #     print i
    #     pdb.set_trace()

    # print load_val(file_in,1,variable['zeta_y'])



if __name__ == '__main__':
    main()
pdb.set_trace()

# with open(file_in,'rb') as f:

#     for element in csv.reader(first_line,delimiter=','):
#         print element


#     pdb.set_trace()
#     start = t.time()
#     data=f.read()
#     print t.time() - start

    
# new_data = data.replace('"','')

# header=[]
# for i_line,row in enumerate(csv.reader(new_data.splitlines(),delimiter=',')):
#     # pdb.set_trace()
#     for element in row:
#         try:
#             float(element)

#         except ValueError:
#             if i_line == 0:
#                 header.append(element)
#             elif i_line > 0:
#                 print "Non float encountered outside of header row -> check format of csv file"



            

    # print dictionary


            
# presidents = ["wash","all"]

# for num, name in enumerate(presidents, start=1):
#     print("President {}: {}".format(num,name))