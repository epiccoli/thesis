import csv
import pdb
import numpy as np
import time as t
# with open('slice_test.58.csv','rb') as csvfile:
#     spamreader = csv.reader(csvfile, delimiter=',', skipinitialspace=True) # quotechar = 'comment_char'
#     for row in spamreader:
#         pdb.set_trace()
#         for element in row:
#             print(element)
        # print(', '.join(row))
        # pdb.set_trace()

with open('slice_test.58.csv','rb') as f:
    start = t.time()
    n_lines = sum(1 for line in f)
    print t.time() - start
    pdb.set_trace()
    data=f.read()
new_data = data.replace('"','')

for row in csv.reader(new_data.splitlines(),delimiter=',', skipinitialspace=True):
    pdb.set_trace()

    dictionary = {}

    for col, element in enumerate(row):
        try:
            float(element)
        except ValueError:
            print "Not a float"
            dictionary[element]=np.zeros(2)




            
# presidents = ["wash","all"]

# for num, name in enumerate(presidents, start=1):
#     print("President {}: {}".format(num,name))