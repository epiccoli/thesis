import numpy as np
import scipy.linalg as LA

A = np.matrix('1,2,3; 0,0,0; 1,2,3')
# A = np.matrix([[1,2,3], [0,0,0], [1,2,3]])

# Verify construction of matrix in python, both methods equivalent
#print A

D, l, r = LA.eig(A, left = True)
D = np.diag(D)

# Verify getting the same results as in matlab 
# print D

# print l

# print r

# outcome: correctly reproduced matlab command [V,D,l] = eig(A)



# TRY TO VERIFY RELATIONS:
#  l'*A = D*l' left eigenvector l
# A*V = V*D right eigenvector V

# print np.subtract(np.matmul(l.conj().T, A),np.matmul(D,l.conj().T)) 
# print np.subtract(np.matmul(A,r),np.matmul(r,D)) # 

# outcome: verified relations, np.matmul for matrix multiplication.

i = 0
a = r[:,i]
b = l[i,:]

aibi=np.zeros(len(a))
for j in range(len(a)):

	aibi[j] = abs(a[j]*b[j])

S = np.sum(aibi)

aibi = aibi/S

print aibi