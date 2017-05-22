import numpy as np
import pdb
import os
import sys
sys.path.append('../../functions/')
from CEMA import EI

def eigenmode(D,M):

	# D is vector of real parts of eigenvalues

	# Identify the index of the maximum eigenvalue
	max_idx = np.argmax(D)
	# Find the corresponding eigenvector matrix
	mode = M[:,max_idx].real

	return mode

def test_ei_parallel(L,R,ei_vec0):
	scores = np.zeros(len(L))

	for k in range(len(L)):

		a = R[:,k]
		b = L[:,k] # changed this according to documentation 

		jac_dim = len(a)


		ab=np.zeros(jac_dim)

		for j in range(jac_dim):

			ab[j] = abs(a[j]*b[j])

		S = np.sum(ab)

		ei = ab/S

		scores[k] = np.dot(ei,ei_vec0)
	idx_CEM = np.argmax(scores)

	return idx_CEM

def follow_ei(D,L,R,ei_prev,idx_track):
	jac_dim = len(L)

	scores = np.zeros(jac_dim)
	err_rel = np.zeros(jac_dim)
	for k in range(jac_dim):

		a = R[:,k]
		b = L[:,k] # changed this according to documentation 

	
		ab=np.zeros(jac_dim)

		for j in range(jac_dim):

			ab[j] = abs(a[j]*b[j])

		S = np.sum(ab)

		ei = ab/S

		# scores[k] = np.dot(ei,ei_prev)
		# pdb.set_trace()
		err_rel[k] = np.sum(abs((ei-ei_prev)/ei_prev))


	# idx_followed = np.argmax(scores)
	idx_followed = np.argmin(err_rel)

	new_ei = EI(D,L,R,idx_followed)

	return new_ei, idx_followed


	

def test_parallel(matrix,v_test):

	# matrix contains column vectors to be tested against v_test vector for parallelism. 
	# returns the index (culumn) of vector in matrix most aligned to v_test
	scores = np.zeros(len(matrix[0,:]))

	for i in range(len(matrix[0,:])):

		scores[i] = np.dot(v_test,matrix[:,i])

	# them = scores[np.argsort(scores)[-10:]]
	# print them
	# pdb.set_trace()
	
	idx_most_aligned = np.argmax(scores)

	return idx_most_aligned


# def test_EI_parallel()



def match_eig(D,prev):

	# use slope to do difference between expected value and D[i]
	
	err_rel = np.zeros(len(D))
	for i in range(len(D)):

		err_rel[i] = abs(D[i]-prev)/abs(prev)

	eig = D[np.argmin(err_rel)]
	print("previous : \t {:f} \t next : {:f}".format(prev,eig))
	
	# pdb.set_trace()

	return eig

