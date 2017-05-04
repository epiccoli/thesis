import numpy as np
import pdb
import os

def eigenmode(D,L,R):

	# D is vector of real parts of eigenvalues

	# Identify the index of the maximum eigenvalue
	max_idx = np.argmax(D)
	# Find the corresponding left and right eigenvectors and display
	left_mode = L[:,max_idx].real
	right_mode = R[:,max_idx].real


	# print left_mode, '\n'
	# print right_mode, '\n'


	# scores = np.zeros(len(L[0,:]))

	# for i in range(len(L[0,:])):
	# 	print 'left modes'

	# 	print L[:,i]

	# 	scores[i] = np.dot(left_mode,L[:,i])

	# print scores

	# idx_most_aligned = np.argmax(scores)

	# print idx_most_aligned

	# pdb.set_trace()

	return left_mode
	

def test_parallel(L,R,CEM):

	scores = np.zeros(len(L[0,:]))

	for i in range(len(L[0,:])):

		scores[i] = np.dot(CEM,L[:,i])

	# print scores
	idx_most_aligned = np.argmax(scores)
	them = scores[np.argsort(scores)[-3:]]
	print them
	pdb.set_trace()
	return idx_most_aligned