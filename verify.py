import numpy as np
import random
import itertools
import math

import pickle

from algorithms import generate_M


with open("results_dict", "rb") as fp:   # Unpickling
   results_dict = pickle.load(fp)


for i, n in enumerate(results_dict):
   
    print('-------n='+str(n)+'-------\n')

    M, bad_idx_1, bad_idx_2 = generate_M(n)

    indp_cols_lists  = results_dict[n]

    indp_cols_total_set = set()

    for j,indp_cols in enumerate(indp_cols_lists):

        indp_cols_total_set  = indp_cols_total_set.union(set( indp_cols))

        M_block = M[ :, indp_cols]

        print('The '+str(j+1)+'-th independent columns:'+str(indp_cols)+'\n')

        print('The rank of the colums:' + str(np.linalg.matrix_rank(M_block))+'\n')

    print('The number of blocks times (2n+1): '+str(len(indp_cols_lists))+'*(2*'+str(n)+'+1)='+str(len(indp_cols_lists)*(2*n+1))+'\n')

    print('The number of columns found after removing potential duplicates:'+str(len(indp_cols_total_set ))+'\n')
