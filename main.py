import numpy as np
import random
import itertools
import math

import pickle

from algorithms import generate_M, findIndpSubmatrix_even,  findIndpSubmatrix_odd



def IndpSubmatrix(n):

    M, bad_idx_1, bad_idx_2  = generate_M(n)

    bad_idx = bad_idx_1.union(bad_idx_2)

    indp_cols_lists = []
    used_bad_idx = set()

    if n%2 == 0: 

        # even case

        for N_1 in range(1,n//2):
    
            N_2 = n//4+N_1//2 -1

            if N_1 < N_2:
    
                indp_cols_list, M, used_bad_idx = findIndpSubmatrix_even( n, N_1, N_2, M, used_bad_idx, bad_idx_1 ,bad_idx_2)

                indp_cols_lists = indp_cols_lists + indp_cols_list

    if n%2 == 1:

        # odd case

        for N_2 in range(2, (n-1)//2 ):

            # for N_1 in range(1, N_2)

            N_1 = N_2//2

            indp_cols_list, M, used_bad_idx = findIndpSubmatrix_odd( n, N_1, N_2, M, used_bad_idx, bad_idx_1 ,bad_idx_2)

            indp_cols_lists = indp_cols_lists + indp_cols_list



    return indp_cols_lists


f = open( 'results.txt','a' )


results_dict = {}


for n in range(50,141):

    indp_cols_lists = []
    used_bad_idx = set()


    M, bad_idx_1, bad_idx_2  = generate_M(n)

    bad_idx = bad_idx_1.union(bad_idx_2)

    indp_cols_lists = IndpSubmatrix(n)

    results_dict[n] = indp_cols_lists

    #print(len( indp_cols_lists ))

    #print(indp_cols_lists)



    indp_cols_total_set = set()

    for indp_cols in indp_cols_lists:

        indp_cols_total_set  = indp_cols_total_set.union(set( indp_cols))

    M_block = M[ :, indp_cols]

    #print(np.linalg.matrix_rank(M_block))

    #print(len(indp_cols_total_set) )

    results_str = str(n) + ' \quad &' + str(len( indp_cols_lists )) + ' \\\\'+'\n'

    f.write(results_str)



with open("results_dict", "wb") as fp:   #Pickling
    
    pickle.dump(results_dict, fp)




