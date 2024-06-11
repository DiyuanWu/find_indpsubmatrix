import numpy as np 
import numpy as np
import random
import itertools


import math

def generate_M(n):
    M = np.zeros( (2*n+1, n**2 ) )

    for i in range(n):
    
        row1_M =  np.zeros((n,n))
        row1_M[i,:] = np.ones( (1,n) )
    
        tilde_I = np.identity(n)
        tilde_I = tilde_I[0:-1 , :]
    
        row_1 = np.zeros( (1,n))
        row_1[:,i] = 1
    
        row_2 = np.zeros( (1,n))
        
        row_2[:,-i-1] = 1
    
        block = np.block(
            [
                [ row1_M ],
                [ tilde_I ],
                [  row_1],
                
                [ row_2 ]
            ]
        )
    
        left_idx = n*i 
        right_idx = n*(i+1)
        M[ :, left_idx:right_idx] = block

    bad_idx_1 = set( [i*n + i for i in range(n)] )
    bad_idx_2 = set( [i*n + n-1-i for i in range(n)] )

    return M, bad_idx_1, bad_idx_2

def isFullrank( M ):

    max_rank = min(M.shape)
    
    return np.linalg.matrix_rank(M) == max_rank 


def findIndpSubmatrix_even( n, N_1, N_2, M, used_bad_idx, bad_idx_1, bad_idx_2):

    indp_cols_list = []

    bad_idx = bad_idx_1.union(bad_idx_2)


    for l in range(N_1,N_2):

        is_fullrank = False
        
        col_idx_1 = { i*n + ((i+1)+ 2*l-1)%n-1 for i in range(n) }
    
            
        col_idx_2 = { i*n + ((i+1)+ 2*l)%n-1 for i in range(n) }
    
        col_idx = col_idx_1.union(col_idx_2)      
    
        col_idx.remove( (n-2*l-1)*n + n-1 )
        
        col_idx = col_idx - bad_idx



        for j_1 in range(2*N_2 +1, int(math.ceil(n/2 + N_1)) ): 

            for j_2 in range(1,  N_1 ):

                temp_col_idx = col_idx
        
                new_col1_idx = n*(n//2-l) + j_1-1
                new_col2_idx = n*(n-l)+j_2-1
    
                temp_col_idx  = col_idx.union({new_col1_idx, new_col2_idx } )

                if np.linalg.matrix_rank( M[:, list(temp_col_idx) ] ) == 2*n -1:

                    col_idx = temp_col_idx

                    is_fullrank = True
                    
                    break

            if is_fullrank: 
                break    


        if is_fullrank:

            # print( len( col_idx ))
            # take two columns such that the last two rows are (1,0) and (0,1)
    
            #c_1 = list( bad_idx_1 - used_bad_idx )[0]
        
            #c_2 = list( bad_idx_2  - used_bad_idx )[0]

            c_1 = random.sample(list( bad_idx_1 - used_bad_idx ), 1)[0]
        
            c_2 = random.sample(list( bad_idx_2  - used_bad_idx ), 1)[0]
        
            col_idx  = col_idx.union({c_1, c_2 } )
            
            M_B = M[:,list(col_idx)]
    
            
        
            #print(M_B.shape)
        
            #print(M_B[(47,48),:])
        
            #print(np.linalg.matrix_rank(M_B) )
         
        
            if np.linalg.matrix_rank(M_B) == 2*n+1:
        
                #print(M_B)
        
                #print( col_idx )
        
                used_bad_idx = used_bad_idx.union({c_1,c_2})
        
                # print(np.linalg.matrix_rank( M[:,list(col_idx) ] ))
                
                M[:,list(col_idx)] =  np.zeros(M_B.shape)
    
                indp_cols_list.append(list(col_idx))
            

    return indp_cols_list, M, used_bad_idx


def findIndpSubmatrix_even_small( n,  M, used_bad_idx, bad_idx_1, bad_idx_2):

    indp_cols_list = []

    good_idx = set([i for i in range(n**2)]) - bad_idx_1

    good_idx = good_idx - bad_idx_2

    bad_idx = bad_idx_1.union(bad_idx_2)


    for l in range(1,n//2):

        is_fullrank = False
        
        col_idx_1 = { i*n + ((i+1)+ 2*l-1)%n-1 for i in range(n) }
    
            
        col_idx_2 = { i*n + ((i+1)+ 2*l)%n-1 for i in range(n) }
    
        col_idx = col_idx_1.union(col_idx_2) 
    
        col_idx.remove( (n-2*l-1)*n + n-1 )

    
        col_idx = col_idx - bad_idx

        #print(np.linalg.matrix_rank(M[:,list(col_idx)]))

        # enumerate all the possible indecies
        for j_1 in list(good_idx): 

            for j_2 in list(good_idx):

                if j_1 == j_2:
                    break;

                temp_col_idx = col_idx
        
                new_col1_idx = j_1 #n*(n//2-l) + j_1-1
                new_col2_idx = j_2 #n*(n-l)+j_2-1
        
                temp_col_idx  = col_idx.union({new_col1_idx, new_col2_idx } )

                if np.linalg.matrix_rank( M[:, list(temp_col_idx) ] ) == 2*n -1:

                    col_idx = temp_col_idx

                    is_fullrank = True
                    
                    break

            if is_fullrank: 
                break 




        if is_fullrank:

            # print( len( col_idx ))
            # take two columns such that the last two rows are (1,0) and (0,1)
    
            c_1 = random.sample(list( bad_idx_1 - used_bad_idx ), 1)[0]
        
            c_2 = random.sample(list( bad_idx_2  - used_bad_idx ), 1)[0]

        
            col_idx  = col_idx.union({c_1, c_2 } )
            
            M_B = M[:,list(col_idx)]
    

        
            if np.linalg.matrix_rank(M_B) == 2*n+1:
        
                #print(M_B)
        
                #print( col_idx )
        
                used_bad_idx = used_bad_idx.union({c_1,c_2})
        
                
                M[:,list(col_idx)] =  np.zeros(M_B.shape)
    
                indp_cols_list.append(list(col_idx))
            

    return indp_cols_list, M, used_bad_idx



def randSubmatrix(n, M, used_idx, num_trials=50 ):

    indp_cols_list = []

    curr_used_idx = used_idx

    unused_col_idx = set( [i*n+j for j in range(n) for i in range(n)] )

    unused_col_idx = unused_col_idx - set(curr_used_idx)
    
    for trial in range(num_trials):

        col_idx  = []

        num_col_block = 2
        
        for i in range(n):

            block_range = set([i*n + j for j in range(n)])

            block_range = block_range  - set(curr_used_idx)

            if len(block_range) >= 2:

                j_1, j_2 = random.sample( list(block_range), 2 )

                col_idx.append(j_1)

                col_idx.append(j_2)

            elif len(block_range) == 1:

                j_1 = list(block_range)[0]

                col_idx.append(j_1)

        if len(col_idx) < 2*n:

            num_rest = 2*n - len(col_idx)

            #print(num_rest, len(list(unused_col_idx)), len( col_idx) )

            new_col_list = list( random.sample( list(unused_col_idx - set( col_idx) ), num_rest ) )

            col_idx = col_idx+new_col_list
      
            


        if np.linalg.matrix_rank( M[ :, col_idx ] ) == 2*n:

            curr_unused_idx = set(unused_col_idx) - set(col_idx )

            for col in list(curr_unused_idx):

                col_idx_temp = col_idx + [col] 

                M_B = M[ :, col_idx_temp ]

                if np.linalg.matrix_rank( M_B ) == 2*n+1:

                    indp_cols_list.append( col_idx_temp )

                    curr_used_idx = curr_used_idx + col_idx_temp

                    unused_col_idx = unused_col_idx - set(curr_used_idx)

                    M[ :,col_idx_temp  ] = np.zeros( M_B.shape)

        if len(unused_col_idx) < 2*n+1:

            break
                  
    return  indp_cols_list, M


def findIndpSubmatrix_odd( n, N_1, N_2, M, used_bad_idx, bad_idx_1, bad_idx_2):

    indp_cols_list = []

    bad_idx = bad_idx_1.union(bad_idx_2)


    for l in range(N_1,N_2):

        is_fullrank = False
        
        col_idx_1 = { i*n + ((i+1)+ 2*l-1)%n-1 for i in range(n) }
    
            
        col_idx_2 = { i*n + ((i+1)+ 2*l)%n-1 for i in range(n) }
    
        col_idx = col_idx_1.union(col_idx_2)      
    
        col_idx.remove( (n-2*l-1)*n + n-1 )
        
        col_idx = col_idx - bad_idx



        for i_1 in range(1, int((n-1)/2+1-N_2) ): #

            for j_2 in range(1,N_1 ): #

                temp_col_idx = col_idx
        
                new_col1_idx = int((i_1-1)*n + (n-1)/2+l)
                new_col2_idx = n*(n-l)+j_2-1
    
                temp_col_idx  = col_idx.union({new_col1_idx, new_col2_idx } )

                if np.linalg.matrix_rank( M[:, list(temp_col_idx) ] ) == 2*n -1:

                    col_idx = temp_col_idx

                    is_fullrank = True
                    
                    break

            if is_fullrank: 
                break    


        if is_fullrank:

            # print( len( col_idx ))
            # take two columns such that the last two rows are (1,0) and (0,1)
    
            #c_1 = list( bad_idx_1 - used_bad_idx )[0]
        
            #c_2 = list( bad_idx_2  - used_bad_idx )[0]

            c_1 = random.sample(list( bad_idx_1 - used_bad_idx ), 1)[0]
        
            c_2 = random.sample(list( bad_idx_2  - used_bad_idx ), 1)[0]
        
            col_idx  = col_idx.union({c_1, c_2 } )
            
            M_B = M[:,list(col_idx)]
    
            
        
            #print(M_B.shape)
        
            #print(M_B[(47,48),:])
        
            #print(np.linalg.matrix_rank(M_B) )
         
        
            if np.linalg.matrix_rank(M_B) == 2*n+1:
        
                #print(M_B)
        
                #print( col_idx )
        
                used_bad_idx = used_bad_idx.union({c_1,c_2})
        
                # print(np.linalg.matrix_rank( M[:,list(col_idx) ] ))
                
                M[:,list(col_idx)] =  np.zeros(M_B.shape)
    
                indp_cols_list.append(list(col_idx))
            

    return indp_cols_list, M, used_bad_idx