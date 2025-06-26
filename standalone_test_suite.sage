def polynomial_compression(E,degree,var_name):
    """
    Returns the corresponding polynomial matrix where the coefficient matrix of degree i is the ith block of ``E``.
    
    INPUT:
    
    - ``degree`` -- integer: the number of blocks in the expanded matrix. If ``E`` does not have zero-blocks,
        this corresponds to the degree of the resulting matrix.
    - ``variable`` -- a symbolic variable (e.g., ``x``) to use in the polynomial.
    
    OUTPUT:
    
    - a polynomial matrix with the same number of rows and E.ncols()//(degree+1) columns
    
    EXAMPLES:
    
    The first 3 x 3 block of N corresponds to the 
    
    sage: R.<x> = GF(97)[]
    sage: S.<y> = ZZ[]
    sage: M = matrix([[x^2+36*x,31*x,0],[3*x+13,x+57,0],[96,96,1]])
    sage: M
    [x^2 + 36*x       31*x          0]
    [  3*x + 13     x + 57          0]
    [        96         96          1]
    sage: N = M.expansion()
    sage: N
    [ 0  0  0 36 31  0  1  0  0]
    [13 57  0  3  1  0  0  0  0]
    [96 96  1  0  0  0  0  0  0]
    sage: polynomial_compression(N,2,x)
    [x^2 + 36*x       31*x          0]
    [  3*x + 13     x + 57          0]
    [        96         96          1]
    sage: polynomial_compression(N,2,y)
    [y^2 + 36*y       31*y          0]
    [  3*y + 13     y + 57          0]
    [        96         96          1]
    """
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    
    poly_ring = PolynomialRing(E.base_ring(),var_name)
    variable = poly_ring.gen()
    
    if E.ncols() % (degree+1) != 0:
        raise ValueError('The column number must be divisible by degree+1.')
    coeff_matrices = [E[:,i*(E.ncols()//(degree+1)):(i+1)*(E.ncols()//(degree+1))] for i in range(degree+1)]
    return sum([coeff_matrices[i]*variable**i for i in range(degree+1)])
    
def striped_krylov_matrix(E, J, degree, shift=None):
    r"""
    Return the block matrix of the following form, with rows permuted according to shift.
    The following uses `E` to refer to ``E``, and `d` to refer to ``degree``.
    
    [     ]
    [  E  ]
    [     ]
    [-----]
    [     ]
    [  EJ ]
    [     ]
    [-----]
    [  .  ]
    [  .  ]
    [  .  ]
    [-----]
    [     ]
    [ EJ^d]
    [     ]
    
    INPUT:
    
    - ``J`` -- a square matrix of size equal to the number of columns of ``E``.
    - ``degree`` -- integer, the maximum exponent for the Krylov matrix.
    - ``shift`` -- list of integers (optional), row priority shifts. If ``None``,
        defaults to all zero.
    
    OUTPUT:

    - A matrix with block rows [E, EJ, EJ^2, ..., EJ^d], row-permuted by shift.
    
    EXAMPLES:
    
    sage: R.<x> = GF(97)[]
    sage: E = matrix([[27,49,29],[50,58,0],[77,10,29]])
    sage: E
    [27 49 29]
    [50 58  0]
    [77 10 29]
    sage: J = matrix(R,[[0,1,0],[0,0,1],[0,0,0]])
    sage: J
    [0 1 0]
    [0 0 1]
    [0 0 0]
    sage: striped_krylov_matrix(E,J,3)
    [27 49 29]
    [50 58  0]
    [77 10 29]
    [ 0 27 49]
    [ 0 50 58]
    [ 0 77 10]
    [ 0  0 27]
    [ 0  0 50]
    [ 0  0 77]
    [ 0  0  0]
    [ 0  0  0]
    [ 0  0  0]
    sage: striped_krylov_matrix(E,J,3,[0,3,6])
    [27 49 29]
    [ 0 27 49]
    [ 0  0 27]
    [ 0  0  0]
    [50 58  0]
    [ 0 50 58]
    [ 0  0 50]
    [ 0  0  0]
    [77 10 29]
    [ 0 77 10]
    [ 0  0 77]
    [ 0  0  0]
    sage: striped_krylov_matrix(E,J,3,[3,0,2])
    [50 58  0]
    [ 0 50 58]
    [ 0  0 50]
    [77 10 29]
    [27 49 29]
    [ 0  0  0]
    [ 0 77 10]
    [ 0 27 49]
    [ 0  0 77]
    [ 0  0 27]
    [ 0  0  0]
    [ 0  0  0]
    
    """
    from sage.combinat.permutation import Permutation
    from sage.matrix.constructor import matrix # for matrix.block
    
    m = E.nrows()
    
    # if no shift, this is equivalent to 0s
    if shift is None:
        shift = [0]*E.nrows()
    
    # size checks
    if len(shift) != m:
        raise ValueError(f"The shift should have the same number of elements as the rows of the matrix ({m}), but had {len(shift)}.")
    if not J.is_square() or J.nrows() != E.ncols():
        raise ValueError(f"The matrix J should be a square matrix and match the number of columns of E ({E.ncols()}), but is of dimension {J.nrows()} x {J.ncols()}.")
    
    # define priority and indexing functions
    # index: [0,1,...,m-1,m,...,2m-1,...,m*degree-1] -> [(0,0),(1,0),...,(m-1,0),(0,1),...,(m-1,1),...,(m-1,degree)]
    priority = lambda c,d : shift[c] + d
    index = lambda i : (i%m,i//m)
    
    # deduce permutation by sorting priorities
    # rows should be ordered by ascending priority, then by associated row number in E (i%m)
    priority_triplets = sorted([[priority(*index(i)),index(i),i] for i in range(m*(degree+1))])
    priority_permutation = Permutation([t[2]+1 for t in priority_triplets])
    
    # build all blocks for the striped matrix
    blocks = []
    
    # for i from 0 to degree (inclusive), add E*J^i
    J_i = matrix.identity(E.ncols())
    for i in range(degree+1):
        blocks.append([E*J_i])
        if i < degree:
            J_i *= J
    
    # return block matrix permuted according to shift
    krylov = matrix.block(blocks,subdivide=False)
    krylov.permute_rows(priority_permutation)
    return krylov

def naive_krylov_rank_profile(E, J, degree, shift=None):
    r"""
    Compute the rank profile (row and column) of the striped Krylov matrix
    built from ``E`` and matrix `J`.

    INPUT:

    - ``J`` -- square matrix used in the Krylov construction.
    - ``degree`` -- maximum power of J in Krylov matrix.
    - ``shift`` -- list of integers (optional): priority row shift.

    OUTPUT:

    - A tuple (row_profile, pivot_matrix, column_profile):
        * row_profile: list of the first r row indices in the striped Krylov matrix `K` corresponding to independent rows
        * pivot: the submatrix of `K` given by those rows
        * column_profile: list of the first r independent column indices in ``pivot_matrix``
      where r is the rank of `K`.
    
    EXAMPLES:
    
    sage: R.<x> = GF(97)[]
    sage: E = matrix([[27,49,29],[50,58,0],[77,10,29]])
    sage: E
    [27 49 29]
    [50 58  0]
    [77 10 29]
    sage: J = matrix(R,[[0,1,0],[0,0,1],[0,0,0]])
    sage: J
    [0 1 0]
    [0 0 1]
    [0 0 0]
    sage: striped_krylov_matrix(E,J,3)
    [27 49 29]
    [50 58  0]
    [77 10 29]
    [ 0 27 49]
    [ 0 50 58]
    [ 0 77 10]
    [ 0  0 27]
    [ 0  0 50]
    [ 0  0 77]
    [ 0  0  0]
    [ 0  0  0]
    [ 0  0  0]
    sage: naive_krylov_rank_profile(E,J,3)
    (
               [27 49 29]
               [50 58  0]
    (0, 1, 3), [ 0 27 49], (0, 1, 2)
    )
    sage: striped_krylov_matrix(E,J,3,[0,3,6])
    [27 49 29]
    [ 0 27 49]
    [ 0  0 27]
    [ 0  0  0]
    [50 58  0]
    [ 0 50 58]
    [ 0  0 50]
    [ 0  0  0]
    [77 10 29]
    [ 0 77 10]
    [ 0  0 77]
    [ 0  0  0]
    sage: naive_krylov_rank_profile(E,J,3,[0,3,6])
    (
               [27 49 29]
               [ 0 27 49]
    (0, 1, 2), [ 0  0 27], (0, 1, 2)
    )
    sage: striped_krylov_matrix(E,J,3,[3,0,2])
    [50 58  0]
    [ 0 50 58]
    [ 0  0 50]
    [77 10 29]
    [27 49 29]
    [ 0  0  0]
    [ 0 77 10]
    [ 0 27 49]
    [ 0  0 77]
    [ 0  0 27]
    [ 0  0  0]
    [ 0  0  0]
    sage: naive_krylov_rank_profile(E,J,3,[3,0,2])
    (
               [50 58  0]
               [ 0 50 58]
    (0, 1, 2), [ 0  0 50], (0, 1, 2)
    )
    """
    from sage.matrix.constructor import matrix
    
    K = striped_krylov_matrix(E,J,degree,shift)
    
    # calculate first independent rows of K
    row_profile = K.transpose().pivots()
    
    # construct submatrix
    pivot = K.matrix_from_rows(row_profile)
    
    # calculate first independent columns of pivot
    col_profile = pivot.pivots()
    
    return row_profile,pivot,col_profile

def krylov_rank_profile(E, J, degree, shift=None, timer=[]):
    r"""
    Compute the rank profile (row and column) of the striped Krylov matrix
    built from ``E`` and matrix `J`.

    INPUT:

    - ``J`` -- square matrix used in the Krylov construction.
    - ``degree`` -- maximum power of J in Krylov matrix.
    - ``shift`` -- list of integers (optional): priority row shift.

    OUTPUT:

    - A tuple (row_profile, pivot_matrix, column_profile):
        * row_profile: list of the first r row indices in the striped Krylov matrix `K` corresponding to independent rows
        * pivot: the submatrix of `K` given by those rows
        * column_profile: list of the first r independent column indices in ``pivot_matrix``
      where r is the rank of `K`.
    
    EXAMPLES:
    
    sage: R.<x> = GF(97)[]
    sage: E = matrix([[27,49,29],[50,58,0],[77,10,29]])
    sage: E
    [27 49 29]
    [50 58  0]
    [77 10 29]
    sage: J = matrix(R,[[0,1,0],[0,0,1],[0,0,0]])
    sage: J
    [0 1 0]
    [0 0 1]
    [0 0 0]
    sage: striped_krylov_matrix(E,J,3)
    [27 49 29]
    [50 58  0]
    [77 10 29]
    [ 0 27 49]
    [ 0 50 58]
    [ 0 77 10]
    [ 0  0 27]
    [ 0  0 50]
    [ 0  0 77]
    [ 0  0  0]
    [ 0  0  0]
    [ 0  0  0]
    sage: krylov_rank_profile(E,J,3)
    (
               [27 49 29]
               [50 58  0]
    (0, 1, 3), [ 0 27 49], (0, 1, 2)
    )
    sage: striped_krylov_matrix(E,J,3,[0,3,6])
    [27 49 29]
    [ 0 27 49]
    [ 0  0 27]
    [ 0  0  0]
    [50 58  0]
    [ 0 50 58]
    [ 0  0 50]
    [ 0  0  0]
    [77 10 29]
    [ 0 77 10]
    [ 0  0 77]
    [ 0  0  0]
    sage: krylov_rank_profile(E,J,3,[0,3,6])
    (
               [27 49 29]
               [ 0 27 49]
    (0, 1, 2), [ 0  0 27], (0, 1, 2)
    )
    sage: striped_krylov_matrix(E,J,3,[3,0,2])
    [50 58  0]
    [ 0 50 58]
    [ 0  0 50]
    [77 10 29]
    [27 49 29]
    [ 0  0  0]
    [ 0 77 10]
    [ 0 27 49]
    [ 0  0 77]
    [ 0  0 27]
    [ 0  0  0]
    [ 0  0  0]
    sage: krylov_rank_profile(E,J,3,[3,0,2])
    (
               [50 58  0]
               [ 0 50 58]
    (0, 1, 2), [ 0  0 50], (0, 1, 2)
    )
    """
    from sage.matrix.constructor import matrix
    import math
    import time
    from sage.combinat.permutation import Permutation
    
    start_time = time.time()
    
    #DEBUG 0
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    m = E.nrows()
    
    if m == 0:
        return (),E,()
    
    #DEBUG 1
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # calculate shift priority function and permutation
    priority = lambda c,d : shift[c] + d
    index = lambda i : (i%m,i//m)
    index_inv = lambda c,d : c + d*m
    
    priority_triplets = sorted([[priority(*index(i)),index(i),i] for i in range(m*(degree+1))])
    priority_permutation = Permutation([t[2]+1 for t in priority_triplets])
    priority_permutation_inv = priority_permutation.inverse()
    
    # maps row c of E*J^d to position i in striped Krylov matrix
    # +/- 1 as permutations are 1-indexed, matrices are 0-indexed
    phi = lambda c,d : priority_permutation_inv(index_inv(c,d) + 1) - 1
    phi_inv = lambda i : index(priority_permutation(i + 1) - 1)
    
    #DEBUG 2
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # calculate row profile of E, with shift applied
    self_permutation = Permutation([pair[1]+1 for pair in sorted([[phi(i,0),i] for i in range(m)])])
    
    #DEBUG 3
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    row_profile_self_permuted = E.with_permuted_rows(self_permutation).transpose().pivots()
    row_profile_self = [self_permutation(i+1)-1 for i in row_profile_self_permuted]
    #DEBUG 4
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # base row_profile_K using row_profile_self
    row_profile_K = sorted([phi(i,0) for i in row_profile_self])
    r = len(row_profile_K)
    
    #DEBUG 5
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    if r == 0:
        return (),E.matrix_from_rows([]),()
    
    M = E.matrix_from_rows([index_inv(*phi_inv(i)) for i in row_profile_K])
    
    #DEBUG 6
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    J_L = None
    for l in range(math.ceil(math.log(degree,2)) + 1):
        # row and degree of profile within unshifted expansion
        c,d = zip(*(phi_inv(i) for i in row_profile_K))
        
        #DEBUG 7
        dbg_num = 7
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        L = pow(2,l)
        # adding 2^l to each degree
        row_extension = sorted([phi(c[i],d[i] + L) for i in range(r) if d[i] + L <= degree])
        
        #DEBUG 8
        dbg_num = 8
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        # concatenate two sequences (order preserved)
        k = row_profile_K + row_extension
        # calculate sorting permutation
        k_perm = Permutation([x[1] for x in sorted([k[i],i+1] for i in range(len(k)))])
        
        #DEBUG 9
        dbg_num = 9
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        # fast calculation of rows formed by indices in k
        if J_L is None:
            J_L = J
        else:
            J_L = J_L * J_L
        
        #DEBUG 10
        dbg_num = 10
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        M = matrix.block([[M],[M*J_L]],subdivide=False)
        
        #DEBUG 11
        dbg_num = 11
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        # sort rows of M, find profile, translate to k (indices of full krylov matrix)
        M.permute_rows(k_perm)
        
        #DEBUG 12
        dbg_num = 12
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        M_t = M.transpose()
        
        #DEBUG 13
        dbg_num = 13
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        row_profile_M = M_t.pivots()
        
        #DEBUG 14
        dbg_num = 14
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        r = len(row_profile_M)
        row_profile_K = [k[k_perm(i+1)-1] for i in row_profile_M]
        
        #DEBUG 15
        dbg_num = 15
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
        
        # calculate new M for return value or next loop
        M = matrix([M.row(i) for i in row_profile_M])
        
        #DEBUG 16
        dbg_num = 16
        time_time = time.time()
        if len(timer) <= dbg_num:
            timer.append(time_time - start_time)
        else:
            time_diff = time.time() - start_time - timer[-1]
            for i in range(dbg_num,len(timer)):
                timer[i] += time_diff
        #DEBUG
    #DEBUG 17
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    col_profile = M.pivots()
    #DEBUG 18
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    return tuple(row_profile_K), M, col_profile

def linear_interpolation_basis(E, J, degree, var_name, shift=None, timer=[]):
    r"""
    Construct a linear interpolant basis for (``E``,`J`) in `s`-Popov form.

    INPUT:

    - ``J`` -- square matrix for Krylov construction.
    - ``degree`` -- upper bound on degree of minpoly(`J`), power of 
    - ``shift`` -- list of integers (optional): controls priority of rows.

    OUTPUT:

    - Matrix whose rows form an interpolation basis.

    """
    import time
    from sage.combinat.permutation import Permutation
    from sage.matrix.constructor import matrix
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
    
    start_time = time.time()
    
    #DEBUG 0
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    poly_ring = PolynomialRing(E.base_ring(),var_name)
    
    m = E.nrows()
    
    # if no shift, this is equivalent to 0s
    if shift is None:
        shift = [0]*E.nrows()
    
    #DEBUG 1
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # calculate shift priority function and permutation
    priority = lambda c,d : shift[c] + d
    index = lambda i : (i%m,i//m)
    index_inv = lambda c,d : c + d*m
    
    priority_triplets = sorted([[priority(*index(i)),index(i),i] for i in range(m*(degree+1))])
    priority_permutation = Permutation([t[2]+1 for t in priority_triplets])
    priority_permutation_inv = priority_permutation.inverse()
    
    # maps row c of EJ^d to position i in striped Krylov matrix
    # +/- 1 as permutations are 1-indexed, matrices are 0-indexed
    phi = lambda c,d : priority_permutation_inv(index_inv(c,d) + 1) - 1
    phi_inv = lambda i : index(priority_permutation(i + 1) - 1)
    
    #DEBUG 2
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # calculate krylov profile
    row_profile, pivot, col_profile = krylov_rank_profile(E,J,degree,shift)
    
    #DEBUG 3
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    if len(row_profile) == 0:
        return matrix.identity(poly_ring,m)
    
    #DEBUG 4
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # (c_k, d_k) = phi^-1 (row_i)
    c, d = zip(*(phi_inv(i) for i in row_profile))
    
    #DEBUG 5
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # degree_c = 0 or max d[k] such that c[k] = c
    inv_c = [None]*m
    degree_c = [0]*m
    for k in range(len(row_profile)):
        if d[k] >= degree_c[c[k]]:
            degree_c[c[k]] = d[k] + 1
            inv_c[c[k]] = k
    
    #DEBUG 6
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    T_absent_indices = [i for i in range(m) if inv_c[i] is None]
    T_present_indices = [inv_c[i] for i in range(m) if inv_c[i] is not None]
    D_absent = E.matrix_from_rows_and_columns(T_absent_indices, col_profile)
    D_present = (pivot.matrix_from_rows(T_present_indices)*J).matrix_from_columns(col_profile)
    D_rows = []
    idx_p = 0
    idx_a = 0
    for i in range(m):
       if inv_c[i] is None:
           D_rows.append(D_absent[idx_a])
           idx_a += 1
       else:
           D_rows.append(D_present[idx_p])
           idx_p += 1
           
    C = pivot.matrix_from_columns(col_profile)
    D = matrix(D_rows)
    
    #DEBUG 7
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    relation = D*C.inverse()
    
    #DEBUG 8
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    # linear interpolation basis in shifted Popov form
    #uncompressed_basis = matrix.block([[-relation,matrix.identity(m)]],subdivide=False)
    
    # construct variable
    variable = poly_ring.gen()
    
    basis_rows = [[[] for j in range(m)] for i in range(m)]
    
    # compression of basis into polynomial form
    for i in range(m):
        basis_rows[i][i] = [0]*degree_c[i] + [1]
    
    #DEBUG 9
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    for col in range(relation.ncols()):
        for row in range(m):
            if len(basis_rows[row][c[col]]) <= d[col]:
                basis_rows[row][c[col]] += [0]*(d[col]+1-len(basis_rows[row][c[col]]))
            basis_rows[row][c[col]][d[col]] -= relation[row][col]
    
    #DEBUG 10
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    output = matrix(poly_ring, basis_rows)
    
    #DEBUG 11
    time_time = time.time()
    timer.append(time_time - start_time)
    #DEBUG
    
    return output

import random
import time

class KrylovTestInstance:
    def manual_init(self, field, m, sigma, E, J, shift, degree):
        # parameter basic validation
        if not(isinstance(field,sage.rings.ring.Field)):
            raise ValueError('field is not a valid field.')
        if not(isinstance(m,int) and 0 <= m):
            raise ValueError('m is not a non-negative integer.')
        if not(isinstance(sigma,int) and 0 <= sigma):
            raise ValueError('sigma is not a non-negative integer.')
        if not(isinstance(E,sage.structure.element.Matrix)):
            raise ValueError('E is not a matrix.')
        if not(isinstance(J,sage.structure.element.Matrix)):
            raise ValueError('J is not a matrix.')
        if shift is not None and (not(isinstance(shift,list)) or any(not(isinstance(elem,int)) for elem in shift)):
            raise ValueError('shift is not a list of integers.')
        if not(isinstance(degree,int) and 1 <= degree and degree & (degree -1) == 0):
            raise ValueError('sigma is not a integer that is a power of two.')

        # field validation
        if E.base_ring() != field:
            raise ValueError('E does not have the correct field.')
        if J.base_ring() != field:
            raise ValueError('J does not have the correct field.')

        # dimension validation
        if E.nrows() != m or E.ncols() != sigma:
            raise ValueError('E does not have the correct dimensions.')
        if J.nrows() != sigma or J.ncols() != sigma:
            raise ValueError('J does not have the correct dimensions.')
        if shift is not None and len(shift) != m:
            raise ValueError('shift does not have the correct length.')

        # field assignment
        self.field = field
        self.m = m
        self.sigma = sigma
        self.E = E
        self.J = J
        self.shift = shift
        self.degree = degree

    def __init__(self, field, m, sigma, options={'E_mode': None, 'J_mode': None, 'shift_mode': 'uniform'}):
        # base parameter validation
        if not(isinstance(field,sage.rings.ring.Field)):
            raise ValueError('field is not a valid field.')
        if not(isinstance(m,int) and 0 <= m):
            raise ValueError('m is not a non-negative integer.')
        if not(isinstance(sigma,int) and 0 <= sigma):
            raise ValueError('sigma is not a non-negative integer.' + str(sigma))
            
        # manual test case
        if 'manual' in options:
            if len(options) != 1:
                raise ValueError('options has extra arguments.')
            manual = options['manual']
            if 'E' not in manual or 'J' not in manual or 'shift' not in manual or 'degree' not in manual:
                raise ValueError('options is ill-formatted, does not contain E, J, shift, and degree.')
            self.manual_init(field, m, sigma, manual['E'], manual['J'], manual['shift'], manual['degree'])
            return

        # generated test case parameter validation
        if 'E_mode' not in options or 'J_mode' not in options or 'shift_mode' not in options:
            raise ValueError('options is ill-formatted, does not contain E_mode, J_mode, and shift_mode.')
        if options['E_mode'] is not None:
            mode = options['E_mode']
            if 'mode' not in mode or mode['mode'] not in ['density','rank']:
                raise ValueError('E_mode does not contain a valid mode.')
            if mode['mode'] == 'density' and not('density' in mode and isinstance(mode['density'],float) and 0 <= mode['density'] <= 1):
                raise ValueError('E_mode does not contain a valid density.')
            if mode['mode'] == 'rank' and not('rank' in mode and isinstance(mode['rank'],int) and 0 <= mode['rank'] <= min(m,sigma)):
                raise ValueError('E_mode does not contain a valid rank.')
        if options['J_mode'] is not None:
            mode = options['J_mode']
            if 'mode' not in mode or mode['mode'] not in ['density','rank','jordan']:
                raise ValueError('J_mode does not contain a valid mode.')
            if mode['mode'] == 'density' and not('density' in mode and isinstance(mode['density'],float) and 0 <= mode['density'] <= 1):
                raise ValueError('J_mode does not contain a valid density.')
            if mode['mode'] == 'rank' and not('rank' in mode and isinstance(mode['rank'],int) and 0 <= mode['rank'] <= sigma):
                raise ValueError('J_mode does not contain a valid rank.')
            if mode['mode'] == 'jordan' and not('block_number' in mode and isinstance(mode['block_number'],int) and 1 <= mode['block_number'] and mode['block_number'] <= sigma):
                raise ValueError('J_mode does not contain a valid block_number.')
        if options['shift_mode'] not in ['uniform','increasing','decreasing','random','hermite','reverse_hermite','plateau','reverse_plateau']:
            raise ValueError('shift_mode is not valid.')
        
        # field assignment
        self.field = field
        self.m = m
        self.sigma = sigma
        self.options = options
        
        # E generation
        if options['E_mode'] is None:
            self.E = matrix.random(field,m,sigma)
        else:
            if options['E_mode']['mode'] == 'density':
                vec = matrix.random(field,1,m*sigma,density=options['E_mode']['density'])
                self.E = matrix(field,m,sigma,vec.list())
            if options['E_mode']['mode'] == 'rank':
                rank = options['E_mode']['rank']
                self.E = matrix.random(field,m,rank) * matrix.random(field,rank,sigma)
            
        # J generation
        if options['J_mode'] is None:
            self.J = matrix.random(field,sigma,sigma)
        else:
            if options['J_mode']['mode'] == 'density':
                vec = matrix.random(field,1,sigma*sigma,density=options['J_mode']['density'])
                self.J = matrix(field,sigma,sigma,vec.list())
            if options['J_mode']['mode'] == 'rank':
                rank = options['J_mode']['rank']
                self.J = matrix.random(field,sigma,rank) * matrix.random(field,rank,sigma)
            if options['J_mode']['mode'] == 'jordan':
                block_boundaries = [0]+sorted(random.sample(range(1,sigma),options['J_mode']['block_number']-1))+[sigma]
                blocks = [matrix.jordan_block(field.random_element(),block_boundaries[i+1]-block_boundaries[i]) for i in range(options['J_mode']['block_number'])]
                self.J = matrix.block_diagonal(blocks)
        
        # shift generation
        if options['shift_mode'] == 'uniform':
            self.shift = [0]*m
        if options['shift_mode'] == 'increasing':
            self.shift = [i for i in range(m)]
        if options['shift_mode'] == 'decreasing':
            self.shift = [i for i in reversed(range(m))]
        if options['shift_mode'] == 'random':
            self.shift = [i for i in range(m)]
            shuffle(self.shift)
        if options['shift_mode'] == 'hermite':
            self.shift = [i*m*sigma for i in range(m)]
        if options['shift_mode'] == 'reverse_hermite':
            self.shift = [i*m*sigma for i in reversed(range(m))]
        if options['shift_mode'] == 'plateau':
            self.shift = [0]*(m//2) + [m*sigma]*(m - m//2)
        if options['shift_mode'] == 'reverse_plateau':
            self.shift = [m*sigma]*(m//2) + [0]*(m - m//2)

        # degree assignment
        if sigma == 0:
            self.degree = 1
        elif sigma & (sigma-1) == 0:
            self.degree = sigma
        else:
            self.degree = sigma
            while self.degree & (self.degree+1) != 0:
                self.degree |= self.degree >> 1
            self.degree += 1

    def generator(self):
        # field
        if self.field.is_finite():
            field_string = f"GF({self.field.order()})"
        elif self.field == QQ:
            field_string = f"QQ"
        else:
            raise ValueError('Field unsupported: ' + str(self.field))
        
        # matrices
        if self.field.is_finite() and self.field.is_prime_field():
            E_string = f"matrix({field_string},{self.m},{self.sigma},{self.E.list()})"
            J_string = f"matrix({field_string},{self.sigma},{self.sigma},{self.J.list()})"
        if self.field.is_finite() and not(self.field.is_prime_field()):
            E_string = f"matrix({field_string},{self.m},{self.sigma},[{','.join([f"{field_string}({elem.list()})" for elem in self.E.list()])}])"
            J_string = f"matrix({field_string},{self.sigma},{self.sigma},[{','.join([f"{field_string}({elem.list()})" for elem in self.J.list()])}])"
        if self.field == QQ:
            E_string = f"matrix({field_string},{self.m},{self.sigma},{self.E.list()})"
            J_string = f"matrix({field_string},{self.sigma},{self.sigma},{self.J.list()})"
        
        # shift
        shift_string = f"[{','.join([f"int({i})" for i in self.shift])}]"
        
        return f"KrylovTestInstance({field_string},int({self.m}),int({self.sigma}),options={{'manual':{{'E':{E_string},'J':{J_string},'shift':{shift_string},'degree':int({self.degree})}}}})"

    def naive_krylov_rank_profile(self):
        return naive_krylov_rank_profile(self.E,self.J,self.degree,self.shift)

    def krylov_rank_profile(self):
        return krylov_rank_profile(self.E,self.J,self.degree,self.shift)

    def linear_interpolation_basis(self):
        return linear_interpolation_basis(self.E,self.J,self.degree,'x',self.shift)

def generate_test_set():
    m_set = [int(i) for i in [0,1,2,3,4,5,6,7,8]]
    sig_set = [int(i) for i in [0,1,2,3,4,5,6,7,8]]
    
    field_set = [GF(q) for q in range(100) if is_prime(q)] + [GF(q**2) for q in range(16) if is_prime(q)] + [GF(q**3) for q in [2,3,5]] + [GF(3**i) for i in [4,5]] + [GF(2**i) for i in [4,5,6,7,8]] + [QQ]
    
    shift_set = ['uniform','increasing','decreasing','random','hermite','reverse_hermite','plateau','reverse_plateau']
    
    density_set = [float(i) for i in [0.01,0.1,0.5,0.99,1.0]]
    
    test_set = []

    for sig in sig_set:
        J_ranks = random.sample(range(sig+1),2) if sig > 1 else []
        jordan_numbers = range(1,max(sig//4+1,2)) if sig > 0 else []

        J_set = []
        if sig > 1:
            for d in density_set:
                J_set.append({'mode':'density','density':d})

            for r in J_ranks:
                J_set.append({'mode':'rank','rank':r})

            for j in jordan_numbers:
                J_set.append({'mode':'jordan','block_number':j})
        elif sig == 0:
            J_set.append({'mode':'density','density':float(1)})
        else:
            J_set.append({'mode':'density','density':float(0)})
            J_set.append({'mode':'density','density':float(1)})

        for m in m_set:
            E_ranks = random.sample(range(min(m,sig)+1),2) if min(m,sig) > 1 else []
            E_set = []
            
            if m > 1:
                for d in density_set:
                    E_set.append({'mode':'density','density':d})

                for r in E_ranks:
                    E_set.append({'mode':'rank','rank':r})
            elif m == 0:
                E_set.append({'mode':'density','density':float(1)})
            else:
                E_set.append({'mode':'density','density':float(0)})
                E_set.append({'mode':'density','density':float(1)})

            for f in field_set:
                for E_opt in E_set:
                    for J_opt in J_set:
                        for shift in shift_set:
                            if random.random() < 0.001: # ~20 seconds
                                test_set.append(KrylovTestInstance(f,m,sig,{'E_mode':E_opt,'J_mode':J_opt,'shift_mode':shift}))

    return test_set

def is_profile_correct(test,profile=None):
    if profile is None:
        profile = test.krylov_rank_profile()
    return profile == test.naive_krylov_rank_profile()

def is_basis_correct(test,basis=None):
    if basis is None:
        basis = test.linear_interpolation_basis()
    if not basis.is_popov(shifts=test.shift):
        return False
    if basis.nrows()*basis.ncols() != 0 and basis.degree() > test.sigma:
        return False
    priority = lambda c,d : test.shift[c] + d
    index = lambda i : (i%test.m,i//test.m)
    priority_triplets = sorted([[priority(*index(i)),index(i),i] for i in range(test.m*(test.degree+1))])
    priority_permutation = Permutation([t[2]+1 for t in priority_triplets])
    try:
        if basis.expansion(test.degree).with_permuted_columns(priority_permutation)*striped_krylov_matrix(test.E,test.J,test.degree,test.shift) != matrix.zero(test.field,test.m,test.sigma):
            return False
    except:
        return False
    return True

def run_test_set(tests):
    success = [{'profile_completed':False,'basis_completed':False,'profile_correct':False,'basis_correct':False} for i in range(len(tests))]
    start = time.time()
    for i in range(len(tests)):
        test = tests[i]
        profile = None
        basis = None
        attempts = 3
        # profile terminates
        for attempt in range(attempts):
            try:
                profile = test.krylov_rank_profile()
                success[i]['profile_completed'] = True
                break
            except:
                if attempt == attempts - 1:
                    success[i]['profile_completed'] = False
        if not success[i]['profile_completed']:
            continue
        # profile is correct
        success[i]['profile_correct'] = is_profile_correct(test,profile)
        # basis terminates
        for attempt in range(attempts):
            try:
                basis = test.linear_interpolation_basis()
                success[i]['basis_completed'] = True
                break
            except:
                if attempt == attempts - 1:
                    success[i]['basis_completed'] = False
        if not success[i]['basis_completed']:
            continue
        # basis is correct
        success[i]['basis_correct'] = is_basis_correct(test,basis)
    running_time = time.time() - start
    correct = sum([sum(result.values()) for result in success])
    total = sum([len(result) for result in success])
    ratio = (float(100)*float(correct))/float(total)
    print(f"Tests passed: {correct} out of {total} ({ratio}%)")
    print(f"Completed in {running_time:.2f}s")
    pc = 5
    bc = 5
    pmn = 5
    bcr = 5
    for i in range(len(tests)):
        failure = False
        if pc > 0 and not success[i]['profile_completed']:
            pc -= 1
            failure = True
        if bc > 0 and not success[i]['basis_completed']:
            bc -= 1
            failure = True
        if pmn > 0 and not success[i]['profile_correct']:
            pmn -= 1
            failure = True
        if bcr > 0 and not success[i]['basis_correct']:
            bcr -= 1
            failure = True
        if failure and tests[i].m*tests[i].sigma != 0:
            print(i)
            print(success[i])
            print(tests[i].generator())

def compare():
    for m in range(10,110,10):
        for sig in range(10,110,10):
            if sig != m:
                continue
            test = KrylovTestInstance(GF(97),int(m),int(sig),{'E_mode':None,'J_mode':None,'shift_mode':'decreasing'})
            start = time.time()
            for i in range(10):
                naive_krylov_rank_profile(test.E,test.J,test.degree,test.shift)
            stop = time.time()
            naive = (stop-start)/10
            start = time.time()
            for i in range(10):
                krylov_rank_profile(test.E,test.J,test.degree,test.shift)
            stop = time.time()
            fast = (stop-start)/10
            print(m,sig)
            print(naive,fast)

def measure_once():
    basis = False
    glob = None
    it = 1
    num_worst = 5
    for i in range(it):
        test = KrylovTestInstance(GF(2**8),int(4096),int(4096),{'E_mode':{'mode':'rank','rank':int(3072)},'J_mode':{'mode':'rank','rank':int(3072)},'shift_mode':'random'})
        timer = []
        if basis:
            linear_interpolation_basis(test.E,test.J,test.degree,'x',test.shift,timer)
        else:
            krylov_rank_profile(test.E,test.J,test.degree,test.shift,timer)
        if glob is None:
            glob = timer
        else:
            for j in range(len(timer)):
                glob[j] += timer[j]
    for i in range(len(glob)):
        glob[i] /= it
    #print(timer[3]-timer[2])
    #print(timer[11]-timer[10])
    print(timer)
    for pair in sorted([(timer[i+1] - timer[i],i) for i in range(len(timer)-1)])[-num_worst:]:
        print(f"between {pair[1]} and {pair[1]+1}: {pair[0]}s")

'''def measure_fields():
    
        test = KrylovTestInstance(QQ,int(64),int(64),{'E_mode':{'mode':'rank','rank':int(48)},'J_mode':{'mode':'rank','rank':int(48)},'shift_mode':'random'})
        timer = []
        linear_interpolation_basis(test.E,test.J,test.degree,'x',test.shift,timer)
        print(timer)'''

def measure():
    m_range = [int(2**i) for i in range(11)]
    s_range = [int(2**i) for i in range(11)]
    q_range = [i for i in range(16) if is_prime(i)]
    p_range = [i+1 for i in range(10)]
    trials = int(10)
    
    """for m in m_range:
        field = GF(256)
        sig = int(128)
        tests = [KrylovTestInstance(field,m,sig,{'E_mode':{'mode':'rank','rank':min(m,sig)//int(2)},'J_mode':{'mode':'rank','rank':sig//int(2)},'shift_mode':'random'}) for i in range(trials)]
        start_time = time.time()
        for test in tests:
            test.krylov_rank_profile()
        mid_time = time.time()
        for test in tests:
            test.linear_interpolation_basis()
        end_time = time.time()
        print(f"m: {m},sigma: {sig},field: GF({256}), krylov_time: {(mid_time-start_time)/trials:.3f}s, basis_time: {(end_time-mid_time)/trials:.3f}s")
    
    for sig in s_range:
        field = GF(256)
        m = int(128)
        tests = [KrylovTestInstance(field,m,sig,{'E_mode':{'mode':'rank','rank':min(m,sig)//int(2)},'J_mode':{'mode':'rank','rank':sig//int(2)},'shift_mode':'random'}) for i in range(trials)]
        start_time = time.time()
        for test in tests:
            test.krylov_rank_profile()
        mid_time = time.time()
        for test in tests:
            test.linear_interpolation_basis()
        end_time = time.time()
        print(f"m: {m},sigma: {sig},field: GF({256}), krylov_time: {(mid_time-start_time)/trials:.3f}s, basis_time: {(end_time-mid_time)/trials:.3f}s")
    
    for q in q_range:
        field = GF(q**5)
        sig = int(128)
        m = int(128)
        tests = [KrylovTestInstance(field,m,sig,{'E_mode':{'mode':'rank','rank':min(m,sig)//int(2)},'J_mode':{'mode':'rank','rank':sig//int(2)},'shift_mode':'random'}) for i in range(trials)]
        start_time = time.time()
        for test in tests:
            test.krylov_rank_profile()
        mid_time = time.time()
        for test in tests:
            test.linear_interpolation_basis()
        end_time = time.time()
        print(f"m: {m},sigma: {sig},field: GF({q**5}), krylov_time: {(mid_time-start_time)/trials:.3f}s, basis_time: {(end_time-mid_time)/trials:.3f}s")"""
    
    for p in p_range:
        field = GF(2**p)
        sig = int(128)
        m = int(128)
        tests = [KrylovTestInstance(field,m,sig,{'E_mode':{'mode':'rank','rank':min(m,sig)//int(2)},'J_mode':{'mode':'rank','rank':sig//int(2)},'shift_mode':'random'}) for i in range(trials)]
        start_time = time.time()
        for test in tests:
            test.krylov_rank_profile()
        mid_time = time.time()
        for test in tests:
            test.linear_interpolation_basis()
        end_time = time.time()
        print(f"m: {m},sigma: {sig},field: GF({2**p}), krylov_time: {(mid_time-start_time)/trials:.3f}s, basis_time: {(end_time-mid_time)/trials:.3f}s")
    
    '''for q in q_range:
        for p in p_range:
            field = GF(q**p)
            for m in m_range:
                for sig in s_range:
                    tests = [KrylovTestInstance(field,m,sig,{'E_mode':{'mode':'rank','rank':min(m,sig)//int(2)},'J_mode':{'mode':'rank','rank':sig//int(2)},'shift_mode':'random'}) for i in range(trials)]
                    start_time = time.time()
                    for test in tests:
                        test.krylov_rank_profile()
                    mid_time = time.time()
                    for test in tests:
                        test.linear_interpolation_basis()
                    end_time = time.time()
                    print(f"m: {m},sigma: {sig},field: GF({q**p}), krylov_time: {(mid_time-start_time)/trials:.3f}s, basis_time: {(end_time-mid_time)/trials:.3f}s")'''
    '''test = KrylovTestInstance(GF(97),int(200),int(200),{'E_mode':None,'J_mode':None,'shift_mode':'decreasing'})
    timing = [0]*20
    krylov_rank_profile(test.E,test.J,test.degree,test.shift)
    print(timing)
    start = time.time()
    linear_interpolation_basis(test.E,test.J,test.degree,test.field[x].gen(),test.shift)
    print(time.time() - start)'''

def reproduce_seg_fault():
    matrix.zero(GF(4),2,0).with_permuted_rows(Permutation([2,1])) # fixed

def fast_perms(shift,deg,c,d):
    # shift[c] - shift[i] + d - number of rows preceding with strictly smaller d
    # (i < c and shift[i] <= shift[c] + d) - +1 if smaller c but degree matches

    return sum(min(max(shift[c] - shift[i] + d + (i < c and shift[i] <= shift[c] + d),0),deg+1) for i in range(len(shift)))

def fast_perms_inv(shift,deg,i):
    # 
    
    # precompute
    n = len(shift)
    min_undepleted_row = 0
    # each entry for a row_c in E is:
    # shift[c]
    # c
    # (internal counter) the number of row_c*J^d counted already in the striped Krylov matrix
    # the number of rows in the striped matrix with priority < shift[c]
    row_data = sorted([shift[j],j,0,0] for j in range(n))
    for j in range(n-1):
        # initialise rows with priority < shift[c] as those with priority < shift[c-1]
        # then calculate the difference
        row_data[j+1][3] = row_data[j][3]
        # how many degrees further can the lowest priority shift go?
        levels_allowed = deg + 1 - row_data[min_undepleted_row][2]
        # what is the degree difference to correct for?
        levels_needed = row_data[j+1][0] - row_data[j][0]
        # while correction is needed
        while levels_needed > 0:
            # number of degrees we add to row_data[j+1][3], each with the same number of rows
            levels_climbed = min(levels_allowed, levels_needed)
            row_data[j+1][3] += levels_climbed * (j+1-min_undepleted_row)
            # count the degrees
            for k in range(min_undepleted_row,j+1):
                row_data[k][2] += levels_climbed
            # update the lowest priority shift not depleted
            while row_data[min_undepleted_row][2] >= deg + 1:
                min_undepleted_row += 1
                levels_allowed = deg + 1 - row_data[min_undepleted_row][2]
            # row_data update
            levels_needed -= levels_climbed
    
    # lookup
    
    # find lowest priority row out of range
    low = 0
    high = n-1
    mid = 0
    while low <= high:
        mid = (low+high)//2
        if i < row_data[mid][3]:
            high = mid - 1
        else:
            low = mid + 1
    l = low
    
    # find lowest priority row not depleted
    low = 0
    high = l-1
    mid = 0
    while low <= high:
        mid = (low+high)//2
        if row_data[mid][0] + deg + 1 < row_data[l-1][0]:
            low = mid + 1
        else:
            high = mid - 1
    
    min_undepleted_row = low
    
    # number of rows required to reach l-1
    threshold = row_data[l-1][3]
    # additional degrees
    d = 0
    while threshold + l - min_undepleted_row <= i:
        # add_d degrees available with l - min_undepleted_row rows
        add_d = (i - threshold) // (l - min_undepleted_row)
        add_d = min(add_d, deg + 1 + row_data[min_undepleted_row][0]  - row_data[l-1][0] - d)
        
        # add extra rows to count, and update additional degrees since threshold
        threshold += add_d * (l - min_undepleted_row)
        d += add_d
        # if the lowest priority row is exhausted, update range
        while row_data[min_undepleted_row][0] + deg + 1 <= row_data[l-1][0] + d:
            min_undepleted_row += 1
    
    # at this point, there is a partial priority block to add
    # priority here is determined by row number in E
    candidates = sorted(row_data[min_undepleted_row:l],key=lambda x: x[1])
    
    # select the correct row
    c = candidates[i - threshold][1]
    # degree is the difference in shifts with the last row, plus any additional degrees since the threshold
    d = row_data[l-1][0] - candidates[i - threshold][0] + d
    return c,d

def check_fast_perms():
    for i in range(100):
        n = random.randint(1,10)
        shift = [random.randint(0,n*n) for j in range(n)]
        
        priority = lambda c,d : shift[c] + d
        index = lambda i : (i%n,i//n)
        index_inv = lambda c,d : c + d*n
        
        priority_triplets = sorted([[priority(*index(i)),index(i),i] for i in range(n*(n+max(shift)-min(shift)+1))])
        priority_permutation = Permutation([t[2]+1 for t in priority_triplets])
        priority_permutation_inv = priority_permutation.inverse()

        # maps row c of self*J^d to position i in striped Krylov matrix
        # +/- 1 as permutations are 1-indexed, matrices are 0-indexed
        phi = lambda c,d : priority_permutation_inv(index_inv(c,d) + 1) - 1
        phi_inv = lambda i : index(priority_permutation(i + 1) - 1)
        
        sorting = [priority_permutation_inv(t+1)-1 for t in range(n*(n+max(shift)-min(shift)+1))]
        ind = [fast_perms(shift,n+max(shift)-min(shift),*index(i)) for i in range(n*(n+max(shift)-min(shift)+1))]
        if sorting != ind:
            print(shift,n+max(shift)-min(shift))
            break
        sorting = [phi_inv(i) for i in range(n*(n+max(shift)-min(shift)+1))]
        ind = [fast_perms_inv(shift,n+max(shift)-min(shift),i) for i in range(n*(n+max(shift)-min(shift)+1))]
        if sorting != ind:
            print(shift,n+max(shift)-min(shift))
            break