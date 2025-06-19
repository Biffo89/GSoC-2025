import math

class Permutation:
    def __init__(self,mapping):
        self.mapping = [0] * len(mapping)
        self.inverse = [0] * len(mapping)
        for i in range(len(mapping)):
            self.mapping[i] = mapping[i]
            self.inverse[mapping[i]] = i
    
    def apply(self,n):
        return self.mapping[n]
    
    def inv(self,n):
        return self.inverse[n]
    
    def matrix(self):
        mat_arr = [[0]*len(self.mapping) for i in range(len(self.mapping))]
        for i in range(len(self.mapping)):
            mat_arr[self.mapping[i]][i] = 1
        return matrix(mat_arr)

def row_degree(M):
    ## exists: M.row_degrees()  |  M.row_degrees(shifts)
    max_deg = [0] * M.nrows()
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if max_deg[i] < M[i][j].degree():
                max_deg[i] = M[i][j].degree()
    return matrix(max_deg,ncols=1)

def expand(P,delta):
    arrs = []
    for idx in range(P.nrows()):
        new_row = []
        for idx2 in range(P.ncols()*(delta+1)):
            if P[idx][idx2%P.ncols()].degree() < idx2//P.ncols():
                new_row.append(0)
            else:
                new_row.append(P[idx][idx2%P.ncols()].list()[idx2//P.ncols()])
        arrs.append(new_row)
    return matrix(arrs)

def compress(M,delta):
    m = M.ncols()//(delta+1)
    arrs = [[0]*m for idx in range(M.nrows())]
    for idx in range(M.ncols()):
        for idx2 in range(M.nrows()):
            arrs[idx2][idx%m] += M[idx2][idx]*x^(idx//m)
    return matrix(arrs)

def expand_shift(delta,E,J,s=None):
    if s is None:
        s = [0]*E.nrows()
    pre_block = []
    sigma = J.nrows()
    J_i = identity_matrix(sigma)
    for i in range(delta + 1):
        pre_block.append([E*J_i])
        J_i *= J
    unshifted = block_matrix(pre_block,subdivide=False)
    r = unshifted.rows()
    unsorted = [[idx//sigma + s[idx%sigma],idx%sigma,x] for idx,x in enumerate(r)]
    unsorted.sort(key=lambda x: [x[0],x[1]])
    return matrix([x[2] for x in unsorted])

def get_shift(s,delta):
    idx = [[s[i%len(s)]+i//len(s),i%len(s),i] for i in range(len(s)*(delta+1))]
    idx.sort()
    return Permutation([x[2] for x in idx])

def row_rank_profile(M):
    ## can be deduced from M.pivots()  ||  M.transpose().pivots()
    ## to add to sage?
    profile = []
    subm = matrix(0,M.ncols())
    r = rank(M)
    for i in range(M.nrows()):
        if subm.stack(matrix(M.rows()[i])).rank() > len(profile):
            profile.append(i)
            if len(profile) == r:
                break
            subm = subm.stack(matrix(M.rows()[i]))
    return r,profile

def col_rank_profile(M):
    ## can be deduced from M.pivots()  ||  M.transpose().pivots()
    ## to add to sage?
    profile = []
    subm = matrix(M.nrows(),0)
    r = rank(M)
    for i in range(M.ncols()):
        if subm.augment(matrix(M.nrows(),1,M.column(i))).rank() > len(profile):
            profile.append(i)
            if len(profile) == r:
                break
            subm = subm.augment(matrix(M.nrows(),1,M.column(i)))
    return r,profile

def krylov_rank_profile(E,J,s,delta):
    m = len(s)
    sigma = J.nrows()
    # phi_s(c,d) = phi_s[c+d*m]
    phi_s = get_shift(s,delta)
    r, c = row_rank_profile(E)
    d = [0] * r
    i = [phi_s.apply(c[idx]) for idx in range(r)]
    K_s = expand_shift(delta,E,J,s)
    M = block_matrix(r,1,[matrix(K_s.rows()[idx]) for idx in i],subdivide=False)
    for l in range(math.floor(math.log(delta,2))):
        if r > len(c):
            c = c + ([0]*(r-len(c)))
            d = d + ([0]*(r-len(d)))
        for idx in range(r):
            inv = phi_s.inv(i[idx])
            d[idx] = inv//m
            c[idx] = inv%r
        i2 = [phi_s.apply(c[idx] + m*(d[idx] + pow(2,l))) for idx in range(r)]
        k = i + i2
        M2 = block_matrix(r,1,[matrix(K_s.rows()[idx]) for idx in i2],subdivide=False)
        M = block_matrix(2,1,[[M],[M2]])
        r2, m2 = row_rank_profile(M)
        M = block_matrix(r2,1,[matrix(M.rows()[idx]) for idx in m2],subdivide=False)
        r = r2
        i = [k[idx] for idx in m2]
    _, j = col_rank_profile(M)
    return i,M,j

def linear_interpolation_basis(E,J,s,delta):
    m = len(s)
    sigma = E.ncols()
    phi_s = get_shift(s,delta)
    i,P_s,j = krylov_rank_profile(E,J,s,delta)
    r = len(i)
    c = []
    d = []
    for idx in range(len(i)):
        inv = phi_s.inv(i[idx])
        d.append(inv//m)
        c.append(inv%r)
    delt = [0]*m
    for idx in range(m):
        for idx2 in range(r):
            if c[idx2] == idx and d[idx2] >= delt[idx]:
                delt[idx] = 1 + d[idx2]
    K_s = expand_shift(delta,E,J,s)
    T_s = matrix(0,sigma)
    for idx in range(m):
        T_s = T_s.stack(matrix(1,sigma,K_s.rows()[phi_s.apply(idx + m*delt[idx])]))
    print('T_s:')
    print(T_s,'\n')
    C = matrix(0,sigma)
    D = matrix(0,sigma)
    for idx in range(r):
        C = C.stack(matrix(P_s.rows()[j[idx]]))
        D = D.stack(matrix(T_s.rows()[j[idx]]))
    R_s = D*C.inverse()
    print("[-R_s|I_m]:")
    unc = (-R_s).augment(identity_matrix(m))
    print(unc,'\n')
    return compress(unc,1)

if __name__ == "__main__":
    R.<x> = GF(97)[]
    print('row degree:')
    M = matrix([[x^2+36*x,31*x,0],[3*x+13,x+57,0],[96,96,1]])
    print(M.row_degrees(),'\n')
    print('expand:')
    exp = expand(M,2)
    print(exp,'\n')
    print('compress:')
    print(compress(exp,2),'\n')
    print('expand_shift(0):')
    E = matrix(R,[[27,49,29],[50,58,0],[77,10,29]])
    m = E.nrows()
    sigma = E.ncols()
    Z = matrix(R,[[0,1,0],[0,0,1],[0,0,0]])
    print(expand_shift(sigma,E,Z),'\n')
    print('expand_shift:')
    s = [3,0,2]
    delta = 4
    print(expand_shift(sigma,E,Z,s),'\n')
    print('row_rank_profile:')
    print(row_rank_profile(M))
    print('naive_krylov_rank_profile:')
    print(naive_krylov_rank_profile(E,Z,s,delta))
    print('krylov_rank_profile:')
    print(krylov_rank_profile(E,Z,s,delta),'\n')
    print('linear_interpolation_basis:')
    print(linear_interpolation_basis(E,Z,s,delta),'\n')