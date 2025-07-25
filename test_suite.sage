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
            if sigma == 0:
                self.shift = [0]*m
            else:
                full_range = [i for i in range(m*sigma)]
                shuffle(full_range)
                self.shift = full_range[:m]
        if options['shift_mode'] == 'hermite':
            self.shift = [i*sigma for i in range(m)]
        if options['shift_mode'] == 'reverse_hermite':
            self.shift = [i*sigma for i in reversed(range(m))]
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
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.naive_krylov_rank_profile(self.J,self.degree,self.shift)

    def krylov_rank_profile(self):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.krylov_rank_profile(self.J,self.degree,self.shift)

    def krylov_rank_profile_early_exit(self):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.krylov_rank_profile_early_exit(self.J,self.degree,self.shift)

    def linear_interpolation_basis(self,polynomial_output=True):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.linear_interpolation_basis(self.J,self.degree,'x',self.shift, polynomial_output=polynomial_output)

    def linear_interpolation_basis_fast_perm(self):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.linear_interpolation_basis_fast_perm(self.J,self.degree,'x',self.shift)

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
                            if random.random() < 0.01: # ~2 minutes
                                test_set.append(KrylovTestInstance(f,m,sig,{'E_mode':E_opt,'J_mode':J_opt,'shift_mode':shift}))

    return test_set

def is_profile_correct(test,profile=None):
    if profile is None:
        profile = test.krylov_rank_profile()
    return profile == test.naive_krylov_rank_profile()

def is_basis_correct(test,basis=None,profile=None):
    if profile is None:
        profile = test.krylov_rank_profile()
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
        if basis.expansion(test.degree).with_permuted_columns(priority_permutation)*test.E.striped_krylov_matrix(test.J,test.degree,test.shift) != matrix.zero(test.field,test.m,test.sigma):
            return False
    except:
        return False
    if len(profile[0]) != sum((basis[i,i].degree() for i in range(basis.nrows()))):
        return False
    return True

def run_test_set(tests):
    success = [{'profile_completed':False,'basis_completed':False,'profile_correct':False,'basis_correct':False} for i in range(len(tests))]
    start = time.time()
    for i in range(len(tests)):
        test = tests[i]
        profile = None
        basis = None
        attempts = 1
        # profile terminates
        try:
            profile = test.krylov_rank_profile()
            success[i]['profile_completed'] = True
        except:
            success[i]['profile_completed'] = False
        if not success[i]['profile_completed']:
            continue
        # profile is correct
        success[i]['profile_correct'] = is_profile_correct(test,profile)
        # basis terminates
        try:
            basis = test.linear_interpolation_basis()
            success[i]['basis_completed'] = True
        except:
            success[i]['basis_completed'] = False
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
        if failure:
            print(i)
            print(success[i])
            print(tests[i].generator())

def measure_once(b):
    test = KrylovTestInstance(GF(97),int(1024),int(1024),{'E_mode':None,'J_mode':None,'shift_mode':'random'})
    print(test.E.base_ring())
    print(timeit('test.J * test.J',globals={'test':test}))
    timer = []
    #test.E.krylov_rank_profile(test.J,test.degree,test.shift,False,timer)
    print(timer)
    timer = []
    test.E.linear_interpolation_basis(test.J,test.degree,'x',test.shift,b,timer)
    print(timer)

def benchmark():
    cases = [(GF(2),int(4096),int(1024),int(4096),int(512)),(GF(3**5),int(64),int(64),int(64),int(64)),(GF(2**8),int(512),int(1024),int(512),int(256)),(GF(257),int(1024),int(512),int(1024),int(512)),(GF(3**10),int(64),int(64),int(64),int(64)),(GF(2**16),int(64),int(64),int(64),int(64)),(GF(65537),int(1024),int(512),int(1024),int(256)),(GF(next_prime(2**26)),int(256),int(256),int(256),int(256))]

    for case in cases:
        field = case[0]
        print (f"GF({field.order()})\n")
        m_range = [int(1),int(4),int(16)]
        for sig in [case[1]//4,case[1]//2,case[1]]:
            print("")
        for m in m_range:
            for sig in [case[1]//int(4),case[1]//int(2),case[1]]:
                print(f"m = {m}, sigma = {sig}, shift = uniform")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'uniform'})

                print(f"krylov_rank_profile: {timeit('test.krylov_rank_profile()',globals={'test':test})}s")

                print(f"m = {m}, sigma = {sig}, shift = hermite")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'hermite'})

                print(f"krylov_rank_profile: {timeit('test.krylov_rank_profile()',globals={'test':test})}s")
            print("\n")
        m_range = [case[2]//int(2),int(15)*(case[2]//int(16)),case[2]]
        for m in m_range:
            for sig in [case[2]//int(4),case[2]//int(2),case[2]]:
                print(f"m = {m}, sigma = {sig}, shift = uniform")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'uniform'})

                print(f"krylov_rank_profile: {timeit('test.krylov_rank_profile()',globals={'test':test})}s")

                print(f"m = {m}, sigma = {sig}, shift = hermite")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'hermite'})

                print(f"krylov_rank_profile: {timeit('test.krylov_rank_profile()',globals={'test':test})}s")
            print("\n")
        m_range = [int(1),int(4),int(16)]
        for m in m_range:
            for sig in [case[3]//int(4),case[3]//int(2),case[3]]:
                print(f"m = {m}, sigma = {sig}, shift = uniform")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'uniform'})

                print(f"linear_interpolation_basis (polynomial): {timeit('test.linear_interpolation_basis(True)',globals={'test':test})}s")

                print(f"m = {m}, sigma = {sig}, shift = hermite")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'hermite'})

                print(f"linear_interpolation_basis (polynomial): {timeit('test.linear_interpolation_basis(True)',globals={'test':test})}s")
            print("\n")
        m_range = [case[4]//int(2),int(15)*(case[4]//int(16)),case[4]]
        for m in m_range:
            for sig in [case[4]//int(4),case[4]//int(2),case[4]]:
                print(f"m = {m}, sigma = {sig}, shift = uniform")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'uniform'})

                print(f"linear_interpolation_basis (polynomial): {timeit('test.linear_interpolation_basis(True)',globals={'test':test})}s")

                print(f"m = {m}, sigma = {sig}, shift = hermite")
                test = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'hermite'})

                print(f"linear_interpolation_basis (polynomial): {timeit('test.linear_interpolation_basis(True)',globals={'test':test})}s")
            print("\n")
