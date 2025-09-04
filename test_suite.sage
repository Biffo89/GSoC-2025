import random
import time
import tabulate

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
        if not(isinstance(degree,list)) or len(degree) != m or any(not(isinstance(elem,int)) or elem < 0 for elem in shift):
            raise ValueError('degree is not a non-negative list of the correct size.')
        # if not(isinstance(degree,int) and 1 <= degree and degree & (degree -1) == 0):
            # raise ValueError('sigma is not a integer that is a power of two.')

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
        maxdegree = 0
        if sigma == 0:
            maxdegree = 1
        elif sigma & (sigma-1) == 0:
            maxdegree = sigma
        else:
            maxdegree = sigma
            while maxdegree & (maxdegree+1) != 0:
                maxdegree |= maxdegree >> 1
            maxdegree += 1
        _, rows = self.E.krylov_basis(self.J,self.shift)
        self.degree = [int(0)] * self.m
        for x in rows:
            self.degree[x[0]] = x[1] + 1 if self.degree[x[0]] <= x[1] else self.degree[x[0]]
        self.degree = [random.randint(self.degree[i], maxdegree+1) for i in range(0, m)]

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
        
        # degree
        degree_string = f"[{','.join([f"int({i})" for i in self.degree])}]"
        
        return f"KrylovTestInstance({field_string},int({self.m}),int({self.sigma}),options={{'manual':{{'E':{E_string},'J':{J_string},'shift':{shift_string},'degree':{degree_string}}}}})"

    #OBSOLETE
    def naive_krylov_profile(self):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.naive_krylov_basis(self.J,self.shift,self.degree)

    def krylov_profile(self, output_rows=True, algorithm=None):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.krylov_basis(self.J, shifts=self.shift, degrees=self.degree, output_rows=output_rows, algorithm=algorithm)

    #OBSOLETE
    def krylov_profile_early_exit(self):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.krylov_profile_early_exit(self.J,self.shift,self.degree)

    def linear_interpolation_basis(self,var=None):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.krylov_kernel_basis(self.J,degrees=self.degree,shifts=self.shift,var=var)

    #OBSOLETE
    def linear_interpolation_basis_fast_perm(self):
        self.E._clear_cache()
        self.J._clear_cache()
        return self.E.linear_interpolation_basis_fast_perm(self.J,'x',self.shift,self.degree)

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
        profile = test.krylov_profile(algorithm='elimination')
    return profile == test.krylov_profile(algorithm='naive')

def is_basis_correct(test,poly_basis=None,coeff_basis=None,profile=None):
    try:
        if profile is None:
            profile = test.krylov_profile()
        if poly_basis is None:
            poly_basis = test.linear_interpolation_basis(var='x')
        if coeff_basis is None:
            coeff_basis = test.linear_interpolation_basis(var=None)
        if not poly_basis[0].is_popov(shifts=test.shift):
            return False
        if poly_basis[0].nrows()*poly_basis[0].ncols() != 0 and poly_basis[0].degree() > test.sigma:
            return False
        # check equality
        coords = test.E._krylov_row_coordinates(test.shift, test.degree)
        poly_basis_ex = matrix.block([[poly_basis[0].coefficient_matrix(i).matrix_from_columns([j for j in range(test.m) if i <= test.degree[j]]) for i in range(max(test.degree,default=0) + 1)]],subdivide=False).with_permuted_columns(Permutation([x[2] + 1 for x in coords]))
        coeff_basis_ex = matrix(test.field, test.m, sum(test.degree)+test.m)
        for i in range(len(coeff_basis[1])):
            idx = coeff_basis[1][i][2]
            for j in range(test.m):
                coeff_basis_ex[j,idx] = coeff_basis[0][j,i]
        if poly_basis_ex != coeff_basis_ex:
            return False
        if poly_basis[1] != coeff_basis[1]:
            return False
        for i,x in enumerate(coeff_basis[1]):
            if x[:2] != coords[x[2]][:2]:
                return False
        # check kernel
        if coeff_basis[0]*test.E.krylov_matrix(test.J,test.shift,test.degree).matrix_from_rows([x[2] for x in coeff_basis[1]]) != 0:
            return False
        if len(profile[1]) != sum((poly_basis[0][i,i].degree() for i in range(poly_basis[0].nrows()))):
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
        # profile terminates
        try:
            profile = test.krylov_profile(algorithm='elimination')
            success[i]['profile_completed'] = True
        except:
            success[i]['profile_completed'] = False
        if not success[i]['profile_completed']:
            continue
        # profile is correct
        success[i]['profile_correct'] = is_profile_correct(test, profile)
        # basis terminates
        try:
            coeff_basis = test.linear_interpolation_basis(var=None)
            poly_basis = test.linear_interpolation_basis(var='x')
            success[i]['basis_completed'] = True
        except:
            success[i]['basis_completed'] = False
        if not success[i]['basis_completed']:
            continue
        # basis is correct
        success[i]['basis_correct'] = is_basis_correct(test,poly_basis,coeff_basis,profile)
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

#OBSOLETE
def measure_once(b):
    test = KrylovTestInstance(GF(97),int(1024),int(1024),{'E_mode':None,'J_mode':None,'shift_mode':'random'})
    print(test.E.base_ring())
    print(timeit('test.J * test.J',globals={'test':test}))
    timer = []
    #test.E.krylov_profile(test.J,test.degree,test.shift,False,timer)
    print(timer)
    timer = []
    test.E.linear_interpolation_basis(test.J,test.degree,'x',test.shift,b,timer)
    print(timer)

#OBSOLETE
def compare_naive():
    cases = [(GF(2),int(4096),int(1024),int(4096),int(512)),(GF(3**5),int(64),int(64),int(64),int(64)),(GF(2**8),int(1024),int(512),int(512),int(256)),(GF(257),int(1024),int(512),int(1024),int(512)),(GF(3**10),int(64),int(64),int(64),int(64)),(GF(2**16),int(64),int(64),int(64),int(64)),(GF(65537),int(1024),int(512),int(1024),int(256)),(GF(next_prime(2**26)),int(256),int(256),int(256),int(256))]

    for case in cases:
        field = case[0]
        print(f"GF({field.order()})\n")
        mul_results = {}
        for sig in sorted(list(set([case[i+1]//int(j) for i in range(4) for j in [1,2,4]]))):
            M1 = matrix.random(field, sig)
            M2 = matrix.random(field, sig)
            t = timeit('M1*M2', globals={'M1':M1, 'M2':M2})
            mul_results[sig] = f"{t.stats[3]:.3g} {t.stats[4]}"
        data = [[key,mul_results[key]] for key in mul_results.keys()]
        print('matrix multiplication')
        print(tabulate.tabulate(data,headers=['n','time'],tablefmt='grid'))
        print()
        for func in ['profile']:
            results_u = {}
            results_h = {}
            columns = set()
            rows = set()
            for m_size in ['small']:
                i = 2 * (func == 'basis') + 1 * (m_size == 'big') + 1
                m_range = [int(1),int(2),int(3),int(4),int(5)] if field != GF(243) else [int(1),int(2),int(3),int(4),int(5),int(8),int(9),int(16),int(17),int(32),int(33)]
                for m in m_range:
                    rows.add(m)
                    if m not in results_u:
                        results_u[m] = {}
                        results_h[m] = {}
                    for sig in [int(1 << i) for i in range(case[i].bit_length()-1)]:
                        if m != sig and (m != 16 or sig != case[i]):
                            continue
                        columns.add(sig)
                        test_u = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'uniform'})
                        test_h = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'hermite'})
                        t_nu = timeit("test_u.krylov_profile(algorithm='naive')" if func == 'profile' else 'test_u.linear_interpolation_basis(output_coefficients=False)',globals={'test_u':test_u},number=125)
                        t_nh = timeit("test_h.krylov_profile(algorithm='naive')" if func == 'profile' else 'test_h.linear_interpolation_basis(output_coefficients=False)',globals={'test_h':test_h},number=125)
                        t_u = timeit("test_u.krylov_profile(algorithm='elimination')" if func == 'profile' else 'test_u.linear_interpolation_basis(output_coefficients=True)',globals={'test_u':test_u},number=125)
                        t_h = timeit("test_h.krylov_profile(algorithm='elimination')" if func == 'profile' else 'test_h.linear_interpolation_basis(output_coefficients=True)',globals={'test_h':test_h},number=125)
                        results_u[m][sig] = f"{t_u.stats[3]:.3g} {t_u.stats[4]} ({t_u.stats[0]})" if t_nu is None else f"{t_nu.stats[3]:.3g} {t_nu.stats[4]} ({t_nu.stats[0]}) / {t_u.stats[3]:.3g} {t_u.stats[4]} ({t_u.stats[0]})"
                        results_h[m][sig] = f"{t_h.stats[3]:.3g} {t_h.stats[4]} ({t_h.stats[0]})" if t_nh is None else f"{t_nh.stats[3]:.3g} {t_nh.stats[4]} ({t_nh.stats[0]}) / {t_u.stats[3]:.3g} {t_u.stats[4]} ({t_h.stats[0]})"
            data_u = [[results_u[m].get(sigma,None) for sigma in sorted(list(columns))] for m in sorted(list(rows))]
            data_h = [[results_h[m].get(sigma,None) for sigma in sorted(list(columns))] for m in sorted(list(rows))]
            print('uniform '+func)
            print(tabulate.tabulate(data_u,headers=['*   sigma\n    *    \nm       *'] + sorted(list(columns)), showindex = sorted(list(rows)),tablefmt='grid'))
            print()
            print('hermite '+func)
            print(tabulate.tabulate(data_h,headers=['*   sigma\n    *    \nm       *'] + sorted(list(columns)), showindex = sorted(list(rows)),tablefmt='grid'))
            print()

def benchmark():
    cases = [(GF(2),int(2048),int(256),int(2048),int(256)),(GF(3**5),int(64),int(64),int(64),int(64)),(GF(2**8),int(512),int(256),int(512),int(256)),(GF(257),int(512),int(256),int(512),int(256)),(GF(3**10),int(32),int(32),int(32),int(32)),(GF(2**16),int(32),int(32),int(32),int(32)),(GF(65537),int(512),int(256),int(512),int(256)),(GF(next_prime(2**26)),int(128),int(128),int(128),int(128))]
    #cases = [(x[0],x[1]//int(4),x[2]//int(4),x[3]//int(4),x[4]//int(4)) for x in cases]
    for case in cases:
        field = case[0]
        print(f"GF({field.order()})\n")
        mul_results = {}
        gau_results = {}
        for sig in sorted(list(set([case[i+1]//int(j) for i in range(4) for j in [1,2,4]]))):
            M1 = matrix.random(field, sig)
            M2 = matrix.random(field, sig)
            t_mul = timeit('M1._clear_cache();M2._clear_cache();M1*M2', globals={'M1':M1, 'M2':M2})
            t_gau = timeit('M1._clear_cache(); M1.echelon_form()', globals={'M1':M1})
            mul_results[sig] = f"{t_mul.stats[3]:.3g} {t_mul.stats[4]}"
            gau_results[sig] = f"{t_gau.stats[3]:.3g} {t_gau.stats[4]}"
        data = [[key,mul_results[key],gau_results[key]] for key in mul_results.keys()]
        print(tabulate.tabulate(data,headers=['n','multiply','gauss'],tablefmt='grid'))
        print()
        for func in ['basis','profile']:
            results_u = {}
            results_h = {}
            columns = set()
            rows = set()
            for m_size in ['small','big']:
                i = 2 * (func == 'basis') + 1 * (m_size == 'big') + 1
                m_range = [int(1),int(4),int(16)] if m_size == 'small' else [case[i]//int(4), case[i]//int(2), case[i]]
                for m in m_range:
                    rows.add(m)
                    if m not in results_u:
                        results_u[m] = {}
                        results_h[m] = {}
                    for sig in [case[i]//int(4),case[i]//int(2),case[i]]:
                        #if sig not in [m,case[i]] or m not in [sig,int(16),case[i]]:
                            #continue
                        columns.add(sig)
                        test_u = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'uniform'})
                        test_h = KrylovTestInstance(field,m,sig,{'E_mode':None,'J_mode':None,'shift_mode':'hermite'})
                        t_nu = None
                        t_nh = None
                        #t_nu = timeit("test_u.krylov_profile(algorithm='naive')" if func == 'profile' else 'test_u.linear_interpolation_basis(output_coefficients=False)',globals={'test_u':test_u})
                        #t_nh = timeit("test_h.krylov_profile(algorithm='naive',var=''x'')" if func == 'profile' else 'test_h.linear_interpolation_basis(var=''x'')',globals={'test_h':test_h})
                        t_u = timeit("test_u.krylov_profile(algorithm='elimination')" if func == 'profile' else 'test_u.linear_interpolation_basis()',globals={'test_u':test_u})
                        t_h = timeit("test_h.krylov_profile(algorithm='elimination')" if func == 'profile' else 'test_h.linear_interpolation_basis()',globals={'test_h':test_h})
                        results_u[m][sig] = f"{t_u.stats[3]:.3g} {t_u.stats[4]}" if t_nu is None else f"{t_u.stats[3]:.3g} {t_u.stats[4]} / {t_nu.stats[3]:.3g} {t_nu.stats[4]}"
                        results_h[m][sig] = f"{t_h.stats[3]:.3g} {t_h.stats[4]}" if t_nh is None else f"{t_h.stats[3]:.3g} {t_h.stats[4]} / {t_nh.stats[3]:.3g} {t_nh.stats[4]}"
            data_u = [[results_u[m].get(sigma,None) for sigma in sorted(list(columns))] for m in sorted(list(rows))]
            data_h = [[results_h[m].get(sigma,None) for sigma in sorted(list(columns))] for m in sorted(list(rows))]
            print('uniform '+func)
            print(tabulate.tabulate(data_u,headers=['*   sigma\n    *    \nm       *'] + sorted(list(columns)), showindex = sorted(list(rows)),tablefmt='grid'))
            print()
            print('hermite '+func)
            print(tabulate.tabulate(data_h,headers=['*   sigma\n    *    \nm       *'] + sorted(list(columns)), showindex = sorted(list(rows)),tablefmt='grid'))
            print()
