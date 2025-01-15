

# INPUT : integer n, integer 0<d <=n
# OUPUT : Number of monomials in n variables of at most degree d.   


def sigma(n,d):
    res =  0 
    for i in range(d+1):
        res = res + binomial(n,i)
    return(res)


# Input : integer n , p_coord  a parsed poly  from BPR(n), list of n poly from BPR(n) and BPR.
# Output : out is the composition of poly_list with p_coord


def compo_BPR_set(p_coord,poly_list,BPR):
    out = Monomial(BPR).set()
    for m in p_coord:
        temp  = Monomial(BPR).set()
        for i in m :
            temp = temp.cartesian_product(poly_list[i])
        out = (out.union(temp)).diff(out.intersect(temp))
    out  = (out.union(Monomial(BPR).set())).diff(out.intersect(Monomial(BPR).set()))
    return(out)

## TEST compo_PBR_set ##
#BPR_test = BooleanPolynomialRing(5,'x')
#M_test = BPR_test.random_element()
#M_parsed = parse_BPR(5, M_test,BPR_test)
#IN_test = []
#for i in range(5):
#    IN_test.append((BPR_test.random_element()).set())
#OUT_test = compo_BPR_set(M_parsed,IN_test,BPR_test)

# Input : integer n , p_coord  a coord from a parsed element from GL(n), list of n poly from BPR(n) and BPR.
# Output : out is the composition of poly_list with p_poly

def compo_GL(p_coord,poly_list):
    out = 0
    for i in p_coord:
        out = out + poly_list[i]
    return(out)


## TEST compo_GL ##
#BPR_test = BooleanPolynomialRing(5,'x')
#M_test = G_test.random_element()
#M_parsed = parse_GL(5, M_test.list())
#IN_test = []
#for i in range(5):
#    IN_test.append(BPR_test.random_element())

#OUT_test = compo_GL(M_parsed[0],IN_test)v
#IN_test,M_parsed[0]

# Input : integer n , polynomial  poly in BPR(n).
# Output : out[i] is the list of indexes of the variables in the i-th monomial of poly.
def parse_BPR(n,poly,BPR):
    var = list(BPR.gens())
    mons = poly.monomials()
    out = []
    for m in mons:
        temp = []
        for i in list(m.variables()):
            temp.append(var.index(i))
        out.append(temp)
    return(out)

## TEST parse_BPR ##
#BPR_test = BooleanPolynomialRing(20,'x')
#P_test = BPR_test.random_element()
#P_test
#parse_BPR(20,P_test,BPR_test)

# Input : integer n ,  matrix m in GL(F2,n).
# Output : out[i] is the list of indexes where mat[i] is 1.
def parse_GL(n,m):
    mat = m.list()
    out = []
    for i in range(n):
        temp = []
        for j in range(n):
            if mat[i][j] == 1 :
                temp.append(j)
        out.append(temp)
    return(out)


# Input : integer n , IN list of n input variables, matrix m in GL(F2,n), a vector v in GF(2)^n.
# Output : polynomials in PTI computing the coordinates of m +v

def AFF_to_PTI(n,IN, m,v):
    mat = m.list()
    out = []
    for i in range(n):
        temp = 0
        for j in range(n):
            if mat[i][j] == 1 :
                temp = temp + IN[j]
        out.append(temp + v[i])
    return(out)



##########################################
##########################################
##########################################

def conv_VS_to_GF(v,Field):
    s = '0b'
    for i in range(len(v)) :
        s = s + str(v[len(v)-i-1])
    a = int(s,2)
    return(Field.fetch_int(a))

def conv_GF_to_VS(n,e,Field):
    s = bin(e.integer_representation())
    v = [int(i) for i in s[2:]]
    v.reverse()
    v =  v  + [0]*(n-len(v)) 
    return(v)
    
#F = GF(2**20)
#a = F.random_element()
#a ==  conv_VS_to_GF(conv_GF_to_VS(20,a,F),F)


def eval_poly_set(IN,P,BPR):
    res = P
    for i in range(len(IN)):
        enc = (BPR.gens()[i] + IN[i]).set()
        res = ll_red_nf_redsb(res, enc)
    return(res)




########################################
########################################
########################################

#Bad conversion, but it works
def aff_list_to_matrix(n,L,PR):
    LL = []
    Vect = []
    for i in L:
        i_copy = i
        temp = []
        for j in range(n):
            m = i_copy.set().multiples_of(PTI_MULT.gens()[j].lm())
            # evaluating these multiples at 1 in j
            m = m.subset1(j)
            temp.append(GF(2)(PTI_MULT(m)))
            #evaluate i_copy at 0 on j
            enc  = ll_encode(PTI_MULT.gens()[j].set())
            i_copy = ll_red_nf_redsb(i_copy, enc ).set()
        # constant term
        temp.append(GF(2)(PTI_MULT(i_copy)))
        LL.append(temp)
    return(Matrix(LL))








    ####################################
    ####################################
    ####################################

def compute_monomial_PTI_GROS_ctrl_size(n,d,tau,MOD,PTI,GROS):
    # Just compute AND over F2^n trough a normal basis.
    # Gives the vector space with respect to the normal basis of the modulus : V = K.vector_space   
    P0 = 0
    for i in range(n):
        P0 = PTI.gens()[i+tau]* (GROS.gens()[0]^i) + P0

    MOD_inj = MOD(GROS.gens()[0])
    
    # Perform Square and Multiply
    RES =  1
    TEMP = P0
    d_bits = [int(i) for i in list(bin(d))[2:]]
    d_bits.reverse
    for i in d_bits:
        if i == 1:
            RES = RES*TEMP
            RES = RES%MOD_inj
            TEMP = TEMP*TEMP
            TEMP = TEMP%MOD_inj
        else:
            TEMP = TEMP*TEMP
            TEMP = TEMP%MOD_inj
        RES = RES%MOD_inj
    return(RES)




def compute_monomial_PTI_GROS(n,d,tau,MOD,PTI,GROS):
    # Just compute AND over F2^n trough a normal basis.
    # Gives the vector space with respect to the normal basis of the modulus : V = K.vector_space   
    P0 = 0
    for i in range(n):
        P0 = PTI.gens()[i+tau]* (GROS.gens()[0]^i) + P0

    MOD_inj = MOD(GROS.gens()[0])
    
    # Perform Square and Multiply
    RES =  1
    d_bits = [int(i) for i in list(bin(d))[2:]]
    for i in d_bits:
        if i == 1:
            RES = RES*RES
            RES = RES * P0
        else:
            RES = RES*RES
        RES = RES%MOD_inj
    return(RES) 



# Computes monomials with specific input in PTIlist whose elements are in PTI
def compute_monomial_PTIlist_GROS_ctrl_size(n,d,MOD,PTIlist,GROS):
    # Just compute AND over F2^n trough a normal basis.
    # Gives the vector space with respect to the normal basis of the modulus : V = K.vector_space   
    P0 = 0
    for i in range(n):
        P0 = PTIlist[i]* (GROS.gens()[0]^i) + P0

    MOD_inj = MOD(GROS.gens()[0])
    
    # Perform Square and Multiply
    RES =  1
    TEMP = P0
    d_bits = [int(i) for i in list(bin(d))[2:]]
    d_bits.reverse()
    for i in d_bits:
        if i == 1:
            RES = RES*TEMP
            RES = RES%MOD_inj
            TEMP = TEMP*TEMP
            TEMP = TEMP%MOD_inj
        else:
            TEMP = TEMP*TEMP
            TEMP = TEMP%MOD_inj
    return(RES) 






# Computes monomials with specific input in PTIlist whose elements are in PTI
def compute_monomial_PTIlist_GROS(n,d,MOD,PTIlist,GROS):
    # Just compute AND over F2^n trough a normal basis.
    # Gives the vector space with respect to the normal basis of the modulus : V = K.vector_space   
    P0 = 0
    for i in range(n):
        P0 = PTIlist[i]* (GROS.gens()[0]^i) + P0

    MOD_inj = MOD(GROS.gens()[0])
    
    # Perform Square and Multiply
    RES =  1
    d_bits = [int(i) for i in list(bin(d))[2:]]
    for i in d_bits:
        if i == 1:
            RES = RES*RES
            RES = RES * P0
        else:
            RES = RES*RES
        RES = RES%MOD_inj
    return(RES)



## Not optimized at all, made to be run once or twice. 
# Input : integer n
# Output : The coordinates of the multiplication in GF(2**n) in the vector space GF(2)**n 
def poly_prod_GF(n,MOD_inj,PTI,GRO):

    # Just compute AND over F2^n trough a normal basis.
    # Gives the vector space with respect to the normal basis of the modulus : V = K.vector_space()

    
    P0,P1 = 0,0
    for i in range(n):
        P0 = PTI.gens()[i]* (GRO.gens()[0]^i) + P0
        P1 = PTI.gens()[i+n]* (GRO.gens()[0]^i) + P1
  
    
    prod = (P0*P1)%MOD_inj 
        
    return(prod.coefficients(sparse = false))


# INPUT : integer n , prod_list are the coordinates of the multiplcation in GF(2*n)
# OUTPUT : XY_list,  XY_list[i][j] is the variable that is multiple of x_j in the ith coordinate of prod_list

# can maybe optimize the 'getting variables part'
def prod_GF_to_XYlist(n,prod_list):
    XY_list = []
    # going through coordinates
    for i in range(n):
        line = []
        for j in range(n):
            # getting the multiples of variable j
            m = prod_list[i].set().multiples_of(PTI_MULT.gens()[j].lm())
            # evaluating these multiples at 1 in j
            m = m.subset1(j)
            # getting the variables left
            var_list = list(m.set().vars().iterindex())
            if var_list == []:
                line.append([])
            else:
                line.append(var_list)
        # getting the constant coefficient by evaluating to 0
        dictio = {}
        for k in range(n):
            dictio[PTI_MULT.gens()[k]] = 0
        temp = prod_list[i].subs(dictio)
        #### The last term is the constant coefficient ####
        line.append(list(temp.set().vars().iterindex())+[temp.constant_coefficient()])
        XY_list.append(line)
    return(XY_list)
            
