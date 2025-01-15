"""    
###############################################################################################################################

Contains the functions to compute an affine multiple of a Polynomial P

Notes : 

    - For smart_mon_red, the ring Ring has to be of the form F(y)[x] 
        For instance :
    
            F = GF(2**10)
            PR = F['y']
            y = PR.gens()[0]
            FF = Frac(PR)
            RRing = PolynomialRing(FF,'x')
            
        To have access to PR after for HW computation for instance
        
        
    - Badly optimized, only try to use for small degrees. (~12-16 max)


###############################################################################################################################
"""

# reduces a bit more efficiently a monomial x^n mod P
# INPUT : integer n, polynomial P in x and y, polynomial ring Ring containing variables x and y
# OUTPUT : reduction of x^n mod P in Ring

def smart_mon_red(n,P,Ring):
    d= P.degree(x)
    I = Ring.ideal(P)
    a = n//d
    b = n%d
    MON = I.reduce(x**a)
    MON2 = I.reduce(x**a)
    L = bin(d)[2:]
    k = len(L)
    for i in range (1,k):
        MON2 = I.reduce(MON2*MON2)
        if L[i] == '1':
            MON2 = I.reduce(MON2*MON)
    MON2 = MON2 * I.reduce(x**b)
    
    return(I.reduce(MON2))


    # Reduces the monomials x^i for i<d where d= deg(P) and 
# computes a matrix M containing the reduction in the base (1,y,....,y^(d-1))
# INPUT: Polynomial P in x, Polynomial Ring P in x and y
# OUTPUT : The matrix M of reductions.

def mons_reduction(P,Ring):
    L = []
    d = P.degree(x)
    for i in range(d):
        L.append(smart_mon_red(2**i,P,Ring))
    M = [[0]*(d-1)+[1]]
    for i in range(d):
        dd = L[i].degree(x)
        Lc = L[i].coefficients(sparse = false)
        Lc.reverse()
        M.append([0]*(d-dd-1)+Lc) 
        #if i == d-1:
            #print(L[i])
    return(M)


    # Computes the Hamming Weight of integer n
def hw(n):
    L = [int(i) for i in bin(n)[2:]]
    return(sum(L))



# Computes max of the Hamming Weight of the degree of the polynomials in y
# in a list of polynomials.
# INPUT: List L of polynomials in y, Field F of their coefficients.
# OUTPUT : The max of the Hamming Weight of the degree of the polynomials in y.
def hw_check(L,F):
    c = 0
    Ring = F['y']
    for i in L:
        Lmon = Ring(i).monomials()
        for j in Lmon :
            if j != 0 and j != 1:
                c = (c >= hw(j.degree()))*c + (c < hw(j.degree()))*hw(j.degree())
    return(c)



    # On the input of the output matrix of mons_reduction, outputs the coefficients of the affine multiple
# INPUT: Matrix M_red output of mons_reduction.
# OUTPUT :The list LL of the coefficients of the polynomial.
def affine_coeffs(M_red):
    
    AA = Matrix(M_red)
    LL = kernel(AA).basis()[0]
    den = LL[-1]
    c= 0
    cd = 0
    for i in LL:
        if c != 0:
            cd = c.degree()
        c = (cd >= i.denominator().degree())*c + (cd < i.denominator().degree())*i.denominator()
    LL = [c*i for i in LL]
        
    return(LL)



    ##########
    ## TEST ##
    ##########


"""
F = GF(2**10)#['a','b','c','d']
#[a,b,c,d] = list(F.gens())
PR = F['y']
y = PR.gens()[0]
FF = Frac(PR)
RRing = PolynomialRing(FF,'x')
x = RRing.gens()[0]

a = F.random_element()
b = F.random_element()
c = F.random_element()
d = F.random_element()

P1 = y + x**3 + d*x**2+  a
P2 = y + x**3 + d*x**2 + b
P3 = y + x**3 + d*x**2 + c

MOD = P1*P2*P3
FACTOR =  x**3 + d*x**2 + a


RES  = mons_reduction(MOD ,RRing)
AFF_COEFF = affine_coeffs(RES)
hw_check(AFF_COEFF)


### CHECK IF TRUE ###
sum([PR(AFF_COEF[i])(FACTOR)* RRing.gens()[0]**(2**(i-1)) for i in range(1,len(AFF_COEFF))])+ PR(AFF_COEFF[0])(FACTOR) == 0

"""