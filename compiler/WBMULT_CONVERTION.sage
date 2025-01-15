

def pk_mon(n):
    L = [1]
    # degree 1
    for i in range(n):
        L.append(PTI_MULT.gens()[i])
        
    # degree 2
    for i in range(n):
        for j in range(i+1,n):
             L.append(PTI_MULT.gens()[i]*PTI_MULT.gens()[j])

    return(L)
        

def mon3(n):
    L = [1]
    # degree 1
    for i in range(n):
        L.append(PTI_MULT.gens()[n+i])
        
    # degree 2
    for i in range(n):
        for j in range(i+1,n):
             L.append(PTI_MULT.gens()[n+i]*PTI_MULT.gens()[n+j])
    # degree 3
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                L.append(PTI_MULT.gens()[n+i]*PTI_MULT.gens()[n+j]*PTI_MULT.gens()[n+k])
    return(L)
        

# INPUT : a polynomial  in PTI_MULT en Y of degree 3
# OUTPUT : a list representing monomials in m bit integers
            
def poly_to_intlist(n,m,P,mon_list):
    P_set = P.set()
    length = sigma(n,3)
    quo = int(length//m)
    rem =  length%m
    # representing as int m bits
    L = []
    RES = []


    # First term is the constant
    L.append(int(PTI_MULT(P_set).constant_coefficient()))
    count = 1
    
    # degree 1
    for i in range(n):
        temp  = P_set.intersect(mon_list[count].set())
        if  temp.stable_hash() == mon_list[count].stable_hash():
            L.append(1)
        else:
            L.append(0)
        count += 1

                
    # degree 2
    if P.degree() > 1:
        for i in range(n):
            for j in range(i+1,n):
                temp  = P_set.intersect(mon_list[count].set())
                if  temp.stable_hash() == mon_list[count].stable_hash():
                    L.append(1)
                else:
                    L.append(0)
                count += 1
       
                        
    # degree 3
    if P.degree() > 2:
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    temp  = P_set.intersect(mon_list[count].set())
                    if  temp.stable_hash() == mon_list[count].stable_hash():
                        L.append(1)
                    else:
                        L.append(0)
                    count += 1
                    
                    
    # Building the list
    qquo = int(len(L)//m)
    rrem = len(L)%m
    for i in range(qquo):
        TEMP = [L[m*i+mm]<<mm for mm in range(m) ]
        res_temp = 0
        for ii in TEMP:
            if ii != 0 :
                res_temp = res_temp + ii
        RES.append(int(res_temp))
    # Adding last int
    TEMP = [L[qquo*m + mm]<<mm for mm in range(rrem) ]
    res_temp = 0
    for ii in TEMP:
        if ii != 0 :
            res_temp = res_temp + ii
    RES.append(int(res_temp))
    # Completing with 0
    FILL = quo + (rem != 0) - qquo - (rrem !=0)
    RES = RES + [0]*FILL
    


    return(RES)


# INPUT : a polynomial  in PTI_MULT en X of degree 2
# OUTPUT : a list representing monomials in m bit integers
            
def pk_to_intlist(n,m,P,mon_list):
    P_set = P.set()
    length = sigma(n,2)
    quo = int(length//m)
    rem =  length%m
    # representing as int m bits
    L = []
    RES = []


    # First term is the constant
    L.append(int(PTI_MULT(P_set).constant_coefficient()))
    count = 1
    
    # degree 1
    for i in range(n):
        temp  = P_set.intersect(mon_list[count].set())
        if  temp.stable_hash() == mon_list[count].stable_hash():
            L.append(1)
        else:
            L.append(0)
        count += 1

                
    # degree 2
    if P.degree() > 1:
        for i in range(n):
            for j in range(i+1,n):
                temp  = P_set.intersect(mon_list[count].set())
                if  temp.stable_hash() == mon_list[count].stable_hash():
                    L.append(1)
                else:
                    L.append(0)
                count += 1
       
                    
                    
    # Building the list
    qquo = int(len(L)//m)
    rrem = len(L)%m
    for i in range(qquo):
        TEMP = [L[m*i+mm]<<mm for mm in range(m) ]
        res_temp = 0
        for ii in TEMP:
            if ii != 0 :
                res_temp = res_temp + ii
        RES.append(int(res_temp))
    # Adding last int
    TEMP = [L[qquo*m + mm]<<mm for mm in range(rrem) ]
    res_temp = 0
    for ii in TEMP:
        if ii != 0 :
            res_temp = res_temp + ii
    RES.append(int(res_temp))
    # Completing with 0
    FILL = quo + (rem != 0) - qquo - (rrem !=0)
    RES = RES + [0]*FILL
    


    return(RES)

    
# INPUT: two lists of int
# OUTPUT : their XOR
def xor_intlist(A,B):
    L = []
    for i in range(len(A)):
        L.append(A[i].__xor__(B[i]))
    return(L)

    
# INPUT : a poly in int format, INmons the list of monomials in an input in int format
# OUTPUT : the evaluation of poly in INmons.
def eval_intlist_poly(poly,INmons):
    res = 0
    Xeval = []
    for i in range(len(poly)):
        Xeval.append(poly[i].__and__(INmons[i]))
    for i in Xeval:
        res = res.__xor__(i)
    res = sum(res.digits(2))%2
    return(res)

    
#Y_MON[3][0],Polist
#INtest = [choice([0,1]) for i in range(n)]
#Mons = mon_eval_int(n,8,INtest)
## eval
#res = 0
#Xeval = []
#for i in range(len(Polist)):
#    Xeval.append(Polist[i].__and__(Mons[i]))
#for i in Xeval:
#    res = res.__xor__(i)
#res = sum(res.digits(2))%2
# teval = substitute_variables(PTI_MULT, [PTI_MULT.gens()[i]for i in range(n)]+[PTI_MULT(i) for i in INtest], \
#                             Y_MON[3][0])
#res == teval



def simple_mon(IN):
    L = []
    for i in range(len(IN)):
        L.append(IN[i])
        for j in range(i+1,len(IN)):
            L.append(IN[i]*IN[j])
    return(L)

# INPUT : integer n, IN list of n bits
# OUTPUT : list of monomials of degree 2 in the bits of IN (in mon_order order), in the int m bits format

def mon_eval_int(n,m,IN):
    L = []
    OUT = []
    # Computing monomials
    for i in range(n):
        L.append(IN[i])
        for j in range(i+1,n):
            L.append(IN[i]*IN[j])
    # Converting into int of m bits
    length = n*(n+1)/2 + 1
    quo = int(length//m)
    rem =  length%m
    for i in range(quo):

            a = sum([ L[m*i+j] << j for j in range(m)] )
            OUT.append( a)
            #print(i,a,L[m*i : m*i +rem])
    
    if rem != 0:
        L.append(1)
        a =  sum([ L[m*quo+j] << j for j in range(rem)])
        OUT.append(a)
    else :
        OUT.append(1)
        
    return(OUT)

# test
#ntest = 6
#Ltest = []
#INtest = [choice([0,1]) for i in range(ntest)]
#for i in range(ntest):
#    Ltest.append(INtest[i])
#    for j in range(i+1,ntest):
#        Ltest.append(INtest[i]*INtest[j])
    

#PLOP,OUTtest = mon_eval_int(ntest,8,INtest)
#bl = []
#for i in OUTtest :
#    temp  = i.digits(2)
#    temp = temp + [0]*(8-len(temp))
#    bl  = bl + temp
    
    
#INtest, Ltest,  PLOP ,OUTtest, bl





# INPUT : a msg of n bits
# OUTPUT : a list representing monomials of deg 3 in msg bits in the form of m bit integers
            
def msg_to_monintlist(n,m,msg):
 
    L = [1]
    # degree 1
    for i in range(n):
        L.append(msg[i])
        
    # degree 2
    for i in range(n):
        for j in range(i+1,n):
             L.append(msg[i]*msg[j])
    # degree 3
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                L.append(msg[i]*msg[j]*msg[k])

                
    length = sigma(n,3)
    quo = int(length//m)
    rem =  length%m
    # representing as int m bits
    RES = []
                                        
    # Building the list
    qquo = int(len(L)//m)
    rrem = len(L)%m
    for i in range(qquo):
        TEMP = [L[m*i+mm]<<mm for mm in range(m) ]
        res_temp = 0
        for ii in TEMP:
            if ii != 0 :
                res_temp = res_temp + ii
        RES.append(int(res_temp))
    # Adding last int
    TEMP = [L[qquo*m + mm]<<mm for mm in range(rrem) ]
    res_temp = 0
    for ii in TEMP:
        if ii != 0 :
            res_temp = res_temp + ii
    RES.append(int(res_temp))
    # Completing with 0
    FILL = quo + (rem != 0) - qquo - (rrem !=0)
    RES = RES + [0]*FILL
    


    return(RES)
 
    
# INPUT : a msg of n bits
# OUTPUT : a list representing monomials of deg 2 in msg bits in the form of m bit integers
            
def sign_to_monintlist(n,m,sign):
 
    L = [1]
    # degree 1
    for i in range(n):
        L.append(sign[i])
        
    # degree 2
    for i in range(n):
        for j in range(i+1,n):
             L.append(sign[i]*sign[j])


                
    length = sigma(n,2)
    quo = int(length//m)
    rem =  length%m
    # representing as int m bits
    RES = []
                                        
    # Building the list
    qquo = int(len(L)//m)
    rrem = len(L)%m
    for i in range(qquo):
        TEMP = [L[m*i+mm]<<mm for mm in range(m) ]
        res_temp = 0
        for ii in TEMP:
            if ii != 0 :
                res_temp = res_temp + ii
        RES.append(int(res_temp))
    # Adding last int
    TEMP = [L[qquo*m + mm]<<mm for mm in range(rrem) ]
    res_temp = 0
    for ii in TEMP:
        if ii != 0 :
            res_temp = res_temp + ii
    RES.append(int(res_temp))
    # Completing with 0
    FILL = quo + (rem != 0) - qquo - (rrem !=0)
    RES = RES + [0]*FILL
    


    return(RES)
 