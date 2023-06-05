



# INPUT : integer n, integer 0<d <=n
# OUPUT : Number of monomials in n variables of at most degree d.   


def sigma(n,d):
    res =  0 
    for i in range(d+1):
        res = res + binomial(n,i)
    return(res)


def file_to_matrix(n,file_name):
    MAT = []
    LINE = []
    File = open(file_name, "rb")
    poly = []
    for i in range(n):
        LINE = []
        for j in range(n+1):
            # Filing lines
            entier = int((sigma(n,3))//8)+((sigma(n,3))%8 != 0) 
            LINE.append(list(File.read(entier)))
        MAT.append(LINE)
    File.close()
    return(MAT)

def file_to_pk(n,a,file_name):
    PK = []
    File = open(file_name, "rb")
    poly = []
    for i in range(n-a):
        entier = int((sigma(n,2))//8)+((sigma(n,2))%8 != 0)
        PK.append(list(File.read(entier)))
    File.close()
    return(PK)


def file_to_test(ntest,file_name):
    TEST_msg = []
    TEST_sign = []
    File = open(file_name, "rb")
    poly = []
    for i in range(ntest):
            TEST_msg.append(list(File.read(n)))
    for i in range(ntest):
            TEST_sign.append(list(File.read(n)))
    File.close()
    return([TEST_msg,TEST_sign])


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