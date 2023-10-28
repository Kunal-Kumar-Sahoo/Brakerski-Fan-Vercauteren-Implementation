import math
from helper import *

def NTT(A, W_table, q):
    n = len(A)
    B = [x for x in A]
    v = int(math.log(n, 2))

    for i in range(v):
        for j in range(2 ** i):
            for k in range(2 ** (v - i - 1)):
                s = j * (2 ** (v - i)) + k
                t = s + (2 ** (v - i - 1))

                w = W_table[((2 ** i) * k)]

                as_temp, at_temp = B[s], B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q
    
    B = index_reverse(B, v)
    return B

def INTT(A, W_table, q):
    n = len(A)
    B = [x for x in A]

    v = int(math.log(n, 2))

    for i in range(v):
        for j in range(2 ** i):
            for k in range(2 ** (v - i - 1)):
                s = j * (2 ** (v - i)) + k
                t = s + (2 ** (v - i - 1))

                w = W_table[((2 ** i) * k)]

                as_temp, at_temp = B[s], B[t]
                
                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q
    
    B = index_reverse(B, v)

    n_inv = modinv(n, q)
    for i in range(n):
        B[i] = (B[i] * n_inv) % q
    
    return B


