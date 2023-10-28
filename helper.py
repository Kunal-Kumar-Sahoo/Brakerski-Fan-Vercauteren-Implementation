from generate_prime import *
from random import randint

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)
    
def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('Modular inverse does not exist')
    else:
        return x % m
    
def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a

def int_reverse(a, n):
    b = ('{:0' + str(n) + 'b}').format(a)
    return int(b[::-1], 2)

def index_reverse(a, r):
    n = len(a)
    b = [0] * n
    for i in range(n):
        rev_idx = int_reverse(i, r)
        b[rev_idx] = a[i]
    return b

def ref_poly_mul(A, B, M):
    C = [0] * (2 * len(A))
    D = [0] * len(A)

    for i, a in enumerate(A):
        for j, b in enumerate(B):
            C[i + j] = (C[i + j] + a * b) % M

    for i in range(len(A)):
        D[i] = (C[i] - C[i + len(A)]) % M
    
    return D

def ref_poly_mul2(A, B):
    C = [0] * (2 * len(A))
    D = [0] * len(A)

    for i, a in enumerate(A):
        for j, b in enumerate(B):
            C[i + j] = (C[i + j] + a * b)
    
    for i in range(len(A)):
        D[i] = C[i] - C[i + len(A)]
    
    return D

def is_root_unity(w, m, q):
    if w == 0:
        return False
    elif pow(w, m // 2, q) == q - 1:
        return True
    return False

def get_proper_prime(n, logq):
    factor = 2 * n
    value = (1 << logq) - factor + 1
    lbound = 1 << (logq - 1)
    while value > lbound:
        if is_prime(value) == True:
            return value
        value -= factor
    raise Exception('Failed to find a proper prime.')

def find_primitive_root(m, q):
    g = (q - 1) // m
    
    if q - 1 != g * m:
        return False
    
    attempt_ctr, attempt_max = 0, 100

    while attempt_ctr < attempt_max:
        a = randint(2, q - 1)
        b = pow(a, g, q)
        if is_root_unity(b, m, q):
            return True
        attempt_ctr += 1
    
    return True, 0

def param_gen(n, logq):
    pfound = False
    while not pfound:
        q = get_proper_prime(n, logq)
        pfound, psi = find_primitive_root(2 * n, q)
    psiv = modinv(psi, q)
    w = pow(psi, 2, q)
    wv = modinv(w, q)

    return q, psi, psiv, w, wv