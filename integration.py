## Courtesy of Joshua Goings <https://joshuagoings.com/2017/04/28/integrals/>

import numpy as np
from scipy.special import hyp1f1

def calc_E(power0, power1, t, Q, a, b):

    p = a + b
    q = a * b / (a + b)

    if (t < 0) or (t > (power0 + power1)):
        return 0.0
        
    elif power0 == power1 == t == 0:
        return np.exp(-q * Q * Q)
        
    elif power1 == 0:
        return (1 / (2 * p)) * calc_E(power0 - 1, power1, t - 1, Q, a, b) - \
                (q * Q / a) * calc_E(power0 - 1, power1, t, Q, a, b)    + \
                (t + 1) * calc_E(power0 - 1, power1, t + 1, Q, a, b)
    else:
        return (1 / (2 * p)) * calc_E(power0, power1 - 1, t - 1, Q, a, b) + \
                (q * Q / b) * calc_E(power0, power1 - 1, t, Q, a, b)    + \
                (t + 1) * calc_E(power0, power1 - 1, t + 1, Q, a, b)
        

def overlap(a, b, N_coord0, N_coord1, powers0, powers1):

    E = 1.0
    for i in range(0, 3):

        Q = N_coord1[i] - N_coord0[i]

        E *= calc_E(powers0[i], powers1[i], 0, Q, a, b)

    result = E * np.power(np.pi / (a + b), 1.5)

    return result


def kinetic(a, b, N_coord0, N_coord1, powers0, powers1):

    l0, m0, n0 = powers0 
    l1, m1, n1 = powers1

    k0 = b * (2 * (l1 + m1 + n1) + 3) * \
         overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1, m1, n1))
    
    k1 = -2 * np.power(b,2) * \
         (overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1 + 2, m1, n1)) + \
          overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1, m1 + 2, n1)) + \
          overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1, m1, n1 + 2)))

    k2 = -0.5 * (l1 * (l1 - 1) * overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1 - 2, m1, n1)) + \
            m1 * (m1 - 1) * overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1, m1 - 2, n1)) + \
            n1 * (n1 - 1) * overlap(a, b, N_coord0, N_coord1, (l0, m0, n0), (l1, m1, n1 - 2)))

    return k0 + k1 + k2


def calc_R(t, u, v, n, p, PCx, PCy, PCz, RPC):

    T = p*RPC*RPC
    val = 0.0
    if t == u == v == 0:
        val += np.power(-2*p,n)*hyp1f1(n+0.5,n+1.5,-T)/(2.0*n+1.0)
    elif t == u == 0:
        if v > 1:
            val += (v-1)*calc_R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
        val += PCz*calc_R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
    elif t == 0:
        if u > 1:
            val += (u-1)*calc_R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
        val += PCy*calc_R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
    else:
        if t > 1:
            val += (t-1)*calc_R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
        val += PCx*calc_R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
    return val


def gaussian_product_center(a, A, b, B):

    P = []

    for i in range(0 , 3):
        P.append((a * A[i] + b * B[i]) / (a + b))

    return P


def nuclear(a, b, coord0, coord1, powers0, powers1, N_coord):

    l0, m0, n0 = powers0 
    l1, m1, n1 = powers1
    
    p = a + b

    P = []
    d = []
    
    for i in range(0 , 3):
        P.append((a * coord0[i] + b * coord1[i]) / (a + b))
        d.append(P[i] - N_coord[i])
    
    RPC = np.linalg.norm(d)

    val = 0.0
    for t in range(l0 + l1 + 1):
        for u in range(m0 + m1 + 1):
            for v in range(n0 + n1 + 1):
                
                val += calc_E(l0, l1, t, coord0[0] - coord1[0], a, b) * \
                       calc_E(m0, m1, u, coord0[1] - coord1[1], a, b) * \
                       calc_E(n0, n1, v, coord0[2] - coord1[2], a, b) * \
                       calc_R(t, u, v, 0, p, d[0], d[1], d[2], RPC)

    val *= 2 * np.pi / p 

    return val


def electron(a,lmn1,A,b,lmn2,B,c,lmn3,C,d,lmn4,D):

    l1,m1,n1 = lmn1
    l2,m2,n2 = lmn2
    l3,m3,n3 = lmn3
    l4,m4,n4 = lmn4
    p = a+b
    q = c+d
    alpha = p*q/(p+q)
    P = gaussian_product_center(a,A,b,B)
    Q = gaussian_product_center(c,C,d,D)
    diff = []
    for i in range(0, 3):
        diff.append(P[i] - Q[i])
    RPQ = np.linalg.norm(diff)

    val = 0.0
    for t in range(l1+l2+1):
        for u in range(m1+m2+1):
            for v in range(n1+n2+1):
                for tau in range(l3+l4+1):
                    for nu in range(m3+m4+1):
                        for phi in range(n3+n4+1):
                            val += calc_E(l1,l2,t,A[0]-B[0],a,b) * \
                                   calc_E(m1,m2,u,A[1]-B[1],a,b) * \
                                   calc_E(n1,n2,v,A[2]-B[2],a,b) * \
                                   calc_E(l3,l4,tau,C[0]-D[0],c,d) * \
                                   calc_E(m3,m4,nu ,C[1]-D[1],c,d) * \
                                   calc_E(n3,n4,phi,C[2]-D[2],c,d) * \
                                   np.power(-1,tau+nu+phi) * \
                                   calc_R(t+tau,u+nu,v+phi,0,\
                                       alpha,P[0]-Q[0],P[1]-Q[1],P[2]-Q[2],RPQ)

    val *= 2*np.power(np.pi,2.5)/(p*q*np.sqrt(p+q))
    return val
