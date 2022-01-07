#Weronika Hebda 311532

import math as m
import numpy as np

#wspolrzedne punktow
# A = lam1,fi1 B = lam2,fi1
#C = lam1,fi2 D = lam2,fi2
lam1 = 20.75
fi1 = 50.25
lam2 = 21.25
fi2 = 50

#bedzie liczona odleglosc punkty C i B


#dane
a=6378137
e2=0.00669437999013
b = a * m.sqrt(1 - e2)
f = 1 - b/a


#algorytm iteracyjny vincent
def vincent(lam1,fi1,lam2,fi2):
    Ua = m.atan((1 - f) * m.tan(m.radians(fi1)))
    Ub = m.atan((1 - f) * m.tan(m.radians(fi2)))
    delta_lam = lam2 - lam1
    L = delta_lam
    while True:
        sino = m.sqrt((m.cos(Ub)*m.sin(L))**2 + (m.cos(Ua)*m.sin(Ub) - m.sin(Ua)*m.cos(Ub)*m.cos(L))**2)
        coso = m.sin(Ua)*m.sin(Ub) + m.cos(Ua)*m.cos(Ub)*m.cos(L)
        o = m.atan(sino/coso)
        sina = (m.cos(Ua)*m.cos(Ub)*m.sin(L)/sino)
        cosa2 = 1 - (sina)**2
        cos2om = coso - (2*m.sin(Ua)*m.sin(Ub))/cosa2
        C = (f/16)*cosa2*(4+f*(4-3*cosa2))
        L1 = m.radians(delta_lam) + (1 - C)*f*sina*(o+C*sino*(cos2om+C*coso*(-1+2*((cos2om)**2))))
        if abs(m.degrees(L1-L)) * 3600 < 0.000001:
            L_last = L1
            break
        else:
            L = L1
    return L_last

L = vincent(lam1,fi1,lam2,fi2)

print(L)