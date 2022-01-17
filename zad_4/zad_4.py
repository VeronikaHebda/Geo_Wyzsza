import math as m
from shapely.geometry.polygon import Polygon
#przeliczanie phi i lambda na x y G-K

a=6378137
e2=0.00669437999013

b2 = a**2*(1-e2)
e2b = (((a ** 2) - (b2)) / (b2))
A0= 1-(e2/4) - ((3*(e2)**2)/64)-((5*(e2)**3)/256)
A2=(3/8)*((e2+((e2**2)/4)+((15*(e2)**3))/128))
A4=(15/256)*(e2**2+((3*((e2)**3))/4))
A6=(35*((e2)**3))/3072
#A = lam1,fi1 C = lam2,fi1
#B = lam1,fi2 D = lam2,fi2
lam1 = 20.75
fi1 = 50.25
lam2 = 21.25
fi2 = 50
#punkt sredniej szerokosci
fi_S = 50.125
lam_S = 21
#punkt srodkowy
fi_E = 50.12527044899897
lam_E = 21.00065108883653


def U1992(fi,lam):
    sigma = a * (A0 * m.radians(fi) - A2 * m.sin(m.radians(2 * fi)) + A4 * m.sin(m.radians(4 * fi)) - A6 * m.sin(
        m.radians(6 * fi)))
    fi = m.radians(fi)
    t = m.tan(fi)
    ni2 = e2b * ((m.cos(fi)) ** 2)
    L0 = 19 * m.pi / 180  # w radianach
    lam = m.radians(lam)
    l = lam - L0
    N = a / (m.sqrt(1 - e2 * (m.sin(fi)) ** 2))
    M = (a * (1 - e2)) / m.sqrt((1 - e2 * m.sin(m.radians(fi)) ** 2) ** 3)

    xgk = sigma + ((l ** 2) / 2) * N * m.sin(fi) * m.cos(fi) * (
                1 + (l ** 2 / 12) * (m.cos(fi) ** 2) * (5 - t ** 2 + 9 * ni2 + 4 * (ni2 ** 2)) + ((l ** 4) / 360) * (
                    (m.cos(fi)) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * (ni2) - 330 * (ni2) * (t ** 2)))
    ygk = l * N * m.cos(fi) * (1 + ((l ** 2) / 6) * ((m.cos(fi)) ** 2) * (1 - (t ** 2) + (ni2)) + ((l ** 4) / 120) * (
                (m.cos(fi)) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * (ni2) - 58 * (ni2) * (t ** 2)))

    m01992 = 0.9993
    X1992 = m01992 * xgk - 5300000
    Y1992 = m01992 * ygk + 500000

    X1992 = round(X1992,3)
    Y1992 = round(Y1992,3)

    R = m.sqrt(M*N)
    mgk = 1 + (ygk**2)/(2*(R**2)) + ygk**4/(24*(R**4))
    m1992 = mgk * m01992

    mgk_2 = mgk**2
    m92_2 = mgk_2 * 0.9993**2

    return X1992,Y1992,round(xgk,3),round(ygk,3),round(mgk,6),round(m1992,6),round(mgk_2,6),round(m92_2,6)


def gk2filam(x,y):
    B = x / (a * A0)
    while True:
        sigma = a * (A0 * B - A2 * m.sin(2 * B) + A4 * m.sin(4 * B) - A6 * m.sin(6 * B))
        B1 = B + (x - sigma) / a * A0
        if abs(m.degrees(B1 - B)) * 3600 < 0.000001:
            break
        else:
            B = B1

    t = m.tan(B1)
    ni2 = e2b * (m.cos(B1) ** 2)
    N = a / m.sqrt(1 - e2 * m.sin(B1) ** 2)
    M = (a * (1 - e2)) / m.sqrt((1 - e2 * m.sin(m.radians(B1)) ** 2) ** 3)
    L0 = 19 * m.pi / 180  # w radianach

    fi = B1 - (y ** 2 * t) / (2 * M * N) * (1 - (y ** 2) / (12 * N ** 2) *
        (5+3*t**2+ni2-9*ni2*t**2-4*ni2**2)+(y**4)/(360*N**4) * (61 + 90*t**2 + 45*t**4))
    lam = L0 + y/(N*m.cos(B1)) * (1 - y**2/(6*N**2) * (1+2*t**2+ni2) +
        y**4/(120*N**4) * (5+28*t**2+24*t**4+6*ni2+8*ni2*t**2))

    fi = m.degrees(fi)
    lam = m.degrees(lam)
    return fi, lam

def strefa(lam):
    if lam >= 22.5:
        return 8,24
    elif lam >= 19.5:
        return 7,21
    elif lam >= 16.5:
        return 6,18
    else:
        return 5,15

def U2000(fi,lam):
    sigma = a * (A0 * m.radians(fi) - A2 * m.sin(m.radians(2 * fi)) + A4 * m.sin(m.radians(4 * fi)) - A6 * m.sin(
        m.radians(6 * fi)))
    fi = m.radians(fi)
    t = m.tan(fi)
    ni2 = (e2b) * ((m.cos(fi)) ** 2)
    nr,x = strefa(lam)
    L0 = m.radians(x)   # w radianach
    lam = m.radians(lam)
    l = lam - L0
    N = a / (m.sqrt(1 - e2 * (m.sin(fi)) ** 2))
    M = (a * (1 - e2)) / m.sqrt((1 - e2 * m.sin(m.radians(fi)) ** 2) ** 3)

    xgk = sigma + ((l ** 2) / 2) * N * m.sin(fi) * m.cos(fi) * (
                1 + (l ** 2 / 12) * (m.cos(fi) ** 2) * (5 - t ** 2 + 9 * ni2 + 4 * (ni2 ** 2)) + ((l ** 4) / 360) * (
                    (m.cos(fi)) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * (ni2) - 330 * (ni2) * (t ** 2)))
    ygk = l * N * m.cos(fi) * (1 + ((l ** 2) / 6) * ((m.cos(fi)) ** 2) * (1 - (t ** 2) + (ni2)) + ((l ** 4) / 120) * (
                (m.cos(fi)) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * (ni2) - 58 * (ni2) * (t ** 2)))

    m02000 = 0.999923
    X2000 = m02000*xgk
    Y2000 = m02000*ygk+1000000*nr+500000
    X2000 = round(X2000,3)
    Y2000 = round(Y2000,3)

    R = m.sqrt(M * N)
    mgk = 1 + (ygk ** 2) / (2 * (R ** 2)) + ygk ** 4 / (24 * (R ** 4))
    m2000 = mgk * m02000

    mgk_2 = mgk**2
    m2000_2 = mgk_2 * 0.999923**2
    return X2000,Y2000,round(m2000,6),round(m2000_2,6)

#A = lam1,fi1 C = lam2,fi1
#B = lam1,fi2 D = lam2,fi2

Ax1992,Ay1992,Axgk,Aygk,Amgk,Am1992,Amgk_2,Am92_2 = U1992(fi1,lam1)
Bx1992,By1992,Bxgk,Bygk,Bmgk,Bm1992,Bmgk_2,Bm92_2 = U1992(fi2,lam1)
Cx1992,Cy1992,Cxgk,Cygk,Cmgk,Cm1992,Cmgk_2,Cm92_2 = U1992(fi1,lam2)
Dx1992,Dy1992,Dxgk,Dygk,Dmgk,Dm1992,Dmgk_2,Dm92_2 = U1992(fi2,lam2)
Ex1992,Ey1992,Exgk,Eygk,Emgk,Em1992,Emgk_2,Em92_2 = U1992(fi_E,lam_E)
Sx1992,Sy1992,Sxgk,Sygk,Smgk,Sm1992,Smgk_2,Sm92_2 = U1992(fi_S,lam_S)

print("Wspołrzędne pkt A w G-K:",Axgk,Aygk)
print("Wspołrzędne pkt B w G-K:",Bxgk,Bygk)
print("Wspołrzędne pkt C w G-K:",Cxgk,Cygk)
print("Wspołrzędne pkt D w G-K:",Dxgk,Dygk)
print("Wspołrzędne pkt E w G-K:",Exgk,Eygk)
print("Wspołrzędne pkt S w G-K:",Sxgk,Sygk,'\n')

print("Wspołrzędne pkt A w 1992:",Ax1992,Ay1992)
print("Wspołrzędne pkt B w 1992:",Bx1992,By1992)
print("Wspołrzędne pkt C w 1992:",Cx1992,Cy1992)
print("Wspołrzędne pkt D w 1992:",Dx1992,Dy1992)
print("Wspołrzędne pkt E w 1992:",Ex1992,Ey1992)
print("Wspołrzędne pkt S w 1992:",Sx1992,Sy1992,'\n')

#A = lam1,fi1 C = lam2,fi1
#B = lam1,fi2 D = lam2,fi2
Ax2000,Ay2000,Am2000,Am2000_2 = U2000(fi1,lam1)
Bx2000,By2000,Bm2000,Bm2000_2 = U2000(fi2,lam1)
Cx2000,Cy2000,Cm2000,Cm2000_2 = U2000(fi1,lam2)
Dx2000,Dy2000,Dm2000,Dm2000_2 = U2000(fi2,lam2)
Ex2000,Ey2000,Em2000,Em2000_2 = U2000(fi_E,lam_E)
Sx2000,Sy2000,Sm2000,Sm2000_2 = U2000(fi_S,lam_S)

print("Wspołrzędne pkt A w 2000:",Ax2000,Ay2000)
print("Wspołrzędne pkt B w 2000:",Bx2000,By2000)
print("Wspołrzędne pkt C w 2000:",Cx2000,Cy2000)
print("Wspołrzędne pkt D w 2000:",Dx2000,Dy2000)
print("Wspołrzędne pkt E w 2000:",Ex2000,Ey2000)
print("Wspołrzędne pkt S w 2000:",Sx2000,Sy2000,'\n')

#Axodwrocone,Ayodwrocone = gk2filam(Axgk,Aygk)
#Bxodwrocone,Byodwrocone = gk2filam(Bxgk,Bygk)
#lolfi92,lollam92 = gk2filam(xgk1992,ygk1992)

#obliczanie pól powierzchni
Pelipsoidalne = 994265196.074311
print("Pole elipsoidalne:", Pelipsoidalne)
PGK = Polygon([(Axgk,Aygk),(Bxgk,Bygk),(Dxgk,Dygk),(Cxgk,Cygk),(Axgk,Aygk)])
print("Pole G-K:",PGK.area)
P2000 = Polygon([(Ax2000,Ay2000),(Bx2000,By2000),(Dx2000,Dy2000),(Cx2000,Cy2000),(Ax2000,Ay2000)])
print("Pole 2000:",P2000.area)
P1992 = Polygon([(Ax1992,Ay1992),(Bx1992,By1992),(Dx1992,Dy1992),(Cx1992,Cy1992),(Ax1992,Ay1992)])
print("Pole 1992:",P1992.area,'\n')

#elementarne skali długości
mgk = [Amgk,Bmgk,Cmgk,Dmgk,Emgk,Smgk]
Kgk = []
m92 = [Am1992,Bm1992,Cm1992,Dm1992,Em1992,Sm1992]
K92 = []
m2000 = [Am2000,Bm2000,Cm2000,Dm2000,Em2000,Sm2000]
K2000 = []
for i in range(0,6):
    K_gk = (1 - mgk[i]) * 1000
    Kgk.append(round(K_gk,2))
    K_92 = (1 - m92[i]) * 1000
    K92.append(round(K_92,2))
    K_20 = (1 - m2000[i]) * 1000
    K2000.append(round(K_20,2))

print("Mgk:",mgk)
print("Kgk(1km):",Kgk)
print("M92:",m92)
print("K92(1km):",K92)
print("m2000:",m2000)
print("K2000(1km):",K2000,'\n')

mgk_2 = [Amgk_2,Bmgk_2,Cmgk_2,Dmgk_2,Emgk_2,Smgk_2]
Kgk_2 = []
m92_2 = [Am92_2,Am92_2,Am92_2,Am92_2,Am92_2,Am92_2]
K92_2 = []
m2000_2 = [Am2000_2,Bm2000_2,Cm2000_2,Dm2000_2,Em2000_2,Sm2000_2]
K2000_2 = []
for i in range(0,6):
    Kgk2 = (1 - mgk_2[i]) * 1000
    Kgk_2.append(round(Kgk2, 2))
    K922 = (1 - m92_2[i]) * 1000
    K92_2.append(round(K922, 2))
    K202 = (1 - m2000_2[i]) * 1000
    K2000_2.append(round(K202, 2))

print("Mgk:",mgk_2)
print("Kgk(1km):",Kgk_2)
print("M92:",m92_2)
print("K92(1km):",K92_2)
print("m2000:",m2000_2)
print("K2000(1km):",K2000_2,'\n')