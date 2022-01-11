import math as m
#przeliczanie phi i lambda na x y G-K

a=6378137
e2=0.00669437999013

b2 = a**2*(1-e2)
e2b = (((a ** 2) - (b2)) / (b2))
A0= 1-(e2/4) - ((3*(e2)**2)/64)-((5*(e2)**3)/256)
A2=(3/8)*((e2+((e2**2)/4)+((15*(e2)**3))/128))
A4=(15/256)*(e2**2+((3*((e2)**3))/4))
A6=(35*((e2)**3))/3072
# A = lam1,fi2 B = lam2,fi2
#C = lam1,fi1 D = lam2,fi1
lam1 = 20.75
fi1 = 50.25
lam2 = 21.25
fi2 = 50


def U1992(fi,lam):
    sigma = a * (A0 * m.radians(fi) - A2 * m.sin(m.radians(2 * fi)) + A4 * m.sin(m.radians(4 * fi)) - A6 * m.sin(
        m.radians(6 * fi)))
    fi = m.radians(fi)
    t = m.tan(fi)
    ni2 = (e2b) * ((m.cos(fi)) ** 2)
    L0 = 19 * m.pi / 180  # w radianach
    lam = m.radians(lam)
    l = lam - L0
    N = a / (m.sqrt(1 - e2 * (m.sin(fi)) ** 2))

    xgk = sigma + ((l ** 2) / 2) * N * m.sin(fi) * m.cos(fi) * (
                1 + (l ** 2 / 12) * (m.cos(fi) ** 2) * (5 - t ** 2 + 9 * ni2 + 4 * (ni2 ** 2)) + ((l ** 4) / 360) * (
                    (m.cos(fi)) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * (ni2) - 330 * (ni2) * (t ** 2)))
    ygk = l * N * m.cos(fi) * (1 + ((l ** 2) / 6) * ((m.cos(fi)) ** 2) * (1 - (t ** 2) + (ni2)) + ((l ** 4) / 120) * (
                (m.cos(fi)) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * (ni2) - 58 * (ni2) * (t ** 2)))

    m01992 = 0.9993
    X1992 = m01992 * xgk - 5300000
    Y1992 = m01992 * ygk + 500000

    B = obliczB(xgk,m0)
    fi2 = B - t/2 * (((ygk/m0*N)**2) * (1 + ni2) - 1/12 * ((ygk/m0*N)**4) * (5 + 3*(t**2) + 6*ni2 -
        6*ni2*(t**2) - 3*(ni2**2) - 9*(t**2)*(ni2**2)) + 1/360 * ())












    #X1992 = round(X1992,3)
    #Y1992 = round(Y1992,3)
    return X1992,Y1992,xgk,ygk

def obliczB(x,m0):
    B = x / (a * A0 * m0)
    while True:
        B1 = x / (a * A0 * m0) + A2 / A0 * m.sin(2 * B) - A4 / A0 * m.sin(4 * B) + A6 / A0 * m.sin(6 * B)
        if abs(m.degrees(B1 - B)) * 3600 < 0.000001:
            return B1
        else:
            B = B1

def gk2fl (xgk,ygk,m0):
    B = obliczB(m0)
    fi = B -


# Ax1992,Ay1992 = U1992(fi2,lam1)
# Bx1992,By1992 = U1992(fi2,lam2)
# Cx1992,Cy1992 = U1992(fi1,lam1)
# Dx1992,Dy1992 = U1992(fi1,lam2)

# print("Wspołrzędne pkt A w 1992:",Ax1992,Ay1992)
# print("Wspołrzędne pkt B w 1992:",Bx1992,By1992)
# print("Wspołrzędne pkt C w 1992:",Cx1992,Cy1992)
# print("Wspołrzędne pkt D w 1992:",Dx1992,Dy1992,'\n')



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

    xgk = sigma + ((l ** 2) / 2) * N * m.sin(fi) * m.cos(fi) * (
                1 + (l ** 2 / 12) * (m.cos(fi) ** 2) * (5 - t ** 2 + 9 * ni2 + 4 * (ni2 ** 2)) + ((l ** 4) / 360) * (
                    (m.cos(fi)) ** 4) * (61 - 58 * (t ** 2) + (t ** 4) + 270 * (ni2) - 330 * (ni2) * (t ** 2)))
    ygk = l * N * m.cos(fi) * (1 + ((l ** 2) / 6) * ((m.cos(fi)) ** 2) * (1 - (t ** 2) + (ni2)) + ((l ** 4) / 120) * (
                (m.cos(fi)) ** 4) * (5 - 18 * (t ** 2) + (t ** 4) + 14 * (ni2) - 58 * (ni2) * (t ** 2)))

    m02000 = 0.999923
    X2000 = m02000*xgk
    Y2000 = m02000*ygk+1000000*nr+500000
    #X2000 = round(X2000,3)
    #Y2000 = round(Y2000,3)
    return X2000,Y2000,xgk,ygk

# Ax2000,Ay2000,nrA = U2000(fi2,lam1)
# Bx2000,By2000,nrB = U2000(fi2,lam2)
# Cx2000,Cy2000,nrC = U2000(fi1,lam1)
# Dx2000,Dy2000,nrD = U2000(fi2,lam2)
lolx1992,loly1992,xgk1992,ygk1992 = U1992(51.70102972777778,18.175462347222222)
lolx2000,loly2000,xgk2000,ygk2000 = U2000(51.70102972777778,18.175462347222222)


# print("Wspołrzędne pkt A w 2000:",Ax2000,Ay2000,nrA)
# print("Wspołrzędne pkt B w 2000:",Bx2000,By2000,nrB)
# print("Wspołrzędne pkt C w 2000:",Cx2000,Cy2000,nrC)
# print("Wspołrzędne pkt D w 2000:",Dx2000,Dy2000,nrD)

print("Wspołrzędne pkt lol w gk92:",xgk1992,ygk1992)
print("Wspołrzędne pkt lol w gk20:",xgk2000,ygk2000)

print("Wspołrzędne pkt lol w 2000:",lolx2000,loly2000)
print("Wspołrzędne pkt lol w 1992:",lolx1992,loly1992)
