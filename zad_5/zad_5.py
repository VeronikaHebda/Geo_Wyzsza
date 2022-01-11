import math as m
import numpy as np
x = 50000
y = 20000
z = 10000

#promien
#elipsoida krasowskiego
a = 6378245
e2 = 0.00669342

def Hirvonen (x,y,z,e2,a):
    r = m.sqrt(x**2+y**2)
    fi = m.atan((x/r) * (1-e2)**(-1))
    while True:
        N = a/(m.sqrt(1 - e2 * m.sin(fi)**2))
        h = r/m.cos(fi) - N
        fi_next = m.atan(z/r * (1 - e2 * N/(N+h))**(-1))
        if abs(fi_next - fi) < 0.000005 * m.pi /(180*3600):
            lam = m.atan(y/x)
            N_last = a/(m.sqrt(1 - e2 * m.sin(fi_next)**2))
            h_last = r/m.cos(fi_next) - N_last
            break
        else:
            fi = fi_next
    return lam,fi_next,h_last,N_last

def konwersja(dziesietne):
    d = int(dziesietne)
    m = int((dziesietne-d) * 60)
    s = (dziesietne - d - m/60)*3600
    s = round(s,5)
    s = str(s)
    d = str(d)
    m = str(m)
    if len(d) == 1:
        d = "0" + d
    if len(m) == 1:
        m = "0" + m
    if len(s) == 1:
        s = "0" + s
    word = d + '° ' + m + "' " + s + "''"
    return word

print("Współrzędne xyz początkowe:", x,y,z)


lam,fi,h,N = Hirvonen(x,y,z,e2,a)
lam2 = konwersja(np.rad2deg(lam))
fi2 = konwersja(np.rad2deg(fi))

print(lam2,fi2,h)

def geo2xyz(fi,lam,h,e2,a):
    N = a/(m.sqrt(1 - e2 * m.sin(fi)**2))
    x = (N + h) * m.cos(fi) * m.cos(lam)
    y = (N + h) * m.cos(fi) * m.sin(lam)
    z = (N * (1 - e2) + h) * m.sin(fi)
    return round(x,3), round(y,3), round(z,3)

x2,y2,z2 = geo2xyz(fi,lam,h,e2,a)
print(x2,y2,z2)