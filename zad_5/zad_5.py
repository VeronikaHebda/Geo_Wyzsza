import math as m
import numpy as np

#elipsoida krasowskiego
a = 6378245
e2 = 0.00669342

#elipsoida grs80
a80=6378137
e280=0.00669437999013


#fi lambda h grs80 z cwiczenia 4
#A = lam1,fi1 B = lam2,fi1
#C = lam1,fi2 D = lam2,fi2
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
h = 0

def geo2xyz(fi, lam, h, a, e2):
    fi = m.radians(fi)
    lam = m.radians(lam)
    N = a/m.sqrt(1-e2*m.sin(fi)**2)
    x = (N+h)* m.cos(fi) * m.cos(lam)
    y = (N + h) * m.cos(fi) * m.sin(lam)
    z = (N*(1-e2)+h)*m.sin(fi)
    return x,y,z
def transformacja(x, y, z):
    x0 = -33.4297
    y0 = 146.5746
    z0 = 76.2865
    kappa = 1 + 0.8407728 * 10 ** -6
    Ex = np.deg2rad(-0.35867 / 3600)
    Ey = np.deg2rad(-0.05283 / 3600)
    Ez = np.deg2rad(0.84354 / 3600)
    macierz = np.array([[kappa, Ez, -Ey],
              [-Ez, kappa, Ex],
              [Ey, -Ex, kappa]])
    macierz_p = np.array([x,
                y,
                z])
    macierz_0 = np.array([x0,
                y0,
                z0])
    macierz_nowychwspol = macierz_p + macierz @ macierz_p + macierz_0

    return macierz_nowychwspol

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
    return lam,fi_next,h_last

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
#A = lam1,fi1 B = lam2,fi1
#C = lam1,fi2 D = lam2,fi2

print("Współrzędne punktu A (phi, lambda, h) na elipsoidzie GRS80:",fi1,lam1,h)
print("Współrzędne punktu B (phi, lambda, h) na elipsoidzie GRS80:",fi1,lam2,h)
print("Współrzędne punktu C (phi, lambda, h) na elipsoidzie GRS80:",fi2,lam1,h)
print("Współrzędne punktu D (phi, lambda, h) na elipsoidzie GRS80:",fi2,lam2,h)
print("Współrzędne punktu E (phi, lambda, h) na elipsoidzie GRS80:",fi_E,lam_E,h)
print("Współrzędne punktu S (phi, lambda, h) na elipsoidzie GRS80:",fi_S,lam_S,h)
#filamh grs80 na xyz grs80
x_A,y_A,z_A = geo2xyz(fi1,lam1,h,a80,e280)
x_B,y_B,z_B = geo2xyz(fi1,lam2,h,a80,e280)
x_C,y_C,z_C = geo2xyz(fi2,lam1,h,a80,e280)
x_D,y_D,z_D = geo2xyz(fi2,lam2,h,a80,e280)
x_E,y_E,z_E = geo2xyz(fi_E,lam_E,h,a80,e280)
x_S,y_S,z_S = geo2xyz(fi_S,lam_S,h,a80,e280)
print("--------------------------------------------------------------------")
print("Współrzędne punktu A (x,y,z) na elipsoidzie GRS80:", x_A,",",y_A,",",z_A)
print("Współrzędne punktu B (x,y,z) na elipsoidzie GRS80:", x_B,",",y_B,",",z_B)
print("Współrzędne punktu C (x,y,z) na elipsoidzie GRS80:", x_C,",",y_C,",",z_C)
print("Współrzędne punktu D (x,y,z) na elipsoidzie GRS80:", x_D,",",y_D,",",z_D)
print("Współrzędne punktu E (x,y,z) na elipsoidzie GRS80:", x_E,",",y_E,",",z_E)
print("Współrzędne punktu S (x,y,z) na elipsoidzie GRS80:", x_S,",",y_S,",",z_S)
print("--------------------------------------------------------------------")

#transformacja
[xk_A,yk_A,zk_A] = transformacja(x_A,y_A,z_A)
[xk_B,yk_B,zk_B] = transformacja(x_B,y_B,z_B)
[xk_C,yk_C,zk_C] = transformacja(x_C,y_C,z_C)
[xk_D,yk_D,zk_D] = transformacja(x_D,y_D,z_D)
[xk_E,yk_E,zk_E] = transformacja(x_E,y_E,z_E)
[xk_S,yk_S,zk_S] = transformacja(x_S,y_S,z_S)
print("Współrzędne punktu A (x,y,z) na elipsoidzie Krakowskiego::", xk_A,",",yk_A,",",zk_A)
print("Współrzędne punktu B (x,y,z) na elipsoidzie Krakowskiego::", xk_B,",",yk_B,",",zk_B)
print("Współrzędne punktu C (x,y,z) na elipsoidzie Krakowskiego::", xk_C,",",yk_C,",",zk_C)
print("Współrzędne punktu D (x,y,z) na elipsoidzie Krakowskiego::", xk_D,",",yk_D,",",zk_D)
print("Współrzędne punktu E (x,y,z) na elipsoidzie Krakowskiego::", xk_E,",",yk_E,",",zk_E)
print("Współrzędne punktu S (x,y,z) na elipsoidzie Krakowskiego::", xk_S,",",yk_S,",",zk_S)
print("--------------------------------------------------------------------")
#xyz krasowkiego na flh krasowskiego
lamk_A,fik_A,hk_A = Hirvonen(xk_A,yk_A,zk_A,e2,a)
#lamk = np.rad2deg(lamk)
#fik = np.rad2deg(fik)

#print("Współrzędne fi lam h krasowskiego:",lamk,fik,hk)

