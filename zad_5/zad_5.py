import math as m
import numpy as np

#elipsoida krasowskiego
a = 6378245
e2 = 0.0066934215520

#elipsoida grs80
a80=6378137
e280=0.00669437999013


#fi lambda h grs80 z cwiczenia 4
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
    kappa = 0.8407728 * 10 ** -6
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
    r = (x ** 2 + y ** 2) ** 0.5
    fi = np.arctan((z / r) * (1 - e2) ** -1)

    n = a / m.sqrt(1 - e2 * np.sin(fi) ** 2)
    h = r / np.cos(fi) - n
    fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    sek = np.deg2rad(0.00005 / 3600)
    while abs(fi2 - fi) >= sek:
        fi = fi2
        n = a / m.sqrt(1 - e2 * np.sin(fi) ** 2)
        h = r / np.cos(fi) - n
        fi2 = np.arctan((z / r) * (1 - e2 * (n / (n + h))) ** -1)
    n = a / m.sqrt(1 - e2 * np.sin(fi2) ** 2)
    h = r / np.cos(fi2) - n
    lam = np.arctan(y / x)
    return fi2,lam,h

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

print("Współrzędne punktu A (phi, lambda, h) na elipsoidzie GRS80:",konwersja(fi1),konwersja(lam1),h)
print("Współrzędne punktu B (phi, lambda, h) na elipsoidzie GRS80:",konwersja(fi2),konwersja(lam1),h)
print("Współrzędne punktu C (phi, lambda, h) na elipsoidzie GRS80:",konwersja(fi1),konwersja(lam2),h)
print("Współrzędne punktu D (phi, lambda, h) na elipsoidzie GRS80:",konwersja(fi2),konwersja(lam2),h)
print("Współrzędne punktu E (phi, lambda, h) na elipsoidzie GRS80:",konwersja(fi_E),konwersja(lam_E),h)
print("Współrzędne punktu S (phi, lambda, h) na elipsoidzie GRS80:",konwersja(fi_S),konwersja(lam_S),h)
#filamh grs80 na xyz grs80
x_A,y_A,z_A = geo2xyz(fi1,lam1,h,a80,e280)
x_B,y_B,z_B = geo2xyz(fi2,lam1,h,a80,e280)
x_C,y_C,z_C = geo2xyz(fi1,lam2,h,a80,e280)
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
print("Współrzędne punktu A (x,y,z) na elipsoidzie Krasowskiego::", xk_A,",",yk_A,",",zk_A)
print("Współrzędne punktu B (x,y,z) na elipsoidzie Krasowskiego::", xk_B,",",yk_B,",",zk_B)
print("Współrzędne punktu C (x,y,z) na elipsoidzie Krasowskiego::", xk_C,",",yk_C,",",zk_C)
print("Współrzędne punktu D (x,y,z) na elipsoidzie Krasowskiego::", xk_D,",",yk_D,",",zk_D)
print("Współrzędne punktu E (x,y,z) na elipsoidzie Krasowskiego::", xk_E,",",yk_E,",",zk_E)
print("Współrzędne punktu S (x,y,z) na elipsoidzie Krasowskiego::", xk_S,",",yk_S,",",zk_S)
print("--------------------------------------------------------------------")
#xyz krasowkiego na flh krasowskiego
fik_A,lamk_A,hk_A = Hirvonen(xk_A,yk_A,zk_A,e2,a)
fik_B,lamk_B,hk_B = Hirvonen(xk_B,yk_B,zk_B,e2,a)
fik_C,lamk_C,hk_C = Hirvonen(xk_C,yk_C,zk_C,e2,a)
fik_D,lamk_D,hk_D = Hirvonen(xk_D,yk_D,zk_D,e2,a)
fik_E,lamk_E,hk_E = Hirvonen(xk_E,yk_E,zk_E,e2,a)
fik_S,lamk_S,hk_S = Hirvonen(xk_S,yk_S,zk_S,e2,a)

lamk_A = np.rad2deg(lamk_A)
lamk_B = np.rad2deg(lamk_B)
lamk_C = np.rad2deg(lamk_C)
lamk_D = np.rad2deg(lamk_D)
lamk_E = np.rad2deg(lamk_E)
lamk_S = np.rad2deg(lamk_S)

fik_A = np.rad2deg(fik_A)
fik_B = np.rad2deg(fik_B)
fik_C = np.rad2deg(fik_C)
fik_D = np.rad2deg(fik_D)
fik_E = np.rad2deg(fik_E)
fik_S = np.rad2deg(fik_S)

print("Współrzędne punktu A (phi, lambda, h) na elipsoidzie Krasowskiego:",konwersja(fik_A),konwersja(lamk_A),hk_A)
print("Współrzędne punktu B (phi, lambda, h) na elipsoidzie Krasowskiego:",konwersja(fik_B),konwersja(lamk_B),hk_B)
print("Współrzędne punktu C (phi, lambda, h) na elipsoidzie Krasowskiego:",konwersja(fik_C),konwersja(lamk_C),hk_C)
print("Współrzędne punktu D (phi, lambda, h) na elipsoidzie Krasowskiego:",konwersja(fik_D),konwersja(lamk_D),hk_D)
print("Współrzędne punktu E (phi, lambda, h) na elipsoidzie Krasowskiego:",konwersja(fik_E),konwersja(lamk_E),hk_E)
print("Współrzędne punktu S (phi, lambda, h) na elipsoidzie Krasowskiego:",konwersja(fik_S),konwersja(lamk_S),hk_S)

