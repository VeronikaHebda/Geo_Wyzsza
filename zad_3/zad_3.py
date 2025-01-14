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

#punkt sredniej szerokosci
fi_aryt = (fi1 + fi2)/2
lam_aryt = (lam1 + lam2)/2


#algorytm iteracyjny vincent
def vincent(lam1,fi1,lam2,fi2):
    Ua = m.atan((1 - f) * m.tan(m.radians(fi1)))
    Ub = m.atan((1 - f) * m.tan(m.radians(fi2)))
    delta_lam = lam2 - lam1
    L = m.radians(delta_lam)
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
            u2 = ((a ** 2 - b ** 2) / b ** 2) * cosa2
            A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
            B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
            delta_sigma = B * sino * (cos2om + 0.25 * B * (coso * (-1 + 2 * (cos2om ** 2)) - (1 / 6) * B * cos2om * (-3 + 4 * (sina ** 2)) * (-3 + 4 * (cos2om ** 2))))
            Scb = b * A * (o - delta_sigma)
            x = m.cos(Ub) * m.sin(L_last)
            y = m.cos(Ua) * m.sin(Ub) - m.sin(Ua) * m.cos(Ub) * m.cos(L_last)
            if (x > 0 and y > 0):
                Aab = m.atan(x / y)
            elif (x > 0 and y < 0):
                Aab = m.atan(x / y) + m.pi
            elif (x < 0 and y < 0):
                Aab = m.atan(x / y) + m.pi
            elif (x < 0 and y > 0):
                Aab = m.atan(x / y) + 2* m.pi

            c = m.cos(Ua) * m.sin(L_last)
            d = -m.sin(Ua) * m.cos(Ub) + m.sin(Ub) * m.cos(Ua) * m.cos(L_last)
            if (c > 0 and d > 0):
                Aba = m.atan(c / d) + m.pi
            elif (c > 0 and d < 0):
                Aba = m.atan(c / d) + 2* m.pi
            elif (c < 0 and d < 0):
                Aba = m.atan(c / d) + 2* m.pi
            elif (c < 0 and d > 0):
                Aba = m.atan(c / d) + 3* m.pi
            break
        else:
            L = L1
    return Scb, Aab, Aba


#algorytm Kivioja
def M(phi):
    return (a*(1 - e2)) / m.sqrt((1 - e2*(m.sin(m.radians(phi))**2))**3)

def N(phi):
    return a / m.sqrt(1 - e2*(m.sin(m.radians(phi))**2))

def Kivioji (S,fi,lam,Az):
    S = S/2
    ds = 1000
    n = int(S/1000)
    ostatni = S%1000

    F = []
    F.append(fi)
    A = []
    A.append(Az)
    L = []
    L.append(lam)

    #ds to jeden fragment
    for i in range(0,n+1):
        dfi = ds*m.cos(A[i])/M(F[i])
        dA = m.sin(A[i])*m.tan(m.radians(F[i]))*ds/N(F[i])

        Bm = F[i] + 0.5 * m.degrees(dfi)
        Aabm = A[i] + 0.5 * dA

        dBm = (m.cos(Aabm) * ds) / M(Bm)
        dlamb=(m.sin(Aabm)*ds)/(N(Bm)*m.cos(m.radians(Bm)))

        dAm=(m.sin(Aabm)*m.tan(m.radians(Bm))*ds)/N(Bm)

        F.append(F[i]+m.degrees(dBm))
        L.append(L[i]+m.degrees(dlamb))
        A.append(A[i]+dAm)
        if i+1 == n:
            ds = ostatni
    fi_koncowe = F[-1]
    lam_koncowe = L[-1]
    return fi_koncowe,lam_koncowe

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

#uzycie fukncji

print("Punkt średniej szerokości:",konwersja(fi_aryt),",",konwersja(lam_aryt))

#środek AD
Sad,Aab,Aba = vincent(lam1,fi1,lam2,fi2)
fi_ad,lam_ad = Kivioji(Sad,fi1,lam1,Aab)
print("Wspolrzedne punktu srodkowego AD - E:",konwersja(fi_ad),",",konwersja(lam_ad))

#srodek CB
Scb,Acb,Abc = vincent(lam1,fi2,lam2,fi1)
fi_cb,lam_cb = Kivioji(Scb,fi2,lam1,Acb)

#azymuty AD
Aab = np.rad2deg(Aab)
Aba = np.rad2deg(Aba)
print("Azymuty AD:", konwersja(Aab),",",konwersja(Aba))


#odleglosc i azymuty miedzy dwoma wyznacoznymi punktami

Se,Ax,Axy = vincent(lam_aryt,fi_aryt,lam_ad,fi_ad)
Se = round(Se,3)

Sed,Aed,Ade = vincent(lam_ad,fi_ad,lam2,fi2)
Aed = np.rad2deg(Aed)
Aed = konwersja(Aed)

Sxd,Axd,Adx = vincent(lam_aryt,fi_aryt,lam2,fi2)
Axd = np.rad2deg(Axd)
Axd = konwersja(Axd)

print("Odleglosc między punktem środkowym AD i punktem średniej szerokości:", Se,"m")

print("Azymut między punktem E i D:", Aed)
print("Azymut między punktem średniej szerokości a D:", Axd)
#obliczanie Pola czworokąta
e = m.sqrt(e2)
fi1 = np.deg2rad(fi1)
fi2 = np.deg2rad(fi2)
lam1 = np.deg2rad(lam1)
lam2 = np.deg2rad(lam2)
P = (b**2*(lam2-lam1)/2)*(((m.sin(fi2)/(1-e2*(m.sin(fi2)**2)))+(1/(2*e))*m.log((1+e*m.sin(fi2))/(1-e*m.sin(fi2))))
- (((m.sin(fi1)/(1-e2*(m.sin(fi1)**2)))+(1/(2*e))*m.log((1+e*m.sin(fi1))/(1-e*m.sin(fi1))))))
P = abs(P)
P = round(P,6)

print("Pole prostokąta:",P,"m²")
print("*Wspolrzedne punktu srodkowego CB:",konwersja(fi_cb),",",konwersja(lam_cb))
