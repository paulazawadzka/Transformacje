# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:30:11 2022

@author: pulka
"""

#INFORMATYKA 2 PROJEKT 1

import math
import numpy as np

class Transformacje:
    def __init__(self, model: str = "grs80"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            f - spłaszczenie
            e2 - mimośród^2
        """
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424528
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.f = (self.a - self.b) / self.a
        self.e2 = 2 * self.f - self.f ** 2          
    
    def d_m_s(self, fi):
        """
        Zamiana radianów na stopnie, minuty, sekundy.

        Parameters
        ----------
        fi : FLOAT
            wartość kąta [rad]

        """
        fi = fi*180/math.pi
        d = np.floor(fi)
        m = np.floor((fi-d)*60)
        s = round((fi-d-m/60)*3600,5)
        print(d,'st',m,'min',s,'sek')
        
    def xyz2blh_hirvonen(self, X, Y, Z):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokość elipsoidalna (fi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotnej iteracji wyznaczenia wsp. fi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie ortokartezjańskim 

        Returns
        -------
        fi_n : FLOAT
            szerokość geodezyjna [stopnie dziesiętne]
        lam : FLOAT
            długość geodezyjna [stopnie dziesiętne]
        h : FLOAT
            wysokość elipsoidalna [m]
        """
        r = math.sqrt(X**2 + Y**2)
        fi_n = math.atan(Z/(r * (1 - self.e2)))
        eps = 0.000001/3600 * math.pi/180
        fi = fi_n * 2
        while np.abs(fi_n-fi) > eps:
            fi = fi_n
            N = self.a/np.sqrt(1 - self.e2 * np.sin(fi_n)**2)
            h = r/np.cos(fi_n) - N
            fi_n = math.atan(Z/(r * (1 - self.e2 * (N/(N+h)))))
        lam = math.atan(Y/X)
        N = self.a/np.sqrt(1 - self.e2 * np.sin(fi_n)**2)
        h = r/np.cos(fi_n) - N
        fi_n = fi_n*180/math.pi
        lam = lam*180/math.pi 
        return fi_n, lam, h
        
        
    def blh2xyz(self, fi, lam, h):
        """
        Zamiana współrzędnych geodezyjnych długości szerokości i wysokości elipsoidalnej (phi, lam, h)
        na współrzędnych ortokartezjańskich (x, y, z). Działanie odwrotne do algorytmu Hirvonena.

        Parameters
        ----------
        fi : FLOAT
            szerokość geodezyjna [stopnie dziesiętne]
        lam : FLOAT
            długość geodezyjna [stopnie dziesiętne]
        h : FLOAT
            wysokość elipsoidalna [m]

        Returns
        -------
        X, Y, Z : FLOAT
             współrzędne w układzie ortokartezjańskim

        """
        fi = fi*math.pi/180
        lam = lam*math.pi/180
        N = self.a / math.sqrt(1 - self.e2 * (math.sin(fi))**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N * (1 - self.e2) + h) * math.sin(fi)
        return X, Y, Z
        
    
    def uklad_1992(self, fi, lam):
        """
        Zamiana współrzędnych geodezyjnych długości szerokości i wysokości elipsoidalnej (phi, lam, h)
        na współrzędnych płaskich układu polskiego 1992 (x, y).

        Parameters
        ----------
        fi : FLOAT
            szerokość geodezyjna [stopnie dziesiętne]
        lam : FLOAT
            długość geodezyjna [stopnie dziesiętne]

        Returns
        -------
        x92, y92 : FLOAT
            współrzędne w układzie 1992

        """
        fi = fi*math.pi/180
        lam = lam*math.pi/180
        m_0 = 0.9993
        N = self.a/(math.sqrt(1 - self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        e_2 = self.e2/(1 - self.e2)
        n2 = e_2 * np.cos(lam)**2
        lam_0 = math.radians(19)
        l = lam - lam_0
        
        A_0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A_2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A_4 = 15/256 * (self.e2**2 + (3*(self.e2**3))/4)
        A_6 = (35*(self.e2**3))/3072
        
        sigma = self.a * ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))

        x92 = round(x * m_0 - 5300000, 3)
        y92 = round(y * m_0 + 500000, 3)   
        
        return x92, y92 
 
    def uklad_2000(self, fi, lam):
        """
        Zamiana współrzędnych geodezyjnych długości szerokości i wysokości elipsoidalnej (phi, lam, h)
        na współrzędnych płaskich układu polskiego 2000 (x, y).

        Parameters
        ----------
        fi : FLOAT
            szerokość geodezyjna [stopnie dziesiętne]
        lam : FLOAT
            długość geodezyjna [stopnie dziesiętne]

        Returns
        -------
        x00, y0 : FLOAT
            współrzędne w układzie 2000

        """
        fi = fi*math.pi/180
        lam = lam*math.pi/180
        m_0 = 0.999923
        N = self.a/(math.sqrt(1 - self.e2 * np.sin(fi)**2))
        t = np.tan(fi)
        e_2 = self.e2/(1 - self.e2)
        n2 = e_2 * np.cos(lam)**2
        lam = math.degrees(lam)
        
        if lam > 13.5 and lam < 16.5:
            s = 5
            lam_0 = 15
        elif lam > 16.5 and lam < 19.5:
            s = 6
            lam_0 = 18
        elif lam > 19.5 and lam < 22.5:
            s = 7
            lam_0 = 21
        elif lam > 22.5 and lam < 25.5:
            s = 8
            lam_0 = 24
            
        lam = math.radians(lam)
        lam_0 = math.radians(lam_0)
        l = lam - lam_0
        
        A_0 = 1 - (self.e2/4) - (3*(self.e2**2))/64 - (5*(self.e2**3))/256
        A_2 = 3/8 * (self.e2 + ((self.e2**2)/4) + ((15*self.e2**3)/128))
        A_4 = 15/256 * (self.e2**2 + (3*(self.e2**3))/4)
        A_6 = (35*(self.e2**3))/3072
        
        
        sigma = self.a* ((A_0*fi) - (A_2*np.sin(2*fi)) + (A_4*np.sin(4*fi)) - (A_6*np.sin(6*fi)))
        
        x = sigma + ((l**2)/2) * (N*np.sin(fi)*np.cos(fi)) * (1 + ((l**2)/12) * ((np.cos(fi))**2) * (5 - t**2 + 9*n2 + (4*n2**2)) + ((l**4)/360) * ((np.cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*np.cos(fi)) * (1 + ((((l**2)/6) * (np.cos(fi))**2) * (1-(t**2) + n2)) +  (((l**4)/(120)) * (np.cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))

        x00 = round(x * m_0, 3)
        y00 = round(y * m_0 + (s*1000000) + 500000, 3)   
        
        return x00, y00 
 
    def xyz2neu(self, X, Y, Z, X0, Y0, Z0):
        """
        Zamiana współrzędnych ortokartezjańskich (x, y, z) na współrzędne topocentryczne (N, E, U).

        Parameters
        ----------
        X, Y, Z : FLOAT
            współrzędne ortokartezjańskie punktu
        X0, Y0, Z0 : FLOAT
            współrzędne ortokartezjańskie drugiego punktu

        Returns
        -------
        N, E, U : FLOAT
            współrzędne topocentryczne

        """
        delta_X = X - X0
        delta_Y = Y - Y0
        delta_Z = Z - Z0
            
        fi, lam, h = self.xyz2blh_hirvonen(X, Y, Z)
        R = np.array([[-np.sin(fi) * np.cos(lam), -np.sin(lam), np.cos(fi) * np.cos(lam)],
                         [-np.sin(fi) * np.sin(lam), np.cos(lam), np.cos(fi) * np.sin(lam)],
                         [np.cos(fi), 0, np.sin(fi)]])     
        d = np.matrix([delta_X, delta_Y, delta_Z])
        d = d.T
        NEU = R * d
            
        N = float(NEU[0])
        E = float(NEU[1])
        U = float(NEU[2])
        return N, E, U
 
    def az_elev_dis(self, X, Y, Z, X0, Y0, Z0):
        """
        Oblicznie azymutu, kąta elewacji, odległości 2D oraz 3D w układzie topocentrycznym (N, E, U)

        Parameters
        ----------
        X, Y, Z : FLOAT
            współrzędne ortokartezjańskie punktu
        X0, Y0, Z0 : FLOAT
            współrzędne ortokartezjańskie drugiego punktu

        Returns
        -------
        Az : FLOAT
            azymut [stopnie dziesiętne]
        elev : FLOAT
            kąt elewacji [stopnie dziesiętne]
        HZdist : FLOAT
            odległość 2D
        ELdist : FLOAT
            odległość 3D

        """
        
        N, E, U = self.xyz2neu(X, Y, Z, X0, Y0, Z0)
        
        HZdist = math.sqrt(E**2 + N**2)
        ELdist = math.sqrt(N**2 + E**2 + U**2)
        Az = math.atan2(E, N)
        elev = math.atan2(U, HZdist)
        if Az < 0:
            Az = Az + 2*math.pi
            
        Az = Az*180/math.pi
        elev = elev*180/math.pi
        return Az, elev, HZdist, ELdist
 
if __name__ == "__main__":
    # utworzenie obiektu
    test = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # transformacje
    fi, lam, h = test.xyz2blh_hirvonen(X, Y, Z)
    print(f"blh: {fi:.7f}, {lam:.7f}, {h:.3f}")
    x, y, z = test.blh2xyz(fi, lam, h)
    print(f"xyz: {x:.3f}, {y:.3f}, {z:.3f}")
    x00, y00 = test.uklad_2000(fi, lam)
    print(f"xy2000: {x00:.3f}, {y00:.3f}")
    x92, y92 = test.uklad_1992(fi, lam)
    print(f"xy1992: {x92:.3f}, {y92:.3f}")   
    N, E, U = test.xyz2neu(X, Y, Z, X+1, Y+1, Z+1)
    print(f"NEU: {N:.6f}, {E:.6f}, {U:.6f}")
    az, el, hor_d, v_d = test.az_elev_dis(X, Y, Z, X+1, Y+1, Z+1)
    print(f"az = {az:.6f}, elev = {el:.6f}, 2D distance = {hor_d:.3f}, 3D distance = {v_d:.3f}")
    