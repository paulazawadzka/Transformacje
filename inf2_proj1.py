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
    
    def s_m_s(self, fi):
        fi = fi*180/math.pi
        d = np.floor(fi)
        m = np.floor((fi-d)*60)
        s = round((fi-d-m/60)*3600,5)
        print(d,'st',m,'min',s,'sek')
        
    def xyz2blh_hirvonen(self, X, Y, Z):
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
        return fi_n, lam, h
        
        
    def blh2xyz(self, fi, lam, h):
        N = self.a / math.sqrt(1 - self.e2 * (math.sin(fi))**2)
        X = (N + h) * math.cos(fi) * math.cos(lam)
        Y = (N + h) * math.cos(fi) * math.sin(lam)
        Z = (N * (1 - self.e2) + h) * math.sin(fi)
        return(X, Y, Z)
        
    def Vincent(self, fi, lam, fi_k, lam_k):
        b = self.a * np.sqrt(1 - self.e2)
        d_lam = lam_k -lam
        U_A = np.arctan((1 - self.f)*np.tan(fi))
        U_B = np.arctan((1 - self.f)*np.tan(fi_k))
        eps = 0.000001/3600 * math.pi/180 # radiany
        L=d_lam
     
        while True:
            
            sin_sgm = np.sqrt((np.cos(U_B)*np.sin(L))**2 + (np.cos(U_A)*np.sin(U_B) - np.sin(U_A)*np.cos(U_B)*np.cos(L))**2)
            cos_sgm = np.sin(U_A)*np.sin(U_B) + np.cos(U_A)*np.cos(U_B)*np.cos(L)
            sgm = np.arctan(sin_sgm/cos_sgm)
            
            sin_alfa = (np.cos(U_A)*np.cos(U_B)*np.sin(L))/sin_sgm
            cos2_alfa = 1 - sin_alfa**2
           
            cos_2sgm_m = cos_sgm - (2*np.sin(U_A)*np.sin(U_B))/cos2_alfa
          
            C = (self.f/16)*cos2_alfa*(4 + self.f*(4 - 3*cos2_alfa))
            
            L1 = L
            L= d_lam + (1 - C)*self.f*sin_alfa*(sgm + C*sin_sgm*(cos_2sgm_m + C*cos_sgm*((-1)+2*(cos_2sgm_m)**2)))
            
            if np.abs(L1 - L)<eps:
                break
        
        u2 = ((self.a**2 - b**2)/b**2)*cos2_alfa
        A = 1 + (u2/16384)*(4096 + u2*(-768+u2*(320-175*u2)))
        B = (u2/1024)*(256 + u2*(-128 + u2*(74-(47*u2))))
        
        d_sgm = B*sin_sgm*(cos_2sgm_m + (1/4)*B*(cos_sgm*(-1+2*(cos_2sgm_m**2)) - (1/6)*B*cos_2sgm_m*(-3+4*(sin_sgm**2))*(-3+4*(cos_2sgm_m**2))))
       
        s_AB = b*A*(sgm - d_sgm)
        s_AB = round(s_AB, 3)
        A_AB = np.arctan((np.cos(U_B)*np.sin(L))/(np.cos(U_A)*np.sin(U_B) - np.sin(U_A)*np.cos(U_B)*np.cos(L)))
        A_BA = np.arctan((np.cos(U_A)*np.sin(L))/(-np.sin(U_A)*np.cos(U_B) + np.cos(U_A)*np.sin(U_B)*np.cos(L)))+np.pi
        return s_AB, A_AB, A_BA
    
    def uklad_1992(self, fi, lam):
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
 
    
#TEST

X = 3664940.500
Y = 1409153.590
Z = 5009571.170
test = Transformacje("grs80")
fi, lam, hel = test.xyz2blh_hirvonen(X, Y, Z)
X1, Y1, Z1 = test.blh2xyz(fi, lam, hel)

XK = 3664000.500
YK = 1409000.590
ZK = 5009000.170
fi_k, lam_k, hel_k = test.xyz2blh_hirvonen(XK, YK, ZK)
s, A_AB, A_BA = test.Vincent(fi, lam, fi_k, lam_k)


XG = 3763917
YG = 1235727
ZG = 4982084+13*15
fiG, lamG, helG = test.xyz2blh_hirvonen(XG, YG, ZG)
X00, Y00 = test.uklad_2000(fiG, lamG)
X92, Y92 = test.uklad_1992(fiG, lamG)
print(X00, Y00)
print(X92, Y92)