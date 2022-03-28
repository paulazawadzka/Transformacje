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
        self.flattening = (self.a - self.b) / self.a
        self.e2 = 2 * self.flattening - self.flattening ** 2          
    
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
        
    
#TEST

X= 3664940.500
Y= 1409153.590
Z= 5009571.170
test = Transformacje("grs80")
fi, lam, hel = test.xyz2blh_hirvonen(X, Y, Z)
X1, Y1, Z1 = test.blh2xyz(fi, lam, hel)
print(X1, Y1, Z1)