# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:27:20 2022

@author: pulka
"""

import numpy as np
from inf2_proj1 import *

el_grs80 = Transformacje(model = 'grs80')

plik = 'wsp_inp.txt'

tablica = np.genfromtxt(plik, delimiter =',', skip_header = 4)

rows, cols = np.shape(tablica)

blh = np.zeros((rows, cols))
xy00 = np.zeros((rows, 2))
xy92 = np.zeros((rows, 2))
neu = np.zeros((rows, cols))
az_elev_dis = np.zeros((rows, 4))

tablica_ze_wsp = np.zeros((rows, 14))

for i in range(rows):
    blh[i] = el_grs80.xyz2blh_hirvonen(tablica[i, 0], tablica[i, 1], tablica[i, 2])
    xy00[i] = el_grs80.uklad_2000(blh[i,0], blh[i, 1])
    xy92[i] = el_grs80.uklad_1992(blh[i,0], blh[i, 1])
    neu[i] = el_grs80.xyz2neu(tablica[i, 0], tablica[i, 1], tablica[i, 2], tablica[i, 0]+1, tablica[i, 1]+1, tablica[i, 2]+1)
    az_elev_dis[i] = el_grs80.az_elev_dis(tablica[i, 0], tablica[i, 1], tablica[i, 2], tablica[i, 0]+1, tablica[i, 1]+1, tablica[i, 2]+1)
    
    tablica_ze_wsp[i] = np.hstack([blh[i], xy00[i], xy92[i], neu[i], az_elev_dis[i]])


np.savetxt("wsp_out.txt", tablica_ze_wsp, delimiter=',', fmt = ['%.7f', '%11.7f', '%8.3f', '%12.3f', '%12.3f', '%11.3f', '%11.3f', '%10.7f', '%10.7f', '%11.7f', '%11.7f', '%12.7f', '%6.3f', '%6.3f'],encoding = 'UTF', header = 'Konwersja współrzednych geodezyjnych - Paula Zawadzka \n-- B -------- L ------- h ------ x2000 ------ y2000 ------ x1992 ------ y1992 ------ N -------- E -------- U -------- Az -------- elev -- 2D dist - 3D dist --')