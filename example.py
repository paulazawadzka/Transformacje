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
xy2000 = np.zeros((rows, 2))
xy92 = np.zeros((rows, 2))

neu = np.zeros((rows, cols))