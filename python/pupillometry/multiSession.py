# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 00:25:39 2020

@author: zhixi
"""
file = r'F:\ZS062.xlsx'
sheet  = 'ZS062'
col = 'all'
cdnn = 'edgeHiddenLEDinPupil-ZS-2021-05-19'


import pandas as pd
import deeplabcut


df = pd.read_excel(file, sheet_name = sheet, usecols = [col])


for i in range(14, len(df)):
    session = df[col][i]; 
    deeplabcut.extractPupil(session, cdnn, label = False);
    
#%%    