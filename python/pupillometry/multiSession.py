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

df = ['mZS066d20211011',
      'mZS068d20211009',
      'mZS068d20211010',
      'mZS068d20211011',
      'mZS069d20211009',
      'mZS069d20211010',
      'mZS069d20211011',
      'mZS070d20211009',
      'mZS070d20211010',
      'mZS070d20211011',
      'mZS071d20211009',
      'mZS071d20211010',
      'mZS071d20211011',
      'mZS072d20211009',
      'mZS072d20211010',
      'mZS072d20211011']
df = ['mZS066d20211013',
      'mZS068d20211013',
      'mZS069d20211013',
      'mZS070d20211013',
      'mZS071d20211013',
      'mZS072d20211013']

for i in range(0, len(df)):
    session = df[i]; 
    deeplabcut.extractPupil(session, cdnn, label = False);
    
#%%    