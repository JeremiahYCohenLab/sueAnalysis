# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 00:25:39 2020

@author: zhixi
"""
file = r'D:\ZS056.xlsx'
sheet  = 'ZS056'
col = 'good'
cdnn = 'combine-withHiddenPoints-2020-03-27'


import pandas as pd
import deeplabcut


df = pd.read_excel(file, sheet_name = sheet, usecols = [col])


for i in range(0, len(df)):
    session = df[col][i]
    deeplabcut.extractPupil(session, cdnn)
    
#%%                