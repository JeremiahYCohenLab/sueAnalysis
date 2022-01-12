# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 00:25:39 2020

@author: zhixi
"""
file = r'F:\inhibitionAll.xlsx'
sheet  = 'allInhibition'
col = 'inhibitionNrwd'
cdnn = 'edgeHiddenLEDinPupil-ZS-2021-05-19'

import pandas as pd
import deeplabcut
#%%
file = r'J:\PSpupil.xlsx'
sheet  = 'PS02pupil'
col = 'all'
cdnn = 'edgeHiddenLEDinPupil-ZS-2021-05-19'

df = pd.read_excel(file, sheet_name = sheet, usecols = [col])

for i in range(0, len(df)):
    session = df.loc[i].at[col];
    if session == session:
        deeplabcut.extractPupil(session, cdnn, label = False);
            
    
#%%
ani = ['ZS066', 'ZS068', 'ZS069', 'ZS070', 'ZS071'];
dates = ['20211027', '20211028', '20211029', '20211030', '20211031'];

for i in range(len(ani)):
    for j in range(len(dates)):
        session = f"m{ani[i]}d{dates[j]}";
        deeplabcut.extractPupil(session, cdnn, label = False)
    
    
    
#%%    