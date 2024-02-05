# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 00:25:39 2020

@author: zhixi
"""
file = r'F:\activationAll.xlsx'
sheet  = 'all'
col = 'activation'
cdnn = 'edgeHiddenLEDinPupil-ZS-2021-05-19'
cdnn = 'pupilEdge-Katie-2021-12-08'
cdnn  = 'acutePupil-ZS-2022-06-21'

import pandas as pd
import numpy as np
import deeplabcut
#%%
file = r'F:\activationAll.xlsx'
sheet  = '689515'
col = 'all'
cdnn = 'edgeHiddenLEDinPupil-ZS-2021-05-19'
cdnn = 'pupilAllen-ZS-2023-12-20'
#cdnn = 'pupilEdge-Katie-2021-12-08'

df = pd.read_excel(file, sheet_name = sheet, usecols = [col])

for i in range(len(df)):
    deeplabcut.extractPupil(df[col][i], cdnn, label = False)
           
    
#%%
ani = ['ZS082', 'ZS083', 'ZS084', 'ZS085', 'ZS086'];
dates = ['20220612'];

for i in range(len(ani)):
    for j in range(len(dates)):
        session = f"m{ani[i]}d{dates[j]}";  
        deeplabcut.extractPupil(session, cdnn, label = False)
    
    
    
#%%    
deeplabcut.add_new_videos(config,[videoNew],copy_videos=False)
#%%
file = r'F:\PSpupil.xlsx'
sheet  = 'PS02pupil'
col = 'all'
cdnn = 'tongueTrackingStraight-ZS-2022-12-01'

df = pd.read_excel(file, sheet_name = sheet, usecols = [col])

for i in range(len(df)):
    if not pd.isna(df[col][i]):
        deeplabcut.extractLick(df[col][i], cdnn, label = False)
#%%    
           



