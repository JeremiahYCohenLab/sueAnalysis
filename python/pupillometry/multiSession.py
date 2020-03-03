# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 00:25:39 2020

@author: zhixi
"""
file = r'Z:\combineRwdDelay.xlsx'
sheet  = 'combine'
col = 'OK'
cdnn =  cdnn = 'rwdDelay-combine-2020-02-20'


from pandas import DataFrame, read_csv
import matplotlib.pyplot as plt
import pandas as pd
import deeplabcut


df = pd.read_excel(file, sheet_name = sheet, usecols = [col])


for i in range(0, len(df)-1):
    session = df[col][i]
    deeplabcut.extractPupil(session, cdnn)
    
#%%                