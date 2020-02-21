# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
def my_function():
    print("Hello From My Function!")


#%%
def extractPupil(session):
    print('processing ' + session)
    
    
    
    
 #%%
import deeplabcut
config_path = r'C:\Users\zhixi\Documents\dlc\rwdDelay-combine-2020-02-20\config.yaml'
videos = [r'C:\Users\zhixi\Documents\dlc\mZS040d20200217\pupil\mZS040d20200217.avi']
deeplabcut.analyzeskeleton(config_path, videos, videotype='avi', shuffle=1, trainingsetindex=0, save_as_csv=True)