import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from pathlib import Path
import argparse
from scipy.spatial import distance
from math import factorial, atan2, degrees, acos, sqrt, pi

def currComputer():
    a = 'Z:\\'
    return a


#%%
def extractPupil(session, cdnn, label = False):        
    import deeplabcut
    root = currComputer()
    video = [root + 'ZS040' + '\\' + session + '\\' + 'pupil' + '\\' + session + '.avi']
    config_path = 'C:\\Users\\zhixi\\Documents\\dlc\\' + cdnn + '\\' + 'config.yaml'
    print('processing ' + session)
    deeplabcut.analyze_videos(config_path, video, save_as_csv=True)
    if label:
        deeplabcut.create_labeled_video(config_path,video, save_frames=False)
    
    deeplabcut.analyzeskeleton(config_path, video, videotype='avi', shuffle=1, trainingsetindex=0, save_as_csv=True)    
    
#%%