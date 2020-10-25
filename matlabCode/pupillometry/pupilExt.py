# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
def currComputer():
    a = 'Z:\\'
    return a


#%%
def extractPupil(session, cdnn, label = False):
    import sys
    if 'deeplabcut' not in sys.modules:
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