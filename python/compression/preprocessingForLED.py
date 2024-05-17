# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 17:04:36 2024

@author: zhixi
"""
import calMean
import numpy as np

if __name__ == '__main__': 
    
    files = [r"Z:\Xinxin\TestVideoRecordings\ledTest\behavior_0_2024-04-05_11-43-09",
            r"Z:\Xinxin\TestVideoRecordings\ledTest\behavior_0_2024-04-05_11-50-43",
            r"Z:\Xinxin\TestVideoRecordings\ledTest\behavior_0_2024-04-05_11-54-27"]
    sessions = ['behavior_0_2024-04-05_11-43-09',
            'behavior_0_2024-04-05_11-50-43',
            'behavior_0_2024-04-05_11-54-27']
    cameras = ['body_camera', 'bottom_camera', 'side_camera_left', 'side_camera_right']
    #%%
    videoDir = r'F:\lickSampleVidoes\ledTest'
    for session in sessions:
        for video in cameras:
            video_path = f'{videoDir}\{session}\\behavior-videos\{video}.avi'
            roiFile =  f'{videoDir}\{session}\\behavior-videos\{video}.txt'
            savingFile = f'{videoDir}\{session}\\behavior-videos\{video}.npy'
            
            roi =  calMean.read_numbers_from_file(roiFile)
            mean_brightness_list = calMean.calculate_mean_brightness_parallel(video_path, roi)
            print(f'Saving {session} {video}')
            np.save(savingFile, mean_brightness_list)
            print(f'Saved {session} {video}')

#%%
# load ROI

    

# Save list to a file in npy format

           
           