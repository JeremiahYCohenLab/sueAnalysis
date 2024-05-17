# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:13:36 2024

@author: zhixi
"""
#%%
import cv2
import numpy as np
import matplotlib.pyplot as plt
from os import path

import pandas as pd
import seaborn as sns
import harp
import csv
#%%
# file = input("Enter the path to the MP4 video file: ")
# four views
files = [r"Z:\Xinxin\FourCameras\behavior_0_2024-04-04_18-39-28",
        r"Z:\Xinxin\FourCameras\behavior_706893_2024-04-04_13-35-57",
        r"Z:\Xinxin\FourCameras\behavior_713379_2024-04-04_14-54-57",
        r"Z:\Xinxin\FourCameras\behavior_715083_2024-04-04_15-56-40",
        r"Z:\Xinxin\FourCameras\behavior_722679_2024-04-04_16-58-05"]
#%% four views after codec change

#%%
# file = input("Enter the path to the MP4 video file: ")
# two views
files = [r"Z:\Xinxin\TestVideoRecordings\behavior_722679_2024-03-29_16-16-12",
        r"Z:\Xinxin\TestVideoRecordings\behavior_706893_2024-03-29_15-06-43",
        r"Z:\Xinxin\TestVideoRecordings\behavior_722679_2024-04-01_16-19-23 - Copy",
        r"Z:\Xinxin\TestVideoRecordings\behavior_706893_2024-04-01_15-03-58",
        r"Z:\Xinxin\TestVideoRecordings\behavior_0_2024-04-01_12-35-16"]
#%% count frames
frameCountBottom = np.ones([len(files),1])
frameCountSide = np.ones([len(files),1])

for i, file in enumerate(files): 
    videoBottom = path.join(file, r'behavior-videos\bottom_camera.avi')  
    videoSide = path.join(file, r'behavior-videos\side_camera.avi')
    
    capBottom = cv2.VideoCapture(videoBottom) 
    capSide = cv2.VideoCapture(videoSide)
    
    frameCountBottom[i] = float(capBottom.get(cv2.CAP_PROP_FRAME_COUNT))
    frameCountSide[i] = float(capSide.get(cv2.CAP_PROP_FRAME_COUNT))
#%% tigger frame time and camera frame time
frameTimesBottom = [None]*len(files)
frameTimesSide = [None]*len(files)
for i, file in enumerate(files): 
    videoBottom = path.join(file, r'behavior-videos\bottom_camera.csv')  
    videoSide = path.join(file, r'behavior-videos\side_camera.csv')
    
    frameTimesBottom[i] = pd.read_csv(videoBottom, header = None)
    frameTimesSide[i] = pd.read_csv(videoSide,  header = None)

frameLensBottom = [len(a) for a in frameTimesBottom]
frameLensSide = [len(a) for a in frameTimesSide]
 #%% trigger counts
# metadata = pd.read_csv(root / r"VideoFolder\side_camera.csv", header=None)
triggers = [None]*len(files)
for i, file in enumerate(files):
    harpfile = path.join(file, r'behavior\raw.harp\BehaviorEvents\Event_94.bin')
    tmp = harp.read(harpfile)
    triggers[i] = tmp.index
    
triggerLens = [len(a) for a in triggers]
#%% plot for each session
# frame time from bottom and side camera so that we can combine them as one
for i, frameTimesCurr in enumerate(frameTimesBottom):
    x = np.diff(frameTimesCurr[0])
    y = np.diff(frameTimesCurr[2])
    sns.jointplot(x = 1000*x, y = y/1000000, kind = 'hist', bins = 10,
                       stat='probability',marginal_ticks=True, 
                       marginal_kws=dict(bins=20, fill=True, stat='probability')).set_axis_labels(xlabel = 'harpTime', ylabel='frameTime')
    plt.suptitle(str(i))
    plt.show()
#%% summary
plt.scatter(frameCountBottom, frameCountSide)
plt.xlabel('bottomCameraCount')
plt.ylabel('sideCameraCount')
plt.show()
#%% time fluctuations
for i, frameTimesCurr in enumerate(frameTimesBottom):
    x = np.diff(frameTimesBottom[i][0])*1000
    y = np.diff(frameTimesBottom[i][2])/1000000
    z = np.diff(frameTimesSide[i][2])/1000000
    frameNum = range(int(frameCountBottom[i])-1)
    fig, axs = plt.subplots(3, 1, figsize=(15, 8))
    
    axs[0].set_xlim(0, frameCountBottom[i])
    # axs[0].set_ylim(-5,5)
    axs[0].plot(frameNum, x, lw = 0.1, color = 'red')
    axs[0].set_title('trigger time')

    
    axs[1].set_xlim(0, frameCountBottom[i])
    # axs[1].set_ylim(-5,5)
    axs[1].plot(frameNum, y, lw = 0.1, color = 'blue')
    axs[1].set_title('frame time Bottom')
    
    axs[2].set_xlim(0, frameCountBottom[i])
    # axs[1].set_ylim(-5,5)
    axs[2].plot(frameNum, z, lw = 0.1, color = 'blue')
    axs[2].set_title('frame time Side')
    plt.show()
    #%% scatter distribution
startNum = 10000
endNum = 10000
alpha = 0.2
for i in range(4):
    x = np.diff(frameTimesBottom[i][0])*1000
    y = np.diff(frameTimesBottom[i][2])/1000000
    z = np.diff(frameTimesSide[i][2])/1000000
    frameNum = range(int(frameCountBottom[i])-1)
    fig, axs = plt.subplots(3, 2, figsize=(15, 8))
    
    
    # axs[0, 0].set_xlim(0, frameCountBottom[i])
    axs[0, 0].set_ylim(np.min(x)-0.0001, np.max(x)+0.0001)
    axs[0, 0].scatter(frameNum[:startNum], x[:startNum], s=3, alpha=alpha)
    axs[0, 0].set_title('trigger time Start')
    axs[0, 0].tick_params(axis='x', labelbottom=False)
      
    axs[0, 1].set_ylim(np.min(x)-0.0001, np.max(x)+0.0001)
    # axs[0].set_ylim(-5,5)
    axs[0, 1].scatter(frameNum[-endNum:], x[-endNum:], s=3, alpha=alpha)
    axs[0, 1].set_title('trigger time End')
    axs[0, 1].tick_params(axis='x', labelbottom=False)


    
    # axs[1, 0].set_xlim(0, frameCountBottom[i])
    # axs[1].set_ylim(-5,5)
    axs[1, 0].set_ylim(np.min(y)-0.0001, np.max(y)+0.0001)
    axs[1, 0].scatter(frameNum[:startNum], y[:startNum], s=3, alpha=alpha)
    axs[1, 0].set_title('frame time Bottom Start')
    axs[1, 0].tick_params(axis='x', labelbottom=False)
    
    # axs[1, 1].set_xlim(0, frameCountBottom[i])
    # axs[1].set_ylim(-5,5)
    axs[1, 1].set_ylim(np.min(y)-0.0001, np.max(y)+0.0001)
    axs[1, 1].scatter(frameNum[-endNum:], y[-endNum:], s=3, alpha=alpha)
    axs[1, 1].set_title('frame time Bottom End')
    axs[1, 1].tick_params(axis='x', labelbottom=False)
    
    
    
    # axs[2, 0].set_xlim(0, frameCountBottom[i])
    # axs[1].set_ylim(-5,5)
    axs[2, 0].set_ylim(np.min(z)-0.0001, np.max(z)+0.0001)
    axs[2, 0].scatter(frameNum[:startNum], z[:startNum], s=3, alpha=alpha)
    axs[2, 0].set_title('frame time Side Start')
    
    # axs[2, 1].set_xlim(0, frameCountBottom[i])
    # axs[1.set_ylim(-5,5)
    axs[2, 1].set_ylim(np.min(z)-0.0001, np.max(z)+0.0001)
    axs[2, 1].scatter(frameNum[-endNum:], z[-endNum:], s=3, alpha=alpha)
    axs[2, 1].set_title('frame time Side End')
    
    plt.suptitle(str(i))
    plt.show()

  #%%


