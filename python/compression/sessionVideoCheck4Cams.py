# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:03:16 2024

@author: zhixi
"""
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
# four views before changing codec
files = [r"Z:\Xinxin\FourCameras\behavior_0_2024-04-04_18-39-28",
        r"Z:\Xinxin\FourCameras\behavior_706893_2024-04-04_13-35-57",
        r"Z:\Xinxin\FourCameras\behavior_713379_2024-04-04_14-54-57",
        r"Z:\Xinxin\FourCameras\behavior_715083_2024-04-04_15-56-40",
        r"Z:\Xinxin\FourCameras\behavior_722679_2024-04-04_16-58-05",
        ]
#%% after codec correction
files = [r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717375\behavior_717375_2024-04-08_10-25-03",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717375\behavior_717375_2024-04-09_09-57-13",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717375\behavior_717375_2024-04-10_07-59-22",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717375\behavior_717375_2024-04-11_09-12-59",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717375\behavior_717375_2024-04-12_08-23-47",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717377\behavior_717377_2024-04-08_11-54-57",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717377\behavior_717377_2024-04-09_11-08-51",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\717377\behavior_717377_2024-04-12_10-04-45",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\722679\behavior_722679_2024-04-10_15-13-53",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\722679\behavior_722679_2024-04-11_16-52-29",
         r"Z:\ephys_rig_behavior_transfer\323_EPHYS3\722679\behavior_722679_2024-04-12_14-01-16"]
#%% led test
files = [r"Z:\Xinxin\TestVideoRecordings\ledTest\behavior_0_2024-04-05_11-43-09",
        r"Z:\Xinxin\TestVideoRecordings\ledTest\behavior_0_2024-04-05_11-50-43",
        r"Z:\Xinxin\TestVideoRecordings\ledTest\behavior_0_2024-04-05_11-54-27"]
#%% count frames from camera and csv
cameraList = ['body_camera', 'bottom_camera', ]
frameCountBottom = np.ones([len(files),1])
frameCountSideLeft = np.ones([len(files),1])
frameCountSideRight = np.ones([len(files),1])
frameCountBody = np.ones([len(files),1])

csvTimesBottom = [None]*len(files)
csvTimesSideLeft = [None]*len(files)
csvTimesSideRight = [None]*len(files)
csvTimesBody = [None]*len(files)
for i, file in enumerate(files): 
    # frames from videos
    videoBottom = path.join(file, r'behavior-videos\bottom_camera.avi')  
    videoSideLeft = path.join(file, r'behavior-videos\side_camera_left.avi')
    videoSideRight = path.join(file, r'behavior-videos\side_camera_right.avi')
    videoBody = path.join(file, r'behavior-videos\body_camera.avi')
    
    capBottom = cv2.VideoCapture(videoBottom)
    capSideLeft = cv2.VideoCapture(videoSideLeft)
    capSideRight = cv2.VideoCapture(videoSideRight)
    capBody = cv2.VideoCapture(videoBody)
    
    frameCountBottom[i] = float(capBottom.get(cv2.CAP_PROP_FRAME_COUNT))
    frameCountSideLeft[i] = float(capSideLeft.get(cv2.CAP_PROP_FRAME_COUNT))
    frameCountSideRight[i] = float(capSideRight.get(cv2.CAP_PROP_FRAME_COUNT))
    frameCountBody[i] = float(capBody.get(cv2.CAP_PROP_FRAME_COUNT))
      
    # frames from csv
    csvBottom = path.join(file, r'behavior-videos\bottom_camera.csv')  
    csvSideLeft = path.join(file, r'behavior-videos\side_camera_left.csv')
    csvSideRight = path.join(file, r'behavior-videos\side_camera_right.csv')
    csvBody = path.join(file, r'behavior-videos\body_camera.csv')
    
    csvTimesBottom[i] = pd.read_csv(csvBottom, header = None)
    csvTimesSideLeft[i] = pd.read_csv(csvSideLeft,  header = None)
    csvTimesSideRight[i] = pd.read_csv(csvSideRight,  header = None)
    csvTimesBody[i] = pd.read_csv(csvBody,  header = None)
    
csvLensBottom = [len(a) for a in csvTimesBottom]
csvLensSideLeft = [len(a) for a in csvTimesSideLeft]   
csvLensSideRight = [len(a) for a in csvTimesSideRight] 
csvLensBody = [len(a) for a in csvTimesBody]   

 # trigger counts in harp
# metadata = pd.read_csv(root / r"VideoFolder\side_camera.csv", header=None)
triggers = [None]*len(files)
for i, file in enumerate(files):
    harpfile = path.join(file, r'behavior\raw.harp\BehaviorEvents\Event_94.bin')
    tmp = harp.read(harpfile)
    triggers[i] = tmp.index
    
triggerLens = [len(a) for a in triggers]
#
allFrameCounts = pd.DataFrame({'harpLen': triggerLens, 
                               'frameCountBottom': frameCountBottom.flatten(), 'csvLensBottom': csvLensBottom, 
                               'frameCountSideLeft': frameCountSideLeft.flatten(), 'csvLensSideLeft': csvLensSideLeft,
                               'frameCountSideRight': frameCountSideRight.flatten(), 'csvLensSideRight': csvLensSideRight,
                               'frameCountBody': frameCountBody.flatten(), 'csvLensBody': csvLensBody})
#%% summary frame dropping
allFrameRatio = allFrameCounts.values
triggerNum = allFrameCounts['harpLen']
triggerNum = triggerNum[:, np.newaxis]
allFrameRatio = allFrameRatio/triggerNum
allFrameDrop = allFrameCounts.values == triggerNum

names = ['bottom', 'sideL', 'sideR', 'body']

plt.imshow(allFrameRatio, cmap = 'Reds_r')
vertical_lines_x = np.arange(0.5, 8.5, 2)
vertical_lines_x_sub = np.arange(1.5, 8.5, 2)
vertical_lines_y = np.arange(0.5, 4.5, 1)
# Plot vertical lines
plt.vlines(vertical_lines_x, ymin=-0.5, ymax=len(allFrameCounts)-0.5, color='k', linestyle='--', lw = 0.3)
plt.vlines(vertical_lines_x_sub, ymin=-0.5, ymax=len(allFrameCounts)-0.5, color='k', linestyle='--', lw = 0.05)
plt.hlines(vertical_lines_y, xmin=-0.5, xmax=8.5, color='k', linestyle='--', lw = 0.3)
plt.xticks(np.array(2*np.arange(0, len(names)) + 1.5), names)
plt.text(-0.5, len(files)+0.25, 'Harp')                        
plt.ylabel('session')
plt.colorbar()  # Add color bar for reference
#%%
plt.imshow(allFrameDrop, cmap = 'Reds_r')
vertical_lines_x = np.arange(0.5, 8.5, 2)
vertical_lines_x_sub = np.arange(1.5, 8.5, 2)
vertical_lines_y = np.arange(0.5, 4.5, 1)
# Plot vertical lines
plt.vlines(vertical_lines_x, ymin=-0.5, ymax=len(allFrameCounts)-0.5, color='k', linestyle='--', lw = 0.3)
plt.vlines(vertical_lines_x_sub, ymin=-0.5, ymax=len(allFrameCounts)-0.5, color='k', linestyle='--', lw = 0.05)
plt.hlines(vertical_lines_y, xmin=-0.5, xmax=8.5, color='k', linestyle='--', lw = 0.3)
plt.xticks(np.array(2*np.arange(0, len(names)) + 1.5), names)
plt.text(-0.5, 5, 'Harp')                        
plt.ylabel('session')
plt.colorbar()
#%% analyze frame dropping over time
plt.clf()
fig, axs = plt.subplots(len(files), 1, figsize=(12, 6))
for ind, file in enumerate(files):
    currCsv = csvTimesSideLeft[ind][0]
    bins = np.arange(np.min(triggers[ind])-0.01, np.max(triggers[ind]) +0.01, 30)
    axs[ind].hist(triggers[ind], bins = bins, alpha = 0.2, color = 'k')
    bins = np.arange(np.min(currCsv)-0.01, np.max(currCsv) +0.01, 30)
    axs[ind].hist((currCsv), bins = bins, alpha = 0.2)

plt.show()
#%%
plt.clf()
fig, axs = plt.subplots(1, len(files), figsize=(12, 6))
for ind, file in enumerate(files):
    currCsv = csvTimesBottom[ind][1]
    axs[ind].hist((np.diff(currCsv)),alpha = 0.2)
    # axs[i].title(ind)
    numInt = len(np.unique(np.diff(currCsv)))
    axs[ind].set_title(numInt)

plt.suptitle('SideL')
plt.show()

#%% plot for each session
b = 0.002
a = 0.02
# frame time from bottom and side camera so that we can combine them as one
for i, frameTimesCurr in enumerate(csvTimesBody):
    
    x = np.diff(frameTimesCurr[0])
    y = np.diff(frameTimesCurr[2])
    # xbins = np.arange(1000*min(x), 1000*max(x)+a, a)
    # ybins = np.arange(min(y)/1000000, (max(y)/1000000)+b, b)
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


