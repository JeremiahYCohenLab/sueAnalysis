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
from pynwb import NWBFile, TimeSeries, NWBHDF5IO
import json
import glob
from aind_ephys_utils import align
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
plt.text(-0.5, 5, 'Harp')                        
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
#%% led on times
session = 'behavior_0_2024-04-05_11-54-27'
videoDir = r'F:\lickSampleVidoes\ledTest'
video = 'side_camera_left'
video_path = f'{videoDir}\{session}\\behavior-videos\{video}.avi'
sigFile = f'{videoDir}\{session}\\behavior-videos\{video}.npy'
roiFile =  f'{videoDir}\{session}\\behavior-videos\{video}.txt'

#%% load signal 
csvFile = f'{videoDir}\{session}\\behavior-videos\{video}.csv'
frameTimes = pd.read_csv(csvFile, header = None)
frameTimesHarp = frameTimes[0] - frameTimes[0][0]
sig = np.load(sigFile)
plt.plot(sig)
#%% detect uprise
thresh = 25
sigUp = (sig[1:]>thresh) & (sig[0:-1]<thresh)
sigUp = np.array((np.where(sigUp))) + 1
sigUpSig = sig[sigUp]
sigPreSig = sig[sigUp-1]

#%% detect led On

binEdges = np.arange(np.min(np.diff(sig)), np.max(np.diff(sig)), 0.1)
plt.clf()
plt.hist(np.diff(sig), bins= binEdges, density=True)
plt.hist(sigUpSig.T-sigPreSig.T, bins= binEdges, density=True)

plt.show()


#%%
json_file = glob.glob(f'{videoDir}\{session}\\behavior\*.json')[0]
behavior_data = json.load(open(json_file))
endTime = behavior_data['B_TrialEndTimeHarp'] - frameTimes[0][0]
goCueTime = behavior_data['B_GoCueTimeBehaviorBoard'] - frameTimes[0][0]
delayStart = behavior_data['B_DelayStartTimeHarp'] - frameTimes[0][0]
allTrigger = np.sort(np.concatenate((endTime, goCueTime, delayStart)))

timeUp = frameTimesHarp[sigUp.flatten()].values

df = align.to_events(timeUp, allTrigger, (-0.01, 0.01), return_df=True)
plt.clf()
plt.scatter(1000*df.time, df.event_index, c='k', s=0.5)
# plt.plot([0,0],[0,df.event_index.max()],'r')
plt.xlabel('Frame time - Event time (ms)')
plt.ylabel('Event No.')
plt.xlim(-5, 5)
plt.show()
#%%
plt.plot(frameTimesHarp, sig, c = 'k')
plt.scatter(endTime, 15.1*np.ones_like(endTime)+15)
plt.scatter(goCueTime, 14.8*np.ones_like(goCueTime)+15)
plt.scatter(delayStart, 14.5*np.ones_like(delayStart)+15)
plt.xlim(endTime[5],endTime[10])
plt.show()
#%%
jitter = df.time.values
plt.scatter(allTrigger, jitter)
plt.show()
#%%
df = align.to_events(timeUp, endTime, (-0.01, 0.01), return_df=True)
jitter = df.time.values
plt.scatter(endTime, 1000*jitter)

df = align.to_events(timeUp, goCueTime, (-0.01, 0.01), return_df=True)
jitter = df.time.values
plt.scatter(goCueTime, 1000*jitter)

df = align.to_events(timeUp, delayStart, (-0.01, 0.01), return_df=True)
jitter = df.time.values
plt.scatter(delayStart, 1000*jitter)

plt.legend(['end', 'goCue', 'delayStart'])
plt.xlabel('Event time in sessions(s)')
plt.ylabel('Frame time - event time (ms)')
plt.show()
#%%
 


















