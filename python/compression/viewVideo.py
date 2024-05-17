# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:34:23 2024

@author: zhixi
"""

# #%%

# import os, ffmpeg
# from pathlib import Path

# input_file = file;
# output_file = Path('F:\lickSampleVidoes\ephysRig/testBottom-01242024233257-0000_test.avi');
# def compress_video(video_full_path, output_file_name, target_size):
#     # Reference: https://en.wikipedia.org/wiki/Bit_rate#Encoding_bit_rate
#     min_audio_bitrate = 32000
#     max_audio_bitrate = 256000

#     probe = ffmpeg.probe(video_full_path)
#     # Video duration, in s.
#     duration = float(probe['format']['duration'])
#     # Audio bitrate, in bps.
#     audio_bitrate = float(next((s for s in probe['streams'] if s['codec_type'] == 'audio'), None)['bit_rate'])
#     # Target total bitrate, in bps.
#     target_total_bitrate = (target_size * 1024 * 8) / (1.073741824 * duration)

#     # Target audio bitrate, in bps
#     if 10 * audio_bitrate > target_total_bitrate:
#         audio_bitrate = target_total_bitrate / 10
#         if audio_bitrate < min_audio_bitrate < target_total_bitrate:
#             audio_bitrate = min_audio_bitrate
#         elif audio_bitrate > max_audio_bitrate:
#             audio_bitrate = max_audio_bitrate
#     # Target video bitrate, in bps.
#     video_bitrate = target_total_bitrate - audio_bitrate

#     i = ffmpeg.input(video_full_path)
#     ffmpeg.output(i, os.devnull,
#                   **{'c:v': 'libx264', 'b:v': video_bitrate, 'pass': 1, 'f': 'mp4'}
#                   ).overwrite_output().run()
#     ffmpeg.output(i, output_file_name,
#                   **{'c:v': 'libx264', 'b:v': video_bitrate, 'pass': 2, 'c:a': 'aac', 'b:a': audio_bitrate}
#                   ).overwrite_output().run()

# # Compress input.mp4 to 50MB and save as output.mp4
# compress_video(input_file, output_file, 500 * 1000)
#%%
import cv2
import numpy as np
import matplotlib.pyplot as plt
# file = input("Enter the path to the MP4 video file: ")
files = ['r"Z:\Xinxin\TestVideoRecordings\behavior_722679_2024-03-29_16-16-12\behavior-videos\bottom_camera.avi"',
        'r"Z:\Xinxin\TestVideoRecordings\behavior_722679_2024-03-29_16-16-12\behavior-videos\bottom_camera.avi"',
        'r"Z:\Xinxin\TestVideoRecordings\behavior_722679_2024-03-29_16-16-12\behavior-videos\bottom_camera.avi"']

for cap = cv2.VideoCapture(file) 

# cv2.namedWindow("Video Player", cv2.WINDOW_NORMAL)
  
while(cap.isOpened()):
    success, frame = cap.read()
    # frame = cv2.resize(frame, (200, 150))
    if success == True:
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        gray = cv2.normalize(gray, None, 0, 256, cv2.NORM_MINMAX)
        cv2.imshow('Video', gray)
        cv2.waitKey(1)
        # plt.imshow(gray);
    else:
        break
    
    
     
cap.release()
cv2.destroyAllWindows()
#%%
# file = "F:\\lickSampleVidoes\\ephysRig\\testBottom-01242024233257-0000.avi";
# import cv2

# from ffpyplayer.player import MediaPlayer
# video=cv2.VideoCapture(file)
# player = MediaPlayer(file)
# while True:
#    ret, frame=video.read()
#    audio_frame, val = player.get_frame()
#    if not ret:
#       print("End of video")
#       break
#    if cv2.waitKey(1) == ord("q"):
#       break
#    cv2.imshow("Video", frame)
#    if val != 'eof' and audio_frame is not None:
#       #audio
#       img, t = audio_frame
# video.release()
# cv2.destroyAllWindows()
#%%
# length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
# print( length )