# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:34:23 2024

@author: zhixi
"""

#%%

import os, ffmpeg
from pathlib import Path
import subprocess

input_file =  r"C:\Users\zhixi\OneDrive - Allen Institute\testCamera\exampleTwoSideVideos\tongueShapeChangeForWater.avi"
output_file = r"C:\Users\zhixi\OneDrive - Allen Institute\testCamera\exampleTwoSideVideos\tongueShapeChangeForWatercompressed.avi"
# def compress_video(video_full_path, output_file_name, target_size):
#     # Reference: https://en.wikipedia.org/wiki/Bit_rate#Encoding_bit_rate
#     min_audio_bitrate = 32000
#     max_audio_bitrate = 256000

#     probe = ffmpeg.probe(video_full_path)
#     # Video duration, in s.
#     duration = float(probe['format']['duration'])
#     # Audio bitrate, in bps.
#     # audio_bitrate = float(next((s for s in probe['streams'] if s['codec_type'] == 'audio'), None)['bit_rate'])
#     # Target total bitrate, in bps.
#     target_total_bitrate = (target_size * 1024 * 8) / (1.073741824 * duration)

#     # Target audio bitrate, in bps
#     # if 10 * audio_bitrate > target_total_bitrate:
#     #     audio_bitrate = target_total_bitrate / 10
#     #     if audio_bitrate < min_audio_bitrate < target_total_bitrate:
#     #         audio_bitrate = min_audio_bitrate
#     #     elif audio_bitrate > max_audio_bitrate:
#     #         audio_bitrate = max_audio_bitrate
#     # Target video bitrate, in bps.
#     audio_bitrate = target_total_bitrate/10;
#     video_bitrate = target_total_bitrate - audio_bitrate

#     i = ffmpeg.input(video_full_path)
#     ffmpeg.output(i, os.devnull,
#                   **{'c:v': 'libx264', 'b:v': video_bitrate, 'pass': 1, 'f': 'avi'}
#                   ).overwrite_output().run()
#     ffmpeg.output(i, output_file_name,
#                   **{'c:v': 'libx264', 'b:v': video_bitrate, 'pass': 2, 'c:a': 'aac', 'b:a': audio_bitrate}
#                   ).overwrite_output().run()

# # Compress input.mp4 to 50MB and save as output.mp4
# # compress_video(input_file, output_file, 100 * 1000)
#%%


def compress_video(input_file, output_file, crf=1, bitrate = '1M'):
    """
    Compresses a video file using ffmpeg.
    
    Parameters:
        input_file (str): Path to the input video file.
        output_file (str): Path to the output compressed video file.
        crf (int): Constant Rate Factor for video compression (default is 23).
    """
    if os.path.exists(output_file):
        os.remove(output_file)  # Remove existing output file if it exists
    
    cmd = [
        'ffmpeg',
        '-i', input_file,
        '-c:v', 'h264_nvenc',
        '-b:v', bitrate,
        '-crf', str(crf),
        '-preset', 'slow',
        '-an',
        output_file
    ]
    
    subprocess.run(cmd)
    
# Example usage:
compress_video(input_file, output_file, crf = 52, bitrate= '3M')

#%%
# import cv2
# import numpy as np
# import matplotlib.pyplot as plt
# # file = input("Enter the path to the MP4 video file: ")
# file = "F:\\lickSampleVidoes\\ephysRig\\12-17-43.243.avi"
# cap = cv2.VideoCapture(file) 

# # cv2.namedWindow("Video Player", cv2.WINDOW_NORMAL)
  
# while(cap.isOpened()):
#     success, frame = cap.read()
#     # frame = cv2.resize(frame, (200, 150))
#     if success == True:
#         gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
#         gray = cv2.normalize(gray, None, 0, 256, cv2.NORM_MINMAX)
#         cv2.imshow('Video', gray)
#         cv2.waitKey(1)
#         # plt.imshow(gray);
#     else:
#         break
    
    
       
# cap.release()
# cv2.destroyAllWindows()

#%%
