# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 18:34:23 2024

@author: zhixi
"""

#%%

import os, ffmpeg
from pathlib import Path
import subprocess


def compress_video(input_file, output_file, crf=23, bitrate = '50M'):
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
input_file =  r"F:\lickSampleVidoes\ephysRig\testBottom-01242024233257-0000.avi"
output_file = r"F:\lickSampleVidoes\ephysRig\testBottom_test_5.avi"
compress_video(input_file, output_file, crf=51, bitrate='100M')
