# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 18:07:27 2024

@author: zhixi

"""
import cv2
import multiprocessing
import numpy as np
import json
import cupy as cp


def calculate_mean_brightness(frame, roi):
    """Function to calculate the mean brightness of a frame."""
    mean_brightness = np.mean(frame[roi[0]:(roi[0]+roi[2]), roi[1]:(roi[1]+roi[3])])
    return mean_brightness



def calculate_mean_brightness_parallel(video_path, roi):
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print("Error: Couldn't open the video file.")
        return None
    
    mean_brightness_list = []

    while True:
        ret, frame = cap.read()
        if not ret:
            break
        if len(mean_brightness_list)%5000 == 0:
            print(len(mean_brightness_list))
        mean_brightness = calculate_mean_brightness(frame, roi)
        mean_brightness_list.append(mean_brightness)
        
            
            
        # if len(mean_brightness_list) >= 30000:
        #     break

    cap.release()
    mean_brightness_list = np.array(mean_brightness_list)
    return mean_brightness_list

def read_numbers_from_file(filename):
    """
    Function to split each line of a text file by '=' and ',' and extract numbers.

    Args:
        filename: The filename of the text file.
    
    Returns:
        A list of numbers extracted from the file.
    """
    numbers = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split(',')  # Split line by '='
            for part in parts:
                elements = part.split('=')  # Split each part by ','
                for element in elements:
                    try:
                        number = int(element)  # Convert the element to a float number
                        numbers.append(number)
                    except ValueError:
                        print(f"Ignoring non-numeric value: '{element}'.")
    return numbers


# Example usage:
if __name__ == '__main__':
  
    session = 'behavior_0_2024-04-05_11-54-27'
    videoDir = r'F:\lickSampleVidoes\ledTest'
    video = 'body_camera'
    video_path = f'{videoDir}\{session}\\behavior-videos\{video}.avi'
    savingFile = f'{videoDir}\{session}\\behavior-videos\{video}.npy'
    roiFile =  f'{videoDir}\{session}\\behavior-videos\{video}.txt'
    # load ROI
    roi =  read_numbers_from_file(roiFile)
        
    mean_brightness_list = calculate_mean_brightness_parallel(video_path, roi)
    # Save list to a file in npy format
    print('Saving')
    np.save(savingFile, mean_brightness_list)
    print('Saved')
