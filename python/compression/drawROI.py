# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 17:04:33 2024

@author: zhixi
"""
import cv2

# Global variables to track mouse events and ROI coordinates
roi_start = None
roi_end = None
roi_selected = False

def draw_roi(image, roi):
    """
    Function to draw a region of interest (ROI) on an image.

    Args:
        image: The input image.
        roi: A tuple (x, y, w, h) representing the coordinates and size of the ROI rectangle.
    Returns:
        An image with the ROI drawn.
    """
    x, y, w, h = roi
    roi_image = image.copy()
    cv2.rectangle(roi_image, (x, y), (x + w, y + h), (0, 255, 0), 2)  # Draw rectangle
    return roi_image

def mouse_callback(event, x, y, flags, param):
    """
    Callback function for mouse events.

    Args:
        event: The type of mouse event (e.g., cv2.EVENT_LBUTTONDOWN, cv2.EVENT_MOUSEMOVE, etc.).
        x, y: The coordinates of the mouse cursor.
        flags: Any flags passed by OpenCV.
        param: Any extra parameters passed to the callback function.
    """
    global roi_start, roi_end, roi_selected

    if event == cv2.EVENT_LBUTTONDOWN:
        roi_start = (x, y)
        roi_end = (x, y)
        roi_selected = False
    elif event == cv2.EVENT_MOUSEMOVE:
        if roi_start is not None and not roi_selected:
            roi_end = (x, y)
    elif event == cv2.EVENT_LBUTTONUP:
        roi_end = (x, y)
        roi_selected = True
        
        
def import_first_frame(video_path):
    """
    Function to import the first frame of a video as an image.

    Args:
        video_path: The path to the video file.
    
    Returns:
        The first frame of the video as an image (numpy array).
    """
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print("Error: Couldn't open the video file.")
        return None

    # Read the first frame
    ret, frame = cap.read()
    if not ret:
        print("Error: Couldn't read the first frame.")
        return None

    # Release the video capture object
    cap.release()

    return frame

def save_roi_location(roi, filename):
    """
    Function to save the ROI location along with the label to a text file.

    Args:
        roi: A tuple (x, y, w, h) representing the coordinates and size of the ROI rectangle.
        label: The label associated with the ROI.
        filename: The filename of the text file to save the ROI location.
    """
    x, y, w, h = roi
    with open(filename, 'a') as f:
        f.write(f"x={x}, y={y}, width={w}, height={h}\n")


# Example usage:
if __name__ == "__main__":
    videoDir = r'F:\lickSampleVidoes\ledTest'
    session = 'behavior_0_2024-04-05_11-54-27';
    video = 'side_camera_right'
    video_path = f'{videoDir}\{session}\\behavior-videos\{video}.avi'
    savingLocation = f'{videoDir}\{session}\\behavior-videos\{video}.txt'
    image = import_first_frame(video_path)
    
    if image is None:
        print("Error: Could not read image.")
    else:
        cv2.namedWindow("ROI Image")
        cv2.setMouseCallback("ROI Image", mouse_callback)

        while True:
            roi_image = image.copy()
            if roi_start is not None and roi_end is not None:
                x, y = min(roi_start[0], roi_end[0]), min(roi_start[1], roi_end[1])
                w, h = abs(roi_start[0] - roi_end[0]), abs(roi_start[1] - roi_end[1])
                cv2.rectangle(roi_image, (x, y), (x + w, y + h), (0, 255, 0), 2)  # Draw rectangle
                
            if roi_selected:
                label = input("Enter label for ROI: ")
                save_roi_location((x, y, w, h), savingLocation)
                roi_selected = False
                
            cv2.imshow("ROI Image", roi_image)
            key = cv2.waitKey(1) & 0xFF
            if key == ord("q"):
                break

        cv2.destroyAllWindows()
        
