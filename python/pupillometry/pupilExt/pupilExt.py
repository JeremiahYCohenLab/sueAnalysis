

def currComputer():
    a = 'Z:\\'
    return a


#%%
def extractPupil(session, cdnn, label = False):        
    import deeplabcut
    root = deeplabcut.currComputer()
    ani = session.split('d')[0];
    ani = ani.split('m')[1];
    video = [root + ani + '\\' + session + '\\' + 'pupil' + '\\' + session + '.avi']
    config_path = 'C:\\Users\\zhixi\\Documents\\dlc\\' + cdnn + '\\' + 'config.yaml'
    print('processing ' + session)
    deeplabcut.analyze_videos(config_path, video, save_as_csv=True)
    if label:
        deeplabcut.create_labeled_video(config_path,video, save_frames=False, shuffle=2)
    
    deeplabcut.analyzeskeleton(config_path, video, videotype='avi', shuffle=1, trainingsetindex=0, save_as_csv=True)    
    
#%%