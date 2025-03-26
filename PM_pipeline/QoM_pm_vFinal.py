# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 13:32:45 2023

@author: fbigand
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 17:15:46 2022

@author: fbigand
"""

#%% 
##############################################################
############  IMPORT LIBRARIES AND SET PARAMETERS ############
##############################################################

from PLmocap.viz import *
from PLmocap.preprocessing import *
from PLmocap.classif import *
from PLmocap.stats import *
import mne_fefe
from ezc3d import c3d
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d import Axes3D, proj3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import signal, interpolate, sparse
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn import linear_model
from sklearn import metrics
import time
import numpy as np
import os
from pylab import *
import seaborn as sns
import pycwt as wavelet


######################## PARAMETERS ########################
fps = 25
DUR = 60
NB_MARKERS = 22
NB_DYADS = 35
NB_TRIALS = 32
NB_CONDITIONS = 4
NB_SONGS = 8
NB_METRICAL_LEVELS = 5

# Structure of the skeleton (markers + links between markers (for visualization))    
list_markers = np.array(['LB Head','LF Head','RF Head','RB Head','Sternum','L Shoulder','R Shoulder',\
                            'L Elbow','L Wrist','L Hand','R Elbow','R Wrist','R Hand','Pelvis','L Hip',\
                            'R Hip','L Knee','L Ankle','L Foot','R Knee','R Ankle','R Foot'])

liaisons = [(0,1),(0,3),(1,2),(2,3),(4,5),(4,6),(5,6),(5,7),(7,8),(8,9),(6,10),(10,11),(11,12), \
            (13,14),(13,15),(14,15),(14,16),(16,17),(17,18),(15,19),(19,20),(20,21)]
dim_name = ["x","y","z"]

# Information about the songs and their musical sections
song_bpms = np.array([111.03,116.07,118.23,118.95,120.47,125.93,128.27,129.06])
musParts_names = ['DRUMS','DRUMS+BASS','DRUMS+BASS+KEYBOARDS','DRUMS+BASS+KEYBOARDS+VOICE']
musParts_beats = [0,16,32,48,80] # start drums , bass, harmony, voice, end
beats_tFrames = np.zeros( (NB_SONGS,81) )
musParts_tFrames = np.zeros( (NB_SONGS,len(musParts_beats)) )
for i in range(NB_SONGS):
    periodbeat = (60/song_bpms[i])
    beats_tFrames[i,:] = np.linspace(0,80*periodbeat,81) # Because 80 beats for each song
    musParts_tFrames[i,:] = beats_tFrames[i,musParts_beats]
    

#%% 
##############################################################
############             IMPORT DATA              ############
##############################################################

print('------------------------')
print('LOADING DATA...')
print('------------------------')

# Input dir and output dir
input_dir = os.path.normpath(( os.getcwd() + "/vicondata/npy" ))
folders = os.listdir(input_dir);    folders=[x for i,x in enumerate(folders) if (x.startswith("Dyad"))]
folders=sorted(folders)

output_dir = os.path.normpath( os.getcwd() + "/QoM_commonPM_vFinal" )
if not (os.path.exists(output_dir)) : os.mkdir(output_dir)

# Create the concatenated multi-subject dataset matrix (all trials concatenated but 1 trial missing for Dyad 32)
LEN_CONCAT_TIME = (musParts_tFrames[:,-1]*fps).astype(int).sum()*2 + (musParts_tFrames[4:,-1]*fps).astype(int).sum()*4
data_subj = np.zeros((NB_DYADS*2 , NB_MARKERS*3 , LEN_CONCAT_TIME))

pmean_tr = np.zeros((NB_DYADS*2 , NB_MARKERS*3 , NB_TRIALS)); dmean_tr = np.zeros(( NB_DYADS*2 , NB_TRIALS))    # mean posture and normalization vector per trial

# Conditions and songs of every trial and dyad
cond_mus = np.zeros( (NB_DYADS , NB_TRIALS) , dtype=int ); cond_vis = np.zeros( (NB_DYADS , NB_TRIALS) , dtype=int ); 
song = np.zeros( (2 , NB_DYADS , NB_TRIALS) , dtype=int ); # first dimension is subj_left (0) and subj_right (1)


for d in range(1, NB_DYADS+1):
    if d < 10 : numDyad = "0" + str(d)
    else : numDyad = str(d)
    print("Dyad " + numDyad)
    
    # Motion data (from Vicon)
    folder = os.path.normpath(( os.getcwd() + "/vicondata/npy/" + folders[d-1] ))
    files=os.listdir(folder + '/subj_LEFT');      files=[x for i,x in enumerate(files) if (x.endswith(".npy"))]
    files=sorted(files)
    
    # Condition & song info (from Presentation)
    logfile = pd.read_excel( os.path.normpath(( os.getcwd() + "/Presentation_log/" + folders[d-1]  + ".xlsx")) , engine='openpyxl')
    logfile_sounds = logfile.loc[logfile.iloc[:,1]== "Sound"]
    log_conditions = logfile_sounds.iloc[np.where(logfile_sounds.iloc[:, 2].str.contains('Start_Trial') == True)[0], 2].to_numpy().astype(str)
    log_songs = logfile_sounds.iloc[np.where(logfile_sounds.iloc[:, 2].str.contains('Start_Trial') == True)[0] + 1, 2].to_numpy().astype(str)
    for tr in range(32):
        cond_mus[(d-1) , tr] = int('Same Music' in log_conditions[tr][12:])
        cond_vis[(d-1) , tr] = int('Yes Visual' in log_conditions[tr][12:])
        song[0,(d-1) , tr] = int(log_songs[tr][16]) - 1
        song[1,(d-1) , tr] = int(log_songs[tr][24]) - 1 
    
    
    # SUBJECT 1 AND 2
    for subj in range(2):   
        iStart=0
        for tr in range(NB_TRIALS):
            if tr+1 < 10 : numTrial = "0" + str(tr+1)
            else : numTrial = str(tr+1)  
            
            # Load motion data
            if subj == 0 : xyz_vec=np.load( os.path.normpath( folder + '/subj_LEFT/' + files[tr] ))
            if subj == 1 : xyz_vec=np.load( os.path.normpath( folder + '/subj_RIGHT/' + files[tr] )) 
            
            # Small last pre-processing step
            # Fill the data if NaN at start and/or at end (as Nexus only allows filling gaps if there is data before and after)
            # (for future: put it as a function of PL MOCAP)
            if False in np.isnan(xyz_vec[:]) : 
                for m in range(0,NB_MARKERS*3,3):
                    # CHECK IF FIRST FRAME IS MISSING
                    if np.isnan(xyz_vec[m,0]):
                        idx_first_val = np.min(np.where(~np.isnan(xyz_vec[m,:])))
                        
                        if m//3 in [0,3]: # if back head markers, do "pattern fill" with front head markers
                            if m//3==0: diff_donor = np.diff(xyz_vec[3:6,:], axis=1) # Left back
                            if m//3==3: diff_donor = np.diff(xyz_vec[6:9,:], axis=1) # Right back
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] + np.cumsum(diff_donor[0,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] + np.cumsum(diff_donor[1,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] + np.cumsum(diff_donor[2,:idx_first_val][::-1])[::-1]
                        
                        elif (subj==0) and (m//3 == 16): # if left knee marker of subject "left" missing, do "pattern fill" with LeftThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] + np.cumsum(diff_donor[0,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] + np.cumsum(diff_donor[1,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] + np.cumsum(diff_donor[2,:idx_first_val][::-1])[::-1]
                            
                        elif (subj==1) and (m//3 == 19): # if right knee marker of subject "right" missing, do "pattern fill" with RightThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] + np.cumsum(diff_donor[0,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] + np.cumsum(diff_donor[1,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] + np.cumsum(diff_donor[2,:idx_first_val][::-1])[::-1]
                        
                        else:
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] 
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] 
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] 
                        
                    # CHECK IF LAST FRAME IS MISSING
                    if np.isnan(xyz_vec[m,-1]):
                        idx_last_val = np.max(np.where(~np.isnan(xyz_vec[m,:])))
                        if m//3 in [0,3]: # if back head markers, do "pattern fill" with front head markers
                            if m//3==0: diff_donor = np.diff(xyz_vec[3:6,:], axis=1) # Left back
                            if m//3==3: diff_donor = np.diff(xyz_vec[6:9,:], axis=1) # Right back
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] + np.cumsum(diff_donor[0,idx_last_val:])
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] + np.cumsum(diff_donor[1,idx_last_val:])
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val] + np.cumsum(diff_donor[2,idx_last_val:])
                        
                        elif (subj==0) and (m//3 == 16): # if left knee marker of subject "left" missing, do "pattern fill" with LeftThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] + np.cumsum(diff_donor[0,idx_last_val:])
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] + np.cumsum(diff_donor[1,idx_last_val:])
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val] + np.cumsum(diff_donor[2,idx_last_val:])
                            
                        elif (subj==1) and (m//3 == 19): # if right knee marker of subject "right" missing, do "pattern fill" with RightThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] + np.cumsum(diff_donor[0,idx_last_val:])
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] + np.cumsum(diff_donor[1,idx_last_val:])
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val] + np.cumsum(diff_donor[2,idx_last_val:])
                            
                        else:
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] 
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] 
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val]
            
            # Downsample data if fps!=fps_orig (250) and there is no nan in the trial (just a sanity check, nans should all have been removed in the former step)
            if (fps != 250) and (False in np.isnan(xyz_vec[:])) : 
                samps = int(DUR*fps)
                xyz_vec_ds=np.zeros((xyz_vec.shape[0],samps))  
                for i in range(xyz_vec_ds.shape[0]): 
                    xyz_vec_ds[i,:]=np.interp(np.linspace(0.0, 1.0, samps, endpoint=False), np.linspace(0.0, 1.0,  xyz_vec.shape[1], endpoint=False), xyz_vec[i,:])
                xyz_vec = xyz_vec_ds
                
            # Reshape: from (Nmarkers*3, Time) to (Nmarkers, 3, Time)
            sz = xyz_vec.shape
            xyz_vec_resh = np.reshape(xyz_vec, (sz[0]//3,3,sz[1]))
            
            # Re-orient the data as both subjects were facing each other
            orient = np.sign(-1 * xyz_vec_resh[ 13 , 1 , 0])
            xyz_vec_resh[:,0,:] = xyz_vec_resh[:,0,:] * orient;   xyz_vec_resh[:,1,:] = xyz_vec_resh[:,1,:] * orient
            
            # Trim to extract PMs only when they listen to music
            if d==1: tBefore=10; tAfter=5
            else: tBefore=8; tAfter=7
            xyz_vec_resh = xyz_vec_resh[:,:,tBefore*fps:DUR*fps - tAfter*fps]
            
            songsTrial = song[:,(d-1),tr]   # the two songs of this dyad and trial
            tStop = int( min(musParts_tFrames[songsTrial , -1]) * fps ) # take the end frame of the shortest song between two subjects    
            xyz_vec_resh = xyz_vec_resh[:,:,:tStop]
            
            # Define the origin of the reference system (average position of point between feet for this trial (to avoid inter-trial offsets))
            Ox = np.mean((xyz_vec_resh[18,0,:]+xyz_vec_resh[21,0,:])/2);   Oy = np.mean((xyz_vec_resh[18,1,:]+xyz_vec_resh[21,1,:])/2);   Oz = np.mean((xyz_vec_resh[18,2,:]+xyz_vec_resh[21,2,:])/2)
            xyz_vec_resh[:,0,:] -= Ox;   xyz_vec_resh[:,1,:] -= Oy;   xyz_vec_resh[:,2,:] -= Oz
            
            # Reshape: come back to initial (Nmarkers*3, Time)
            sz_resh = xyz_vec_resh.shape
            xyz_vec = np.reshape(xyz_vec_resh, (sz_resh[0]*sz_resh[1] , sz_resh[2]))
            
            # Remove IdThighAdd markers (used only to distinguish between the two skeletons in the room)
            xyz_vec = xyz_vec[:66,:]
            
            ############## NORMALIZATION FOR MULTI-SUBJECT PM ANALYSIS ##############
            # De-mean by the average posture of the trial
            pmean_tr[(d-1)*2+subj,:,tr] = np.mean(xyz_vec,1)
            xyz_vec -= pmean_tr[(d-1)*2+subj,:,tr].reshape((-1,1))
            
            # Divide by the general std over all markers (this way, every subject contributes equally to the variance captured by PCA)
            # dmean_tr[(d-1)*2+subj,tr] = np.mean( np.linalg.norm(xyz_vec,axis=0) )
            dmean_tr[(d-1)*2+subj,tr] = np.std( xyz_vec[:] )
            xyz_vec /= dmean_tr[(d-1)*2+subj,tr].reshape((-1,1))
            
            # Store data
            data_subj[(d-1)*2+subj,:,iStart:iStart+tStop] = xyz_vec.copy()
            
            # Sanity check (only one trial missing for one dyad)
            if True in np.isnan(xyz_vec):
                print('NAN TRIAL: Participant ' + str((d-1)*2+subj) + ' - Trial ' + str(tr+1) + ' - Length ' + str(tStop))
            
            iStart+=tStop


#%% 
##############################################################
############     EXTRACT PRINCIPAL MOVEMENTS      ############
##############################################################

print('------------------------')
print('EXTRACTING PRINCIPAL MOVEMENTS...')
print('------------------------')

# Combine data into a matrix usable for PCA (the 66 channels by time; concatenated over trials and participants)
pos_mat=data_subj[0,:,:]
for i in range(1,NB_DYADS*2):
    pos_mat= np.hstack((pos_mat,data_subj[i,:,:][:,np.where(~np.isnan(data_subj[i,0,:]))[0]]))
pos_mat = pos_mat.T
del data_subj

# Apply PCA using Singular Value Decomposition
U, S, V = np.linalg.svd(pos_mat, full_matrices=False)
e=S**2
common_nrj = np.cumsum(e) / np.sum(e);     nbEigen = [i for (i, val) in enumerate(common_nrj) if val>0.95][0];
common_PC_scores = (U*S)
common_eigen_vects = V
# del pos_mat

#%% 
##############################################################
############        VALIDATION OF THE PMs         ############
##############################################################

print('------------------------')
print('VALIDATING THE PMs...')
print('------------------------')

################ SPECTRAL ANALYSIS OF THE PMs ################
# Compute Power Spectral Density (PSD) using Welch method
NB_PM = 14

Nwin=fps
freqs, psd = signal.welch(common_PC_scores[:,:NB_PM],fps,axis=0,nperseg=Nwin,noverlap=3*Nwin//4)

# PLOT PSD
fig = plt.figure(figsize=(8,20));
st = fig.suptitle('Frequency content of the first ' + str(NB_PM) + ' common PMs')
for pm in range(NB_PM):
    ax = fig.add_subplot(NB_PM//2,2,pm+1)
    plt.plot(freqs,psd[:,pm],c='tab:blue',label='mean');  
    plt.xticks(np.arange(min(freqs), max(freqs)+1, 5));  ax.set_xlim(0,min(20,fps//2)); 
    plt.title('PM ' + str(pm+1))
    if (i+1)%2==1: plt.ylabel('PSD (amp**2/Hz)')
    if i==6 or i == 7 : plt.xlabel('Frequency (Hz)')
fig.tight_layout()
fig.savefig(output_dir + '/spectral_analysis.png', dpi=600, bbox_inches='tight'); plt.close()  

# FILTER the PC_scores
print('LOW-PASS FILTERING...')
fc = 6  # Cut-off frequency of the filter
w = fc / (fps / 2) # Normalize the frequency
b, a = signal.butter(4, w, 'low')  # 4th-order Butterwoth filter

iStart=0
for i in range(NB_DYADS*2):
    if i in [62,63]: nbTr = NB_TRIALS - 1
    else: nbTr = NB_TRIALS
    for tr in range(nbTr):
        tStop = int( min(musParts_tFrames[songsTrial , -1]) * fps ) # take the end frame of the shortest song between two subjects    
        common_PC_scores[iStart:iStart+tStop,:] = signal.filtfilt(b, a, common_PC_scores[iStart:iStart+tStop,:], axis=0)
        
        iStart += tStop

#%% 
##############################################################
############           VIZ PMs AS 2D PLOTS         ############
##############################################################
            
pos_viz=0
if pos_viz==1:
    print('------------------------')
    print('PLOTTING THE PMs AS 2D PLOTS...')
    print('------------------------')
    
    ##### VISUALIZE THE N FIRST PMs (2-post graph with the min and max PM postures across signers) #####
    NB_PM = 15
    for pm in range(NB_PM):
        EXAG = 1
        if pm in [2]: EXAG = 1.5
        if pm in [3,6,10]: EXAG = 2
        if pm in [7]: EXAG = 3
        if pm in [8]: EXAG = 4
        common_eigenmov = np.outer(EXAG*common_PC_scores[:,pm] , common_eigen_vects[pm,:]).T
        iStart=0; 
        for i in range(NB_DYADS*2) :
            if i in [62,63]: nbTr = NB_TRIALS - 1
            else: nbTr = NB_TRIALS
            for tr in range(nbTr):
                songsTrial = song[:,i//2,tr]
                tStop = int( min(musParts_tFrames[songsTrial , -1]) * fps ) # take the end frame of the shortest song between two subjects    
                common_eigenmov[:,iStart:iStart+tStop] *= dmean_tr[i,tr]
                
                iStart += tStop
        
        common_eigenmov += np.mean(np.nanmean(pmean_tr,2),0).reshape((-1,1))
        common_eigenmov = np.reshape( common_eigenmov ,(NB_MARKERS , 3 , common_PC_scores.shape[0]) )
        
        # plot postures at min and max of the PM weightings
        i_min = argmin(common_PC_scores[:,pm]); i_max = argmax(common_PC_scores[:,pm])
    
        times=[i_min,i_max]
        print('PM ' + str(pm+1) + ' - min: ' + str(i_min) + ' - max: ' + str(i_max))
        
        plot_2frames(common_eigenmov,times,"XZ",liaisons=liaisons,center_sens=13,save_dir=output_dir + '/PM' + str(pm+1) + '_XZ.png'); plt.close()
        plot_2frames(common_eigenmov,times,"YZ",liaisons=liaisons,center_sens=13,save_dir=output_dir + '/PM' + str(pm+1) + '_YZ.png'); plt.close()

#%% 
##############################################################
############         VIZ PMs AS PL VIDEOS         ############
##############################################################

video=0
if video==1:
    print('------------------------')
    print('VISUALIZING THE PMs AS PL VIDEOS...')
    print('------------------------')
    
    SUBJ = 10
    DUR_VID = 5
    NB_PM = 15
    for pm in range(NB_PM) :  
        EXAG = 1
        if pm in [2]: EXAG = 1.5
        if pm in [3,6,10,13]: EXAG = 2
        if pm in [7]: EXAG = 3
        if pm in [8]: EXAG = 4
        common_eigenmov = np.outer(EXAG*common_PC_scores[:,pm] , common_eigen_vects[pm,:]).T
        iStart=0
        for i in range(NB_DYADS*2):
            if i in [62,63]: nbTr = NB_TRIALS - 1
            else: nbTr = NB_TRIALS
            for tr in range(nbTr):
                songsTrial = song[:,i//2,tr]
                tStop = int( min(musParts_tFrames[songsTrial , -1]) * fps ) # take the end frame of the shortest song between two subjects    
                common_eigenmov[:,iStart:iStart+tStop] *= dmean_tr[i,tr]
                
                iStart += tStop
        
        common_eigenmov += np.mean(np.nanmean(pmean_tr,2),0).reshape((-1,1))
        common_eigenmov = common_eigenmov[: , (SUBJ+1)*LEN_CONCAT_TIME-DUR_VID*fps : (SUBJ+1)*LEN_CONCAT_TIME]
        
        # Downsample for video
        fps_VID = 25
        if fps != fps_VID : 
            samps = int(DUR_VID*fps_VID)
            common_eigenmov_ds=np.zeros((common_eigenmov.shape[0],samps))  
            for i in range(common_eigenmov_ds.shape[0]): 
                common_eigenmov_ds[i,:]=np.interp(np.linspace(0.0, 1.0, samps, endpoint=False), np.linspace(0.0, 1.0,  common_eigenmov.shape[1], endpoint=False), common_eigenmov[i,:])
            common_eigenmov = common_eigenmov_ds
        
        common_eigenmov = np.reshape( common_eigenmov , (NB_MARKERS, 3, DUR_VID*fps_VID ) )
        
        # Synthesize PL video
        minDim = np.array([(common_eigenmov[:,0,:]).min() , (common_eigenmov[:,1,:]).min() , 0])
        maxDim = np.array([(common_eigenmov[:,0,:]).max() , (common_eigenmov[:,1,:]).max() , (common_eigenmov[:,2,:]).max()])
        if not (os.path.exists(output_dir + '/PM' + str(pm+1) + '_subj-' + str(SUBJ+1) + '.mp4')) :
            video_PL(common_eigenmov, output_dir + '/PM' + str(pm+1) + '_subj-' + str(SUBJ+1) + '.mp4', minDim=minDim, maxDim=maxDim, fps=fps_VID, title='PM ' + str(pm+1) + '\n(subject ' + str(SUBJ+1) + ')')
            

#%% 
##############################################################
#########       COMPUTE QUANTITY OF MOVEMENT         #########
##############################################################

print('------------------------')
print('COMPUTING QUANTITY OF MOVEMENT...')
print('------------------------')

# Choose how many PMs to run the analysis on
NB_PM = 15

# Find number of frames for shortest song (during song + during silence before & after (3 bars of silence))
tStop_min = int( min(musParts_tFrames[: , -1]) * fps )
Nsamps_before_min = int( 3 * tStop_min / 20 )
NB_T = tStop_min + Nsamps_before_min*2

# Convert time frames in metrical frames
t_norm = np.arange(0,NB_T)/fps;        dt_norm = t_norm[1]
bpm_max = max(song_bpms)
t_norm_beat = t_norm/(60/bpm_max) + 1

# Cross-wavelet transform (XWT) parameters
NB_OCT = 6;   NB_SCALES_PER_OCT=16
dt=1/fps; dj=1/NB_SCALES_PER_OCT; J=NB_OCT/dj
NB_FREQ = int(J)+1

# Matrices of XWT for ANOVA averaged over movement frequencies
QoM_formatJASP = np.zeros( (NB_T, NB_PM, NB_DYADS , 4) ) # for each marker: Nparticipants x 4 conditions

for d in range(1, NB_DYADS+1):
    if d < 10 : numDyad = "0" + str(d)
    else : numDyad = str(d)
    print("Dyad " + numDyad)
    
    # Motion data (from Vicon)
    folder = os.path.normpath(( os.getcwd() + "/vicondata/npy/" + folders[d-1] ))
    files=os.listdir(folder + '/subj_LEFT');      files=[x for i,x in enumerate(files) if (x.endswith(".npy"))]
    files=sorted(files)
    
    # Condition & song info (from Presentation)
    logfile = pd.read_excel( os.path.normpath(( os.getcwd() + "/Presentation_log/" + folders[d-1]  + ".xlsx")) , engine='openpyxl')
    logfile_sounds = logfile.loc[logfile.iloc[:,1]== "Sound"]
    log_conditions = logfile_sounds.iloc[np.where(logfile_sounds.iloc[:, 2].str.contains('Start_Trial') == True)[0], 2].to_numpy().astype(str)
    log_songs = logfile_sounds.iloc[np.where(logfile_sounds.iloc[:, 2].str.contains('Start_Trial') == True)[0] + 1, 2].to_numpy().astype(str)
    for tr in range(32):
        cond_mus[(d-1) , tr] = int('Same Music' in log_conditions[tr][12:])
        cond_vis[(d-1) , tr] = int('Yes Visual' in log_conditions[tr][12:])
        song[0,(d-1) , tr] = int(log_songs[tr][16]) - 1
        song[1,(d-1) , tr] = int(log_songs[tr][24]) - 1 
    
    QoM_dyad = np.zeros(( NB_T, NB_PM , NB_TRIALS ))
    for tr in range(NB_TRIALS):
        if tr+1 < 10 : numTrial = "0" + str(tr+1)
        else : numTrial = str(tr+1)
        songsTrial = song[:,(d-1),tr]
        freqsBeat = 1/(60/song_bpms[songsTrial])
        
        # SUBJECT 1 AND 2
        for subj in range(2):   
                
            # Load motion data
            if subj == 0 : xyz_vec=np.load( os.path.normpath( folder + '/subj_LEFT/' + files[tr] ))
            if subj == 1 : xyz_vec=np.load( os.path.normpath( folder + '/subj_RIGHT/' + files[tr] )) 
            
            # Small last pre-processing step
            # Fill the data if NaN at start and/or at end (as Nexus only allows filling gaps if there is data before and after)
            # (for future: put it as a function of PL MOCAP)
            if False in np.isnan(xyz_vec[:]) : 
                for m in range(0,NB_MARKERS*3,3):
                    # CHECK IF FIRST FRAME IS MISSING
                    if np.isnan(xyz_vec[m,0]):
                        idx_first_val = np.min(np.where(~np.isnan(xyz_vec[m,:])))
                        
                        if m//3 in [0,3]: # if back head markers, do "pattern fill" with front head markers
                            if m//3==0: diff_donor = np.diff(xyz_vec[3:6,:], axis=1) # Left back
                            if m//3==3: diff_donor = np.diff(xyz_vec[6:9,:], axis=1) # Right back
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] + np.cumsum(diff_donor[0,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] + np.cumsum(diff_donor[1,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] + np.cumsum(diff_donor[2,:idx_first_val][::-1])[::-1]
                        
                        elif (subj==0) and (m//3 == 16): # if left knee marker of subject "left" missing, do "pattern fill" with LeftThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] + np.cumsum(diff_donor[0,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] + np.cumsum(diff_donor[1,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] + np.cumsum(diff_donor[2,:idx_first_val][::-1])[::-1]
                            
                        elif (subj==1) and (m//3 == 19): # if right knee marker of subject "right" missing, do "pattern fill" with RightThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] + np.cumsum(diff_donor[0,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] + np.cumsum(diff_donor[1,:idx_first_val][::-1])[::-1]
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] + np.cumsum(diff_donor[2,:idx_first_val][::-1])[::-1]
                        
                        else:
                            xyz_vec[m,:idx_first_val] = xyz_vec[m,idx_first_val] 
                            xyz_vec[m+1,:idx_first_val] = xyz_vec[m+1,idx_first_val] 
                            xyz_vec[m+2,:idx_first_val] = xyz_vec[m+2,idx_first_val] 
                        
                    # CHECK IF LAST FRAME IS MISSING
                    if np.isnan(xyz_vec[m,-1]):
                        idx_last_val = np.max(np.where(~np.isnan(xyz_vec[m,:])))
                        if m//3 in [0,3]: # if back head markers, do "pattern fill" with front head markers
                            if m//3==0: diff_donor = np.diff(xyz_vec[3:6,:], axis=1) # Left back
                            if m//3==3: diff_donor = np.diff(xyz_vec[6:9,:], axis=1) # Right back
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] + np.cumsum(diff_donor[0,idx_last_val:])
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] + np.cumsum(diff_donor[1,idx_last_val:])
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val] + np.cumsum(diff_donor[2,idx_last_val:])
                        
                        elif (subj==0) and (m//3 == 16): # if left knee marker of subject "left" missing, do "pattern fill" with LeftThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] + np.cumsum(diff_donor[0,idx_last_val:])
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] + np.cumsum(diff_donor[1,idx_last_val:])
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val] + np.cumsum(diff_donor[2,idx_last_val:])
                            
                        elif (subj==1) and (m//3 == 19): # if right knee marker of subject "right" missing, do "pattern fill" with RightThighAdd marker
                            diff_donor = np.diff(xyz_vec[66:69,:], axis=1) # Right Thigh Add
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] + np.cumsum(diff_donor[0,idx_last_val:])
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] + np.cumsum(diff_donor[1,idx_last_val:])
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val] + np.cumsum(diff_donor[2,idx_last_val:])
                            
                        else:
                            xyz_vec[m,idx_last_val+1:] = xyz_vec[m,idx_last_val] 
                            xyz_vec[m+1,idx_last_val+1:] = xyz_vec[m+1,idx_last_val] 
                            xyz_vec[m+2,idx_last_val+1:] = xyz_vec[m+2,idx_last_val]
            
            # If there is no nan in the trial (just a sanity check, nans should all have been removed in the former step)
            if False in np.isnan(xyz_vec[:]): 
                # Downsample data if fps!=fps_orig (250)
                if fps!= 250:
                    samps = int(DUR*fps)
                    xyz_vec_ds=np.zeros((xyz_vec.shape[0],samps))  
                    for i in range(xyz_vec_ds.shape[0]): 
                        xyz_vec_ds[i,:]=np.interp(np.linspace(0.0, 1.0, samps, endpoint=False), np.linspace(0.0, 1.0,  xyz_vec.shape[1], endpoint=False), xyz_vec[i,:])
                    xyz_vec = xyz_vec_ds
                    
                # # Low-pass filter
                # fc = 6  # Cut-off frequency of the filter
                # w = fc / (fps / 2.) # Normalize the frequency
                # b, a = signal.butter(4, w, 'low')
                # xyz_vec = signal.filtfilt(b, a, xyz_vec, axis=1)
                
            # Reshape: from (Nmarkers*3, Time) to (Nmarkers, 3, Time)
            sz = xyz_vec.shape
            xyz_vec_resh = np.reshape(xyz_vec, (sz[0]//3,3,sz[1]))
            
            # Re-orient the data as both subjects were facing each other
            orient = np.sign(-1 * xyz_vec_resh[ 13 , 1 , 0])
            xyz_vec_resh[:,0,:] = xyz_vec_resh[:,0,:] * orient;   xyz_vec_resh[:,1,:] = xyz_vec_resh[:,1,:] * orient
            
            # Trim to measure XWT of the PMs only 3 bars silence before - 20 bars music - 3 bars silence after
            # Trim 7s before + 45s + 5s after (to have enough time covering 3 bars before and after for every song, this trimming is done further below)
            if d==1: tBefore=3; tAfter=0
            else: tBefore=1; tAfter=2
            xyz_vec_resh = xyz_vec_resh[:,:,tBefore*fps:DUR*fps - tAfter*fps]
            
            
            # Define the origin of the reference system (average position of point between feet for this trial (to avoid inter-trial offsets))
            Ox = np.mean((xyz_vec_resh[18,0,:]+xyz_vec_resh[21,0,:])/2);   Oy = np.mean((xyz_vec_resh[18,1,:]+xyz_vec_resh[21,1,:])/2);   Oz = np.mean((xyz_vec_resh[18,2,:]+xyz_vec_resh[21,2,:])/2)
            xyz_vec_resh[:,0,:] -= Ox;   xyz_vec_resh[:,1,:] -= Oy;   xyz_vec_resh[:,2,:] -= Oz
            
            # Reshape: come back to initial (Nmarkers*3, Time)
            sz_resh = xyz_vec_resh.shape
            xyz_vec = np.reshape(xyz_vec_resh, (sz_resh[0]*sz_resh[1] , sz_resh[2]))
            
            # Remove IdThighAdd markers (used only to distinguish between the two skeletons in the room)
            xyz_vec = xyz_vec[:66,:]
            
            ############## NORMALIZATION FOR MULTI-SUBJECT PM ANALYSIS (Federolf et al. 2013) ##############
            # De-mean by the average posture of the trial
            xyz_vec -= pmean_tr[(d-1)*2+subj,:,tr].reshape((-1,1))
            
            # Normalize by the std
            # xyz_vec /= dmean_tr[(d-1)*2+subj,tr].reshape((-1,1))
            
            # Compute velocity of the PMs
            # Compute number of frames for each song that represent 3 bars (to trim 3 bars silence - 20 bars music - 3 bars silence)
            lensong_diff = np.around(musParts_tFrames[songsTrial , -1] * fps).astype(int)
            lensong = min(lensong_diff)
            Nsampes_before_diff = np.around(3 * lensong_diff / 20 ).astype(int)     # 3 bars
            Nsampes_before = int(3 * lensong / 20 )    # 3 bars
            if subj==0: 
                pos_subj1 = np.dot(xyz_vec.T , common_eigen_vects.T).T          # Get PM score by projecting the trajectory onto the eigen vector
                
                
                
                traj_subj1 = np.gradient( pos_subj1 , 1/fps , axis=1 )          # Compute velocity
                
                # Low-pass filter the PM
                fc = 12  # Cut-off frequency of the filter
                w = fc / (fps / 2.) # Normalize the frequency
                b, a = signal.butter(2, w, 'low')
                traj_subj1 = signal.filtfilt(b, a, traj_subj1, axis=1) 
                

            if subj==1:   # same but for subject 2
                pos_subj2 = np.dot(xyz_vec.T , common_eigen_vects.T).T
                
                traj_subj2 = np.gradient( pos_subj2 , 1/fps , axis=1 ) 
                
                # Low-pass filter the PM
                fc = 12  # Cut-off frequency of the filter
                w = fc / (fps / 2.) # Normalize the frequency
                b, a = signal.butter(2, w, 'low')
                traj_subj2 = signal.filtfilt(b, a, traj_subj2, axis=1) 


        ####### summed QoM #######
        if True not in np.isnan(traj_subj1[:]):
            Nsamps_before_min = int( 3 * tStop_min / 20 )
            Nsamps = tStop_min + Nsamps_before_min * 2      # total number of frames of shortest trial
            print('XWT... Trial ' + str(numTrial))
            sumQoM = abs(traj_subj1) + abs(traj_subj2)
            for pm in range(NB_PM):
                QoM_dyad[:,pm,tr] = np.interp(np.linspace(0.0, 1.0, Nsamps), np.linspace(0.0, 1.0,  len(sumQoM[pm,:])), sumQoM[pm,:])
                
                # SMOOTH with rolling average window of 3 bars
                kernel_size = Nsamps_before_min # 3 bars
                kernel = np.ones(kernel_size) / kernel_size
                QoM_dyad[:,pm,tr] = np.convolve(QoM_dyad[:,pm,tr], kernel, mode='same')
  
        else:       # for only one trial that is missing: put nans
            QoM_dyad[:,:,tr] = np.nan
      
        
    ########### MEANS OF CONDITIONS, OVER PARTICIPANTS, FOR ANOVA ############ 
    # Create mask of conditions
    YesVis_mask = (cond_vis[(d-1),:]==1); NoVis_mask = (cond_vis[(d-1),:]==0)
    SameMus_mask = (cond_mus[(d-1),:]==1); DiffMus_mask = (cond_mus[(d-1),:]==0)
    
    # IMS data of this dyad for each condition
    QoM_dyad_YesVisSameMus = QoM_dyad[:,:,YesVis_mask & SameMus_mask]
    QoM_dyad_YesVisDiffMus = QoM_dyad[:,:,YesVis_mask & DiffMus_mask]
    QoM_dyad_NoVisSameMus = QoM_dyad[:,:,NoVis_mask & SameMus_mask]
    QoM_dyad_NoVisDiffMus = QoM_dyad[:,:,NoVis_mask & DiffMus_mask]
    
    # Average over trials
    QoM_dyad_YesVisSameMus_mean = np.nanmean(QoM_dyad_YesVisSameMus,axis=2)
    QoM_dyad_YesVisDiffMus_mean = np.nanmean(QoM_dyad_YesVisDiffMus,axis=2)
    QoM_dyad_NoVisSameMus_mean = np.nanmean(QoM_dyad_NoVisSameMus,axis=2)
    QoM_dyad_NoVisDiffMus_mean = np.nanmean(QoM_dyad_NoVisDiffMus,axis=2)

    # Store for ANOVA
    QoM_formatJASP[:,:,(d-1),0] = QoM_dyad_YesVisSameMus_mean
    QoM_formatJASP[:,:,(d-1),1] = QoM_dyad_YesVisDiffMus_mean
    QoM_formatJASP[:,:,(d-1),2] = QoM_dyad_NoVisSameMus_mean
    QoM_formatJASP[:,:,(d-1),3] = QoM_dyad_NoVisDiffMus_mean
    
    IMS_trials[(d-1),:,:,:] = IMS_dyad.mean(axis=1)
    

#%% 
##############################################################
######  ANOVA ANALYSIS WITH CLUSTER-BASED PERMUTATION   ######
##############################################################

print('------------------------')
print('RUNNING THE ANOVA ANALYSIS...')
print('------------------------')

# Choose number of PMs to analyse
NB_PM=15
NB_T = tStop_min + Nsamps_before_min*2 

# Create adjacency matrix for the temporal clusters (neighbours are adjacent bins in time)
nei_mask_time = sparse.diags([1., 1.], offsets=(-1, 1), shape=(NB_T - 2*Nsamps_before_min, NB_T - 2*Nsamps_before_min))
adjacency = mne_fefe.stats.combine_adjacency( nei_mask_time )

# Average the XWT power over frequencies before running ANOVA on the time-series
X_lin = QoM_formatJASP[:,:NB_PM,:,:]
 

X_lin = np.transpose(X_lin, [2, 0, 1, 3])   # reshape to good format for MNE library ()
X_lin_norm = X_lin.copy()

# Normalize to reduce inter-dyad differences
for pm in range(14):
    for d in range(35):
        X_lin_norm[d,:,pm,:] /= X_lin_norm[d,:,pm,:].std()
X=X_lin_norm

X = [np.squeeze(x) for x in np.split(X, 4, axis=-1)]
factor_levels = [2, 2]

# Difference of means between the main factors (so "YesVis vs. NoVis"; "SameMus vs DiffMus"; "(SameMus vs DiffMus) if YesVis vs. if NoVis")
diffMeans = np.zeros( (3 , NB_T, NB_PM))
diffMeans[0,:,:] = ( (X_lin[:,:,:,0].mean(axis=0) + X_lin[:,:,:,1].mean(axis=0))/2 ) - ( (X_lin[:,:,:,2].mean(axis=0) + X_lin[:,:,:,3].mean(axis=0))/2 )
diffMeans[1,:,:] = ( (X_lin[:,:,:,0].mean(axis=0) + X_lin[:,:,:,2].mean(axis=0))/2 ) - ( (X_lin[:,:,:,1].mean(axis=0) + X_lin[:,:,:,3].mean(axis=0))/2 )
diffMeans[2,:,:] = ( X_lin[:,:,:,0].mean(axis=0) - X_lin[:,:,:,1].mean(axis=0) ) - ( X_lin[:,:,:,2].mean(axis=0) - X_lin[:,:,:,3].mean(axis=0) )

# Store the cluster locations (tStart and tStop) + the pValue of the cluster for each of the 3 effects
cluster_start_vis = []; cluster_start_mus = []; cluster_start_int = []; clusterSIG_start_vis = []; clusterSIG_start_mus = []; clusterSIG_start_int = [];
cluster_stop_vis = []; cluster_stop_mus = []; cluster_stop_int = []; clusterSIG_stop_vis = []; clusterSIG_stop_mus = []; clusterSIG_stop_int = [];
p_vis = []; p_mus = []; p_int = [];

# Run cluster permutation independently on each PM (then correct for this)
for pm in range(NB_PM):     
    print('\n-------\nPRINCIPAL MOVEMENT PM ' + str(pm+1) + '\n-------\n')
    X_pm = [x[:,Nsamps_before_min:-Nsamps_before_min,pm] for x in X]        # analyse only when they listen to music

    for effect in range(3): 
        if effect==0:
            print('-------\nMain effect Visual \n-------')
            effects = 'A'
            def stat_fun(*args):
                # get f-values only + weight the Fvalue by the sign of the difference (to avoid clusters of different sign)
                diffMeansVis = (args[0].mean(axis=0) + args[1].mean(axis=0))/2 - (args[2].mean(axis=0) + args[3].mean(axis=0))/2 
                return mne_fefe.stats.f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,
                                 effects=effects, return_pvals=False)[0] * np.sign(diffMeansVis)
        
        if effect==1:
            print('-------\nMain effect Music \n-------')
            effects = 'B'
            def stat_fun(*args):
                # get f-values only + weight the Fvalue by the sign of the difference (to avoid clusters of different sign)
                diffMeansMus = (args[0].mean(axis=0) + args[2].mean(axis=0))/2 - (args[1].mean(axis=0) + args[2].mean(axis=0))/2 
                return mne_fefe.stats.f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,
                                 effects=effects, return_pvals=False)[0] * np.sign(diffMeansMus)
            
        if effect==2:
            print('-------\nInteraction Visual x Music \n-------')
            effects = 'A:B'
            def stat_fun(*args):
                # get f-values only + weight the Fvalue by the sign of the difference (to avoid clusters of different sign)
                diffMeans_IMS_mus_YesVis = args[0].mean(axis=0) - args[1].mean(axis=0)
                diffMeans_IMS_mus_NoVis = args[2].mean(axis=0) - args[3].mean(axis=0)
                diff_interact = diffMeans_IMS_mus_YesVis - diffMeans_IMS_mus_NoVis
                return mne_fefe.stats.f_mway_rm(np.swapaxes(args, 1, 0), factor_levels=factor_levels,
                                  effects=effects, return_pvals=False)[0] * np.sign(diff_interact)
        
        
        ######## CLUSTERING ########
        pthresh = 0.05 / NB_PM          # Threhsold pValue to create initial clusters
        # pthresh = 0.05
        n_subjects=NB_DYADS
        f_thresh = mne_fefe.stats.f_threshold_mway_rm(n_subjects, factor_levels, effects, pthresh)      # Threhsold Fvalue
        
        # Define number of permutations
        n_permutations = 10001 
        
        print('Clustering.')
        F_obs, clusters, cluster_p_values, H0 = clu = \
            mne_fefe.stats.cluster_level.spatio_temporal_cluster_test(X_pm, adjacency=adjacency, n_jobs=None,
                                         threshold=f_thresh, stat_fun=stat_fun,
                                         n_permutations=n_permutations,t_power=1,
                                         buffer_size=None,out_type='indices',stat_cluster='sum')
        
        # Retain clusters of at least 2 temporal bins
        clusters_filt = [x for i, x in enumerate(clusters) if len(x[x == True]) > 1]
        idx = [i for i, x in enumerate(clusters) if len(x[x == True]) > 1]
        cluster_p_values_filt = cluster_p_values[idx]

        # Define pValue threshold for significant clusters by contrast with random distribution of n_permutations clusters
        final_thresh = 0.05 / NB_PM
        # final_thresh = 0.05
        idx_cluster_sig = np.where(cluster_p_values_filt<final_thresh)[0]
        
        # Store tStart, tStop and pValue of each significant cluster
        if len(idx_cluster_sig)>0:
            idx_cluster_sig = np.where(cluster_p_values_filt<final_thresh)[0]
            first_vals = np.zeros(len(idx_cluster_sig),dtype=int); last_vals = np.zeros(len(idx_cluster_sig),dtype=int); p_corr = np.zeros(len(idx_cluster_sig));
            for c in range(len(idx_cluster_sig)):
                first_vals[c] = clusters_filt[idx_cluster_sig[c]][0][0] + Nsamps_before_min
                last_vals[c] =  clusters_filt[idx_cluster_sig[c]][0][-1] + Nsamps_before_min
                p_corr[c] = cluster_p_values_filt[idx_cluster_sig[c]] 
            
            if effect==0: clusterSIG_start_vis.append(first_vals); clusterSIG_stop_vis.append(last_vals); p_vis.append(p_corr)
            if effect==1: clusterSIG_start_mus.append(first_vals); clusterSIG_stop_mus.append(last_vals); p_mus.append(p_corr)
            if effect==2: clusterSIG_start_int.append(first_vals); clusterSIG_stop_int.append(last_vals); p_int.append(p_corr)
            
        else: 
            if effect==0: cluster_start_vis.append(np.empty(0)); cluster_stop_vis.append(np.empty(0)); p_vis.append(np.empty(0))
            if effect==1: cluster_start_mus.append(np.empty(0)); cluster_stop_mus.append(np.empty(0)); p_mus.append(np.empty(0))
            if effect==2: cluster_start_int.append(np.empty(0)); cluster_stop_int.append(np.empty(0)); p_int.append(np.empty(0))
            
            if effect==0: clusterSIG_start_vis.append(np.empty(0)); clusterSIG_stop_vis.append(np.empty(0))
            if effect==1: clusterSIG_start_mus.append(np.empty(0)); clusterSIG_stop_mus.append(np.empty(0))
            if effect==2: clusterSIG_start_int.append(np.empty(0)); clusterSIG_stop_int.append(np.empty(0))


#%% 
##############################################################
######  PLOT THE RESULTS OF THE ANOVA WITH CLUSTERS     ######
######                                                  ######
######  1. PREPARE PM VIZ BY DETECTING MIN AND MAX      ######
##############################################################



# Adjust the exaggeration of the rendering of the PM min/max manually, for PMs that are hard to visualize/interpret
EXAG = np.zeros(NB_PM)
PM_min = np.zeros( (NB_PM , NB_MARKERS , 3) ); PM_max = np.zeros( (NB_PM , NB_MARKERS , 3) )
for pm in range(NB_PM):
    EXAG[pm] = 1
    # if pm in [13]: EXAG[pm] = 0.5
    if pm in [2]: EXAG[pm] = 1.5
    if pm in [3,6,10]: EXAG[pm] = 2
    if pm in [7]: EXAG[pm] = 3
    if pm in [8]: EXAG[pm] = 4
    common_eigenmov = np.outer(EXAG[pm]*common_PC_scores[:,pm] , common_eigen_vects[pm,:]).T
    iStart=0; 
    for i in range(NB_DYADS*2) :
        if i in [62,63]: nbTr = NB_TRIALS - 1
        else: nbTr = NB_TRIALS
        for tr in range(nbTr):
            songsTrial = song[:,i//2,tr]
            tStop = int( min(musParts_tFrames[songsTrial , -1]) * fps ) # take the end frame of the shortest song between two subjects    
            common_eigenmov[:,iStart:iStart+tStop] *= dmean_tr[i,tr]
            
            iStart += tStop
    
    common_eigenmov += np.mean(np.nanmean(pmean_tr,2),0).reshape((-1,1))
    common_eigenmov = np.reshape( common_eigenmov ,(NB_MARKERS , 3 , common_PC_scores.shape[0]) )
    
    # Find min and max of the PM weightings
    i_min = argmin(common_PC_scores[:,pm]) ; i_max = argmax(common_PC_scores[:,pm]) 
    
    PM_min[pm,:,:] = common_eigenmov[:,:,i_min];  PM_max[pm,:,:] = common_eigenmov[:,:,i_max]
    
#%% THE MAGIC PLOT THAT SUMMARIZES ALL OUR RESULTS

matplotlib.use('Agg')

NB_PM = 14

# Keep only one plan (2D)
minDim = np.array([ min((PM_min[:,:,0].min(),PM_max[:,:,0].min())) , min((PM_min[:,:,1].min(),PM_max[:,:,1].min())) , 0])
maxDim = np.array([ max((PM_min[:,:,0].max(),PM_max[:,:,0].max())) , max((PM_min[:,:,1].max(),PM_max[:,:,1].max())) , max((PM_min[:,:,2].max(),PM_max[:,:,2].max()))])

bigtitle=['VISUAL','MUSIC']
for effect in range(2):   
    fig = plt.figure(figsize=(25,20))
    gs = fig.add_gridspec(NB_PM//2 , 5, width_ratios=[3,10,0.5,3,10])
    
    for pm in range(NB_PM):
        m = pm%7; col=(pm//7)*3
        
        if effect==0: clusterSIG_start = clusterSIG_start_vis[pm]; clusterSIG_stop = clusterSIG_stop_vis[pm]; p_corr = p_vis[pm]; colors=['cornflowerblue','tomato']; label='visual'
        if effect==1: clusterSIG_start = clusterSIG_start_mus[pm]; clusterSIG_stop = clusterSIG_stop_mus[pm]; p_corr = p_mus[pm]; colors=['mediumpurple','forestgreen']; label='music'
        
        ### FIRST BLOCK: PMs DESCRIPTION ###
        ax1 = fig.add_subplot(gs[m, col+0])
        
        # Adjust your scale along the 3 axis 
        Kx = (maxDim[0] - minDim[0])*0.2; Ky = (maxDim[1] - minDim[1])*0.2; Kz = (maxDim[2] - minDim[2])*0.1;
        
        ax1.set_aspect('equal')
        ax1.get_xaxis().set_visible(False)  
        ax1.get_yaxis().set_visible(False)
        if pm in [3,7,9]: dim1=1; dim2=2 #YZ plane
        else:   dim1=0; dim2=2 #XZ plane
        for i in range(NB_MARKERS) :
            # MIN
            ax1.scatter(PM_min[pm,i,dim1], PM_min[pm,i,dim2],  c='gray', marker='o', s=55,  alpha=0.6)
            
            # MAX
            ax1.scatter(PM_max[pm,i,dim1], PM_max[pm,i,dim2],  c='black', marker='o', s=55, alpha=0.6)
            
        ax1.set_xlim(minDim[dim1]-Kx,maxDim[dim1]+Kx); ax.set_ylim(minDim[dim2]-Kz,maxDim[dim2]+Kz); ax.invert_xaxis()
        
        for l in liaisons :
            c1 = l[0];   c2 = l[1]  # get the two joints
            # MIN
            ax1.plot([PM_min[pm,c1,dim1], PM_min[pm,c2,dim1]], [PM_min[pm,c1,dim2], PM_min[pm,c2,dim2]], 'k-', lw=1.5, c='gray')   
            
            # MAX
            ax1.plot([PM_max[pm,c1,dim1], PM_max[pm,c2,dim1]], [PM_max[pm,c1,dim2], PM_max[pm,c2,dim2]], 'k-', lw=1.5, c='black')   
        ax1.set_axis_off()
        ax1.set_title('PM ' + str(pm+1) + '\na=' + format(EXAG[pm],'g'),rotation='horizontal',x=-0.3,y=0.5,fontsize=20)
        
        if m == 0 and col==0: 
            ax2 = fig.add_subplot(gs[m, col+2]); ax2.set_axis_off()
            ax2.set_title(r'$\bf{' + bigtitle[effect] + '}$\n' + 'ANOVA CLUSTERS ON QoM OVER TIME',fontsize=20,pad=30)
        
        ax3 = fig.add_subplot(gs[m, col+1])
        if effect==0:
            mean_high = (X_lin_norm[:,:,pm,0].mean(axis=0) + X_lin_norm[:,:,pm,1].mean(axis=0))/2 ; label_high = 'YesVis'
            mean_low = (X_lin_norm[:,:,pm,2].mean(axis=0) + X_lin_norm[:,:,pm,3].mean(axis=0))/2 ; label_low = 'NoVis'
        
        if effect==1:
            mean_high = (X_lin_norm[:,:,pm,0].mean(axis=0) + X_lin_norm[:,:,pm,2].mean(axis=0))/2 ; label_high = 'SameMus'
            mean_low = (X_lin_norm[:,:,pm,1].mean(axis=0) + X_lin_norm[:,:,pm,3].mean(axis=0))/2 ; label_low = 'DiffMus'
            
        ax3.plot(t_norm_beat,mean_high,label=label_high,c='black')
        ax3.plot(t_norm_beat,mean_low,label=label_low,c='black', linestyle='dashed')
        ymin = 0; ymax = max(X_lin_norm[:,:,pm,:].mean(axis=0).flatten())*1.1
        # if m==0: ax3.set_ylabel('Cross power' , fontsize = 15)
        
        for c in range(len(clusterSIG_start)):
            assert( min(np.sign(diffMeans[effect,clusterSIG_start[c]:clusterSIG_stop[c],pm])) == max(np.sign(diffMeans[effect,clusterSIG_start[c]:clusterSIG_stop[c],pm])) ) # Sign of a cluster should be the same
            sign = np.sign(diffMeans[effect,clusterSIG_start[c]:clusterSIG_stop[c],pm])[0]
            if sign == -1: color=colors[0]
            if sign == 1: color=colors[1]
            # ax3.axhline(0.1,t_norm_beat[clusterSIG_start[c]] / (t_norm_beat.max()+1), t_norm_beat[clusterSIG_stop[c]] / (t_norm_beat.max()+1) , c=color, linewidth=13, alpha=0.7, label=label)
            ax3.axvspan(t_norm_beat[clusterSIG_start[c]], t_norm_beat[clusterSIG_stop[c]], facecolor=color, alpha=0.4)
            # ax3.text((t_norm_beat[clusterSIG_start[c]] + t_norm_beat[clusterSIG_stop[c]])/2, ymax*0.95, sign_str, horizontalalignment='center', fontsize=11)
            ax3.text(t_norm_beat[clusterSIG_start[c]], ymax/2, "p = " + str(p_corr[c]), rotation=90, verticalalignment='center',style='italic')
        if (m==0):    # for the legend
            ax3.axvspan(0,0,0 , facecolor=colors[0], alpha=0.6, label=label_high+' < '+label_low) 
            ax3.axvspan(0,0,0 , facecolor=colors[1], alpha=0.6, label=label_high+' > '+label_low) 
        ax3.axvspan(t_norm_beat[0], t_norm_beat[Nsamps_before_min], facecolor='lightgray',zorder=10,alpha=0.6); ax3.axvspan(t_norm_beat[-Nsamps_before_min], t_norm_beat[-1], facecolor='lightgray',zorder=10,alpha=0.6);
        
        ax3.tick_params(axis='y', labelsize=12)

        ax3.set_ylim(ymin,ymax); 
        ax3.set_xticks( np.arange(1,105,4) ) 
        ax3.set_xlim( ([min(t_norm_beat),max(t_norm_beat)]))
       
        if m == 0 : 
            x_labels = np.empty(np.arange(1,105,4).shape, dtype('<U21'))
            x_labels[1]='     SILENCE'; x_labels[-2]='       SILENCE';
            x_labels[5]+='DRUMS';  x_labels[9]+='+BASS';  x_labels[13]+='+KEYBOARD';  x_labels[17]+='+VOICE'; x_labels[21]+='+VOICE BIS';
            ax3.set_xticklabels( x_labels , fontsize=11 )
            ax3.tick_params(axis='x', which='both', length=0)
            ax3.xaxis.set_ticks_position('top')
            ax3.vlines(np.arange(13,94,16) , ymin,ymax*1.1, color='gray', alpha=0.6, zorder=20 , clip_on=False)
        elif m == 6 : 
            x_labels = (np.arange(1,105,4)//4 - 2).astype('<U21')
            x_labels[0] = x_labels[2] = x_labels[-2] = ''; x_labels[1]=''; x_labels[-1]=''
            ax3.set_xticklabels( x_labels , fontsize=11 ); ax3.set_xlabel('Bar', fontsize=11)
            ax3.vlines(np.arange(13,94,16) , ymin,ymax, color='gray', alpha=0.6, zorder=20)
        else : 
            ax3.set_xticks([])
            ax3.vlines(np.arange(13,94,16) , ymin,ymax, color='gray', alpha=0.6, zorder=20)

        # ax3.vlines(93+6.8 , ymin,ymax, color='gray', linestyles='dashed',zorder=20,alpha=0.5)
        
        if m == 0: ax3.legend(loc='upper left' , fontsize = 12).set_zorder(30)
    
    fig.tight_layout()
    fig.savefig(output_dir + '/SUPP_MAGICPLOT-QoM_' + bigtitle[effect] + '.png', dpi=600, bbox_inches='tight'); plt.close() 

