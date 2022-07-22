# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 20:58:27 2021

@author: AnalysisPC
"""

#-----------------------------------------------------------------------------
# Prophesee camera data filtering and display/export
#-----------------------------------------------------------------------------
# September 29th, 2021
# AUTHOR:
# Clément Cabriel
# Institut Langevin, ESPCI Paris / CNRS
# clement.cabriel@espci.fr , cabriel.clement@gmail.com
#-----------------------------------------------------------------------------
# Python 3.8.0
# Set the parameters below, then run the program
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# INPUT FORMAT
# Events .raw file or events .npy file: (columns) 'x': x position (in pix), 'y': y position (in pix), 'p': event polarity/algebric intensity, 't': timestamp (in us)
#-----------------------------------------------------------------------------
# OUTPUT FORMAT
# Events .npy file: (columns) 'x': x position (in pix), 'y': y position (in pix), 'p': event polarity/algebric intensity, 't': timestamp (in us)
# Frames .tif image stack
#-----------------------------------------------------------------------------

# import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import matplotlib
import imageio
from PIL import Image
from metavision_core.event_io.raw_reader import RawReader

#-----------------------------------------------------------------------------
# USER DEFINED PARAMETERS
#-----------------------------------------------------------------------------

# Import parameters
filepath='C:/absolute/path/to/data/AF647_coverslip.raw'
buffer_size=1e8                 # increase the buffer size if the input file is too large to be loaded. Decrease the buffer size if the PC is out of memory. Tested successfully with buffer size = 4e9 for 128 Gb of memory

# General filtering parameters
time_limits=[0,1]             # Time restriction of the dataset (in seconds) [t_min,t_max]
xy_limits='auto'                  # ROI restriction (in pixels) [x_min,y_min,x_max,y_max]. =None not to filter. ='auto' to automatically detect the limits
discard_up=False                # =True to discard positive events (rising edges)
discard_down=False              # =True to discard negative events (falling edges)

# Spatio-temporal filtering parameters
skip_filtering=True                 # =True to skip the filtering step; =False not to
space_gating=3                  # Half width of the space gating (in pix)
time_gating=50.                 # Half width of the time gating (in ms)
detection_threshold=4           # Minimum number of events per detection window

# Export parameters
export_events_data=False       # =True to export the (filtered) data; =False not to
export_frames=True              # =True to export the frames as a .tif stack; =False not to
xy_bin=1                        # Reduction ratio for the xy coordinates (set N to bin NxN pixels to 1 pixel)
time_bin=10.                    # Time bin used for the time sampling of the frames (in ms)
integrate_signal=False          # =True to integrate over time the events on the reconstructed maps. Can be used to generate camera-like intensity frames
sign_display="both"         # ="absolute" to display both positive and negative events, both counting as positive; ="both" to display both positive and negative events, each counting as their own sign; ="positive" to display only positive events; ="negative" to display only negative events
sum_all_frames=False             # =True to project all frames on one; =False not to

#-----------------------------------------------------------------------------
# END OF USER DEFINED PARAMETERS
#-----------------------------------------------------------------------------

# Initializations
time_limits_raw=[time_limits[0],time_limits[1]]
output_depth=16                 # =16 to export as 16-bit images; =8 to export as 8-bit images
if output_depth!=16 and output_depth!=8:
    output_depth=8
    print('Output image depth invalid; setting output_depth=8')
on_off_bias=0.                  # Compensation factor to normalise between positive and negative events, which may be acquired with different thresholds. =0 for no bias, in ]-1,0[ interval to favor negative events, in ]0,1[ interval to favor positive events
frame_size=[1280,720]            # Size of the sensor array (in pixels) [x,y]

# Data loading
print('Loading data...')
if filepath[-4:]=='.raw':
    record_raw = RawReader(filepath,max_events=int(buffer_size))
    sums = 0
    while not record_raw.is_done() and record_raw.current_event_index() < buffer_size:
        events = record_raw.load_delta_t(50000)
        sums += events.size
    record_raw.reset()
    events = record_raw.load_n_events(sums)
elif filepath[-4:]=='.npy':
    events=np.load(filepath)
print('Data loaded')
print('')

# Statistics display
if filepath[-4:]=='.npy' or filepath[-4:]=='.raw':
    print('Number of events:',len(events))
    msk=(events['p']==0)
    print('Number of positive events:',np.sum(msk==False))
    print('Number of negative events:',np.sum(msk))
    print('Fraction of negative events: ',np.sum(msk)*100.0/len(msk),'%')
    print('Total acquisition time: ',np.max(events['t'])/1e6,'s')
    print('')

# Basic filtering
time_limits=[time_limits[0],np.min([time_limits[1],np.max(events['t']+1)*1./1e6])]
events=events[(events['t']<time_limits[1]*1e6)*(events['t']>=time_limits[0]*1e6)]
events['t']-=int(np.floor(time_limits[0]*1e6))
if not xy_limits==None:
    if xy_limits=='auto':
        xy_limits=[np.min(events['x']),np.min(events['y']),np.max(events['x'])+1,np.max(events['y'])+1]
    events=events[(events['x']>=xy_limits[0])*(events['x']<xy_limits[2])*(events['y']>=xy_limits[1])*(events['y']<xy_limits[3])]
    events['x']-=xy_limits[0]
    events['y']-=xy_limits[1]
    frame_size=[xy_limits[2]-xy_limits[0],xy_limits[3]-xy_limits[1]]
if discard_up:events=events[(events['p']!=1)]
if discard_down:events=events[(events['p']!=0)]
print('Filtering done')
print('Number of events after filtering:',len(events))
print('')

# Spatio-temporal filtering
if not skip_filtering==True:
    print('Starting spatio-temporal filtering...')
    time_gating*=1000
    nb_slices=[1+int(np.ceil(np.max(events['t'])*1./time_gating)),1+int(np.ceil(np.max(events['x'])*1./space_gating)),1+int(np.ceil(np.max(events['y'])*1./space_gating))]
    events_slices=np.empty((nb_slices[0],nb_slices[1],nb_slices[2]),dtype=object)
    print('Slicing data...')
    for t in range(nb_slices[0]):
        print(t*100./nb_slices[0],'%')
        msk=(events['t']>=(t-1)*time_gating)*(events['t']<(t+2)*time_gating)
        sub_data=events[msk]
        for x in range(nb_slices[1]):
            msk=(sub_data['x']>=(x-1)*space_gating)*(sub_data['x']<(x+2)*space_gating)
            sub_sub_data=sub_data[msk]
            for y in range(nb_slices[2]):
                msk=(sub_sub_data['y']>=(y-1)*space_gating)*(sub_sub_data['y']<(y+2)*space_gating)
                sub_sub_sub_data=sub_sub_data[msk]
                events_slices[t,x,y]=sub_sub_sub_data
                
    valid_neighbors=np.zeros(len(events))
    print('Counting neighbours...')
    for k in np.arange(len(events)):
        if k%1000==0:print(k*100./len(events),'%')
        ind_k=[int(np.floor(events['t'][k]*1./time_gating)),int(np.floor(events['x'][k]*1./space_gating)),int(np.floor(events['y'][k]*1./space_gating))]
        sub_data=events_slices[ind_k[0],ind_k[1],ind_k[2]]
        msk=(sub_data['t']>=events['t'][k]-time_gating)*(sub_data['t']<=events['t'][k]+time_gating)
        msk*=(sub_data['x']>=events['x'][k]-space_gating)*(sub_data['x']<=events['x'][k]+space_gating)
        msk*=(sub_data['y']>=events['y'][k]-space_gating)*(sub_data['y']<=events['y'][k]+space_gating)
        msk_same_xy=(sub_data['x']==events['x'][k])*(sub_data['y']==events['y'][k])
        msk*=(msk_same_xy==False)
        valid_neighbors[k]=np.sum(msk)
    
    events=events[valid_neighbors>=detection_threshold]
    print('Number of events conserved',len(events),'/',np.shape(valid_neighbors)[0],'( '+str(len(events)*100./np.shape(valid_neighbors)[0])+" % )")
    print('Spatio-temporal filtering done')
    print('')
if export_events_data:
    if skip_filtering==False:
        name=filepath[:-4]+'_filtered_'+str(space_gating)+'-'+str(int(time_gating/1000))+'-'+str(detection_threshold)+'.npy'
    elif skip_filtering==True:
        name=filepath[:-4]+'_unfiltered.npy'
    np.save(name,events)
    print('Data saved')
    print('')

# Data processing
events['p']=events['p']*2-1.+on_off_bias
print('Data processing done')
print('')

# Frame generation
time_bin*=1e3
xy_bin=int(np.round(xy_bin))
if xy_bin<1:xy_bin=1;print('Invalid xy binning; settting xy_bin=1')
if not xy_bin==1:
    frame_size=[int(np.ceil(frame_size[0]/xy_bin)),int(np.ceil(frame_size[1]/xy_bin))]
    events['x']=np.floor(events['x']/xy_bin).astype(int)
    events['y']=np.floor(events['y']/xy_bin).astype(int)
nb_frames=int(np.ceil(np.max(events['t'])/time_bin))
if sum_all_frames:nb_frames=1;time_bin=1e12
frames=np.zeros((frame_size[1],frame_size[0],nb_frames))
print('Generating frames...')
for k in np.arange(nb_frames):
    if k%250==0:print(k*100./np.max(np.arange(nb_frames)))
    time_limits_k=[int(np.floor(k*time_bin)),int(np.floor((k+1)*time_bin))]
    msk=(events['t']>=time_limits_k[0])*(events['t']<time_limits_k[1])
    sub_data=events[msk]
    if sign_display=='absolute':sub_data['p']=np.abs(sub_data['p'])
    elif sign_display=='positive':sub_data=sub_data[sub_data['p']>0]
    elif sign_display=='negative':sub_data=sub_data[sub_data['p']<0];sub_data['p']=np.abs(sub_data['p'])
    for l in range(len(sub_data)):
        if integrate_signal:
            frames[sub_data['y'][l],sub_data['x'][l],k:]+=sub_data['p'][l]
        else:
            frames[sub_data['y'][l],sub_data['x'][l],k]+=sub_data['p'][l]
if export_frames:
    if sum_all_frames==False and sign_display=='both' and output_depth==8:imageio.mimwrite(filepath[:-4]+'_frames_'+str(time_bin/1000)+'ms-'+str(xy_bin)+'pix_'+sign_display+'.tif', np.swapaxes(np.swapaxes(frames+127,0,2),1,2).astype('uint8'))
    elif sum_all_frames==False and sign_display!='both' and output_depth==8:imageio.mimwrite(filepath[:-4]+'_frames_'+str(time_bin/1000)+'ms-'+str(xy_bin)+'pix_'+sign_display+'.tif', np.swapaxes(np.swapaxes(frames,0,2),1,2).astype('uint8'))
    elif sum_all_frames==False and output_depth==16:imageio.mimwrite(filepath[:-4]+'_frames_'+str(time_bin/1000)+'ms-'+str(xy_bin)+'pix_'+sign_display+'.tif', np.swapaxes(np.swapaxes(frames,0,2),1,2).astype('int16'))
    else:Image.fromarray(frames[:,:,0]).save(filepath[:-4]+'_frames_sum_'+sign_display+'.tif')
print('Frames generation done')
print('')
