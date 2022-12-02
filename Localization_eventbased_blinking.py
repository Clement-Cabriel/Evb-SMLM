# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 19:10:13 2022

@author: Clement
"""

#-----------------------------------------------------------------------------
# Prophesee camera data processing for blinking analysis
#-----------------------------------------------------------------------------
# January 12th, 2022
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
# Events .raw file or events .npy file: (columns) 'x': x position (in pix), 'y': y position (in pix), 'p': event polarity, 't': timestamp (in us)
#-----------------------------------------------------------------------------
# OUTPUT FORMAT
# Events .npy file: (columns) [id,t,y,x,nb_up,nb_down,t_up,t_down] with all times in us and all distances in object plane pixels
#-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.ndimage
import scipy.signal
from scipy import optimize
import heapq
import time
from datetime import datetime
import gc
from joblib import Parallel, delayed
import multiprocessing

#-----------------------------------------------------------------------------
# USER DEFINED PARAMETERS
#-----------------------------------------------------------------------------

# Import parameters
filepath='C:/absolute/path/to/data/AF647_tubulin.npy'
buffer_size=1e9                 # increase the buffer size if the input file is too large to be loaded. Decrease the buffer size if the PC is out of memory. Tested successfully with buffer size = 4e9 for 128 Gb of memory

# Optical parameters
pixel_size=67.                  # Size of the pixel in the object plane in nm

# Computation parameters
multithread=True                # =True to use multithread parallelization; =False not to. Note: multithread parallelization increases calculation speed

# General filtering parameters
time_limits=[0.,1e9]             # Time restriction of the dataset (in seconds) [t_min,t_max]
xy_limits=None                  # ROI restriction (in pixels) [x_min,y_min,x_max,y_max]. =None not to filter. ='auto' to automatically detect the limits

# PSF detection parameters
threshold_detection=3.         # Wavelet detection threshold
exclusion_radius=4.              # Radius of the exclusion area (if two or more PSFs are closer than twice the value, they will all be discarded) (in pixels)
min_diameter=1.25               # Minimum radius of the thresholded area (in pixels)
max_diameter=4.                 # Maximum radius of the thresholded area (in pixels)
time_bin_frames=20              # Time bin (in ms) of the frames

# Localization parameters
area_radius=4                   # Radius of the fitting area (in pixels)
time_area_limits=[60.,200.]      # Relative time limits (in ms) of the events to consider within an ROI. =[t_minus,t_plus]. The events considered for the localization are those within [t0-t_minus,t0+t_plus[ where t0 is the detected time of the ROI
localization_method='COM'  # Method for the PSF localization. ='COM' for center of mass calculation, ='Gaussian' for Gaussian fitting

# Display parameters
display_frames=True             # =True to display the frames, =False not to
number_frames_display=0       # Number of frames to display

# Output parameters
export_localized_results=True   # =True to save the localized results; =False not to

#-----------------------------------------------------------------------------
# END OF USER DEFINED PARAMETERS
#-----------------------------------------------------------------------------

# Initializations - general
start_time=datetime.now()
font = {'size'   : 16,
        'weight' : 'normal'}
mpl.rc('font', **font)
col=['b','r','g','k','c','m','y']
time_bin_frames*=1000
t_min_display=0.
time_area_limits=[time_area_limits[0]*1000,time_area_limits[1]*1000]
batch_size_PSFdetection=50000       # Slice sizes when loading the events data on the fly
batch_size_localization=50000
if number_frames_display<1:display_frames=False
if display_frames:multithread=False;export_localized_results=False
if multithread==True:num_cores = multiprocessing.cpu_count()
else:num_cores = 1

# Initializations - wavelet
kernel1=np.array([0.0625,0.25,0.375,0.25,0.0625])
kernel2=np.array([0.0625,0,0.375,0,0.25,0,0.0625])
kernel_size=8
kernel=np.ones((kernel_size,kernel_size))/(kernel_size**2.)

# Initializations - localization
initial_guess_w=150./pixel_size
yrange=np.arange(2*area_radius+1)*1.
xrange=np.arange(2*area_radius+1)*1.

# Logs
def Save_parameters_init(events,start_time):
    
    parameters_dictionary={"filepath":filepath,
    "pixel_size":pixel_size,
    "buffer_size":buffer_size,"multithread":multithread,
    "time_limits":time_limits,"xy_limits":xy_limits,
    "threshold_detection":threshold_detection,"exclusion_radius":exclusion_radius,"min_diameter":min_diameter,"max_diameter":max_diameter,
    "time_bin_frames":time_bin_frames,"area_radius":area_radius,"time_area_limits":time_area_limits,"localization_method":localization_method,
    "display_frames":display_frames,"number_frames_display":number_frames_display,
    "export_localized_results":export_localized_results}
    
    start_time=start_time.strftime("%d/%m/%Y %H:%M:%S")
    run_parameters_dictionary={"Start time":start_time,"Min_t":np.min(events['t']),"Max_t":np.max(events['t']),"Min_y":np.min(events['y']),"Max_y":np.max(events['y']),"Min_x":np.min(events['x']),"Max_x":np.max(events['x'])}
    
    with open(filepath[:-4]+'_log.txt', 'w') as f:
        print(__file__,file=f)
        print('',file=f)
        print('PARAMETERS',file=f)
        for k in list(parameters_dictionary.keys()):
            print('        '+k+'='+str(parameters_dictionary[k]), file=f)
    with open(filepath[:-4]+'_log.txt', 'a') as f:
        print('',file=f)
        print('INITIALIZATION LOG',file=f)
        for k in list(run_parameters_dictionary.keys()):
            print(k+'='+str(run_parameters_dictionary[k]), file=f)
    
def Save_parameters_final():
    
    end_time=datetime.now()
    end_time=end_time.strftime("%d/%m/%Y %H:%M:%S")
    final_parameters_dictionary={"End time":end_time}
    
    with open(filepath[:-4]+'_log.txt', 'a') as f:
        print('',file=f)
        print('COMPLETION LOG',file=f)
        for k in list(final_parameters_dictionary.keys()):
            print(k+'='+str(final_parameters_dictionary[k]), file=f)
        print('Done',file=f)

# Frames generation
def generate_single_frame(events,times):
    # frame=np.zeros((frame_size[0],frame_size[1]),dtype=float)
    frame=np.zeros((int(np.ceil(frame_size[0])),int(np.ceil(frame_size[1]))))
    msk=(events['t']>=times[0])*(events['t']<times[1])
    sub_events=events[msk]
    for k in np.arange(len(sub_events)):
        xycoords=[int(np.floor(sub_events['y'][k])),int(np.floor(sub_events['x'][k]))]
        frame[xycoords[0],xycoords[1]]+=1
    return frame
    
# PSF detection
def wavelet_detection(frame):
    V1=scipy.ndimage.convolve1d(scipy.ndimage.convolve1d(frame,kernel1,axis=1),kernel1,axis=0)
    V2=scipy.ndimage.convolve1d(scipy.ndimage.convolve1d(V1,kernel2,axis=1),kernel2,axis=0)
    Wavelet2nd=V1-V2
    Wavelet2nd-=scipy.ndimage.convolve(Wavelet2nd,kernel)
    Wavelet2nd*=(Wavelet2nd>=0)
    image_to_label=(Wavelet2nd>=threshold_detection*np.std(Wavelet2nd))*Wavelet2nd
    return image_to_label

def detect_PSFs(frame):# verify exceptions for zero detections
    
    # Generate image to label
    image_to_label=wavelet_detection(frame)
    
    # Create labels
    labels,nb_labels=scipy.ndimage.label(image_to_label)
    label_list=np.arange(nb_labels)+1
    area_sizes=scipy.ndimage.sum((image_to_label>0),labels,index=label_list)
    
    # Filter labels by size
    msk=(area_sizes>(min_diameter*2+1)**2)*(area_sizes<(max_diameter*2+1)**2)
    label_list=label_list[msk]
    labels*=np.isin(labels,label_list)
    nb_labels=len(label_list)
    
    # Calculate center of mass
    coordinates_COM=scipy.ndimage.center_of_mass(image_to_label,labels,label_list)
    coordinates_COM=np.asarray(coordinates_COM)
    if len(np.shape(coordinates_COM))==1:
        if np.shape(coordinates_COM)[0]==0:
            coordinates_COM=np.zeros((0,2))
    
    # Filter coordinates by distance to each other
    msk=np.ones(nb_labels,dtype=bool)
    mgy1,mgy2=np.meshgrid(coordinates_COM[:,0],coordinates_COM[:,0].transpose())
    mgx1,mgx2=np.meshgrid(coordinates_COM[:,1],coordinates_COM[:,1].transpose())
    msk2=(np.abs(mgy1-mgy2)<2*exclusion_radius)*(np.abs(mgx1-mgx2)<2*exclusion_radius)
    np.fill_diagonal(msk2,False)
    ind=np.nonzero(msk2)
    ind=ind[0].tolist()+ind[1].tolist()
    ind=list(set(ind))
    msk[ind]=False
    label_list=label_list[msk]
    labels*=np.isin(labels,label_list)
    nb_labels=len(label_list)
    coordinates_COM=coordinates_COM[msk]
    
    # Filter coordinates by distance to the frame edges
    if nb_labels>0:
        msk=(coordinates_COM[:,0]>area_radius*2)*(coordinates_COM[:,1]>area_radius*2)*(coordinates_COM[:,0]<np.shape(frame)[0]-area_radius*2)*(coordinates_COM[:,1]<np.shape(frame)[1]-area_radius*2)
        label_list=label_list[msk]
        labels*=np.isin(labels,label_list)
        nb_labels=len(label_list)
        coordinates_COM=coordinates_COM[msk]
        
    if len(np.shape(coordinates_COM))==1:
        coordinates_COM=np.zeros((0,2))
        
    x=coordinates_COM[:,1]
    y=coordinates_COM[:,0]
    ROIs=np.zeros((np.shape(x)[0],3))
    ROIs[:,0]=y
    ROIs[:,1]=x
    
    return ROIs

def process_frame(frame,times):    
    ROIs=detect_PSFs(frame)
    ROIs[:,2]=times[0]
    return ROIs

# Display
def display_frame_PSFs(frame,frame_number):    
    plt.figure('Frame '+str(frame_number))
    figManager=plt.get_current_fig_manager()
    figManager.window.showMaximized()
    
    ROIs=detect_PSFs(frame)
    sbp=plt.subplot(1,1,1)
    sbp.imshow(frame,cmap='hot',interpolation='none')
    plt.plot(ROIs[:,1],ROIs[:,0],'sb',ms=12,mew=1.5,markerfacecolor='none')
    plt.gca().set_aspect(1.)
    
# Localization
def process_ROI(events,ROI):
    msk=(events['t']>=ROI[2]-time_area_limits[0])*(events['t']<ROI[2]+time_area_limits[1])
    rho=(events['y']-ROI[0])**2+(events['x']-ROI[1])**2
    msk*=(rho<=area_radius**2)
    sub_events=events[msk]
    localized_data=localization_PSF(sub_events,ROI)
    events_ROI=np.array([])
    return localized_data,events_ROI
    
def localization_PSF(sub_events,ROI):
    if localization_method=='COM':
        y=np.mean(sub_events['y'])
        x=np.mean(sub_events['x'])
    elif localization_method=='Gaussian':
        edge=[int(np.round(ROI[0]))-area_radius,int(np.round(ROI[1]))-area_radius]
        opt,err=gaussian_fitting(sub_events,edge)
        y=opt[0]+edge[0]
        x=opt[1]+edge[1]
    t=np.mean(sub_events['t'])
    msk=(sub_events['p']==1)
    nb_up=np.sum(msk)
    t_up=np.mean(sub_events['t'][msk])
    msk=(sub_events['p']==0)
    nb_down=np.sum(msk)
    t_down=np.mean(sub_events['t'][msk])
    return np.array([t,y,x,nb_up,nb_down,t_up,t_down])

def gaussian_fitting(sub_events,edge):
    sub_image=np.zeros((area_radius*2+1,area_radius*2+1))
    for k in np.arange(len(sub_events)):
        sub_image[sub_events['y'][k]-edge[0],sub_events['x'][k]-edge[1]]+=1
    initial_guess=area_radius,area_radius,initial_guess_w,np.max(sub_image),np.median(sub_image)
    errorfunction=lambda p:np.ravel(gaussian_symmetrical(*p)-sub_image)
    p_full = optimize.leastsq(errorfunction,initial_guess,gtol=1e-4,ftol=1e-4,full_output=True)
    p=p_full[0]
    fv=p_full[2]['fvec']
    ss_err=(fv**2).sum()
    ss_tot=((sub_image-sub_image.mean())**2).sum()
    err=ss_err/ss_tot    
    return p,err

def gaussian_symmetrical(y0,x0,width,height,offset):
    fX=np.exp(-(xrange-x0)**2/(2.*width**2))
    fY=np.exp(-(yrange-y0)**2/(2.*width**2))
    fY=fY.reshape(len(fY),1)
    return offset+height*fY*fX

def load_events(ev_list,ev_loaded,sl_size,ind_loaded):
    ind_max=np.min([len(ev_list),ind_loaded[1]+sl_size])
    ev_loaded=np.append(ev_loaded,ev_list[ind_loaded[1]:ind_max])
    return ev_loaded,[ind_loaded[0],ind_max]

def unload_events(ev_loaded,t_min,ind_loaded):
    msk=(ev_loaded['t']>=t_min)
    if np.sum(msk)==0:msk[-1]=True
    ev_loaded=ev_loaded[msk]
    ind_min=np.argmax(msk)
    return ev_loaded,[ind_loaded[0]+ind_min,ind_loaded[1]]

# Data loading
print('Loading data...')
if filepath[-4:]=='.raw':
    from metavision_core.event_io.raw_reader import RawReader
    record_raw = RawReader(filepath,max_events=int(buffer_size))
    sums = 0
    while not record_raw.is_done() and record_raw.current_event_index() < buffer_size:
        events = record_raw.load_delta_t(50000)
        sums += events.size
    record_raw.reset()
    events = record_raw.load_n_events(sums)
elif filepath[-4:]=='.npy':
    events=np.load(filepath)
nb_ev_init=len(events)
print('Data loaded')
print('')

# Statistics display
print('Number of events:',len(events))
if filepath[-4:]=='.npy':
    msk=(events['p']==0)
    print('Number of positive events:',np.sum(msk==False))
    print('Number of negative events:',np.sum(msk))
    print('Fraction of negative events: ',np.sum(msk)*100.0/len(events),'%')
    print('Total acquisition time: ',np.max(events['t'])/1e6,'s')
    print('')
    
# Basic filtering
print('Filtering data...')
events=events[(events['t']<time_limits[1]*1e6)*(events['t']>=time_limits[0]*1e6)]
if not xy_limits==None:
    if xy_limits=='auto':
        xy_limits=[np.min(events['x']),np.min(events['y']),np.max(events['x'])+1,np.max(events['y'])+1]
    events=events[(events['x']>=xy_limits[0])*(events['x']<xy_limits[2])*(events['y']>=xy_limits[1])*(events['y']<xy_limits[3])]
    events['x']-=xy_limits[0]
    events['y']-=xy_limits[1]
time_max=np.max(events['t'])
frame_size=[np.max(events['y'])+1,np.max(events['x'])+1]
if export_localized_results:
    Save_parameters_init(events,start_time)
time.sleep(0.5)
gc.collect()
print('Filtering done')
print('Number of events after filtering:',len(events))
print('')
    

print('Checking data sorting...')
diff=events['t']-np.roll(events['t'],1)
if np.min(diff[1:])<0:
    print('Sorting data...')
    events=np.sort(events,order='t')
print('Data sorted')
print('')
    
def slice_data(events,nb_slices):
    
    slice_size=1.*len(events)/nb_slices
    slice_size=np.int64(np.ceil(slice_size))
    data_split=[]
    for k in np.arange(nb_slices):
        ind=[np.compat.long(k*slice_size),np.compat.long((k+1)*slice_size)]
        data_split.append(events[ind[0]:ind[1]])
    return data_split
                
def compute_thread(i,sub_events):
    
    print('Detecting PSFs (thread '+str(i)+')...')
    time_max=np.max(sub_events['t'])
    time_min=np.min(sub_events['t'])
    nb_frames=int(np.ceil((time_max-time_min)/time_bin_frames))
    cnt=0
    ROIs=[]
    events_loaded=sub_events[0:2]
    indices_loaded=[0,2]
    for k in np.arange(nb_frames):
        times=[time_min+k*time_bin_frames,time_min+(k+1)*time_bin_frames]
        # Display
        if display_frames:
            if times[0]<t_min_display or cnt==number_frames_display:
                pass
            else:
                frame=generate_single_frame(sub_events,times)
                if display_frames:display_frame_PSFs(frame,k)
                ROIs.append(process_frame(frame,times))
                cnt=cnt+1
        # Processing
        else:
            if k%100==0:print('Detecting PSFs (thread '+str(i)+'): '+str(k)+' / '+str(nb_frames)+' ( '+str(int(np.round(100.*k/nb_frames)))+' % )');
            # Load data slices on the fly
            cnt_load=0
            while events_loaded['t'][-1]<=times[1]:
                if indices_loaded[1]>=len(sub_events)-1:break
                else:
                    events_loaded,indices_loaded=load_events(sub_events,events_loaded,batch_size_PSFdetection,indices_loaded)
                    if cnt_load==0:
                        events_loaded,indices_loaded=unload_events(events_loaded,times[0],indices_loaded)
                    cnt_load+=1
            # Detect PSFs
            frame=generate_single_frame(events_loaded,times)
            ROIs.append(process_frame(frame,times))
    ROIs=np.vstack(ROIs)
    time.sleep(0.5)
    gc.collect()
    
    print('Localizing (thread '+str(i)+')...')
    localized_data,deprecated_array=process_ROI(sub_events,ROIs[0,:])
    localized_data=np.zeros((np.shape(ROIs)[0],np.shape(localized_data)[0]))
    deprecated_array=np.array([])
    if display_frames:sub_events=sub_events[sub_events['t']>=t_min_display]
    cnt=0
    events_loaded=sub_events[0:2]
    indices_loaded=[0,2]
    for k in np.arange(np.shape(ROIs)[0]):
        # if number_frames_display>0:break
        if k%1000==0:print('Localizing (thread '+str(i)+'): '+str(k)+' / '+str(np.shape(ROIs)[0])+' ( '+str(int(np.round(100.*k/np.shape(ROIs)[0])))+' % )')
        # Loading data slices on the fly
        cnt_load=0
        while events_loaded['t'][-1]<=ROIs[k,2]+time_bin_frames+time_area_limits[1]:
            if indices_loaded[1]>=len(sub_events)-1:break
            else:
                events_loaded,indices_loaded=load_events(sub_events,events_loaded,batch_size_localization,indices_loaded)
                if cnt_load==0:
                    events_loaded,indices_loaded=unload_events(events_loaded,ROIs[k,2]-time_area_limits[0],indices_loaded)
                cnt_load+=1
        # Localization
        localized_data[k,:],deprecated_array_k=process_ROI(events_loaded,ROIs[k,:])
    
    return localized_data,deprecated_array
    
print('Detecting PSFs and localizing...')
events_split=slice_data(events,num_cores)
RES = Parallel(n_jobs=num_cores,backend="loky")(delayed(compute_thread)(i,events_split[i]) for i in range(len(events_split)))
localized_data=[]
for i in np.arange(np.shape(RES)[0]):
    localized_data.append(RES[i][0])
localized_data=np.vstack(localized_data)
msk=(localized_data[:,2]>=0)*(localized_data[:,1]>=0)*(localized_data[:,2]<frame_size[1])*(localized_data[:,1]<frame_size[0])
localized_data=localized_data[msk,:]
time.sleep(0.5)
gc.collect()
print('PSF detection and localization done')
print('Number of molecules detected:',np.shape(localized_data)[0])
print('')

if export_localized_results:
    print('Saving results...')
    localized_data_out=np.zeros(np.shape(localized_data)[0],dtype=[('id',np.int64), ('t',np.float32), ('y',np.float32), ('x',np.float32), ('N_up', np.int16), ('N_down', np.int16), ('t_up', np.float32), ('t_down', np.float32)])
    localized_data_out['id']=np.arange(np.shape(localized_data)[0])
    localized_data_out['t']=localized_data[:,0]
    localized_data_out['y']=localized_data[:,1]
    localized_data_out['x']=localized_data[:,2]
    localized_data_out['N_up']=localized_data[:,3]
    localized_data_out['N_down']=localized_data[:,4]
    localized_data_out['t_up']=localized_data[:,5]
    localized_data_out['t_down']=localized_data[:,6]
    np.save(filepath[:-4]+'_localized.npy',localized_data_out)
    Save_parameters_final()
    print('Results saved')
    print('')
            
print('')
print('Done')
print('')
