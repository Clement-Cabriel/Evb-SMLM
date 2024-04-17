# -*- coding: utf-8 -*-
"""
Created on Thu Jan 20 21:26:01 2022

@author: Clement
"""

#-----------------------------------------------------------------------------
# 2D Direct Cross-Correlation image-based drift correction
#-----------------------------------------------------------------------------
# January 20th, 2022
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
# Localization .npy file. Required fields: 'y','x','t' in arbitrary units. Note: t should start at 0
#-----------------------------------------------------------------------------
# OUTPUT FORMAT
# Drift-corrected .npy file: same format as input, same units as input
# Drift .npy file: (columns) y,x,t with y,x in input pixels and t in input time unit
#-----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.ndimage
from scipy.optimize import leastsq

class DCC_workflow():
    
    def __init__(self):
        
        #-----------------------------------------------------------------------------
        # USER DEFINED PARAMETERS
        #-----------------------------------------------------------------------------
        
        # Acquisition parameters
        self.filepath='//oaxaca/Data Backup/Clement/2022-09-13/COS7 alpha-tub-DNAPAINT-Atto655 permeab 0.5% Nikon EMCCD-gen4 8nM/B_gen4_SMLM60%_62-105/recording_2022-09-13T18-17-36.106Z_pos_localized.npy'
        self.filepath='Tubulin_DCC_test.npy'
        self.initial_pixel=80.              # Localization pixel size (in nm). Set =1 if the coordinates are already in nm
        
        # Correction parameters
        self.correction_pixel=50            # Pixel size of the correction (in nm)
        self.smoothing_frame=1.             # Sigma (in zoomed correction_pixel coordinates) of the Gaussian filtering of the frames
        self.smoothing_correlation=1.       # Sigma (in zoomed correction_pixel coordinates) of the Gaussian filtering of the correlation functions (recommanded with maxima detection)
        self.frame_ref_t=[166e6,333e6]      # Reference slice t limits (in us). =[t_min,t_max]
        self.slice_size_t=150e6             # Slice width in t (in us)
        self.max_shift=1000.                # Max value of the allowed shift relative to the reference image (in nm)
        self.ROI=[0,0,1e9,1e9]              # ROI limits (in initial localization pixels). =[x_min,y_min,x_max,y_max]
        
        # Computation parameters
        self.maximum_detection_method='Gaussian'              # maximum_detection_method for the peak detection. ='Gaussian' or ='Max' (maximum detection, i.e. taking the highest pixel). Gaussian yields more precise results, but the fitting might fail in some cases
        self.width_ROI_gaussian=100.        # Half window of the fitting around the maximum value (in nm)
        
        # Display parameters
        self.display_ref_frame=True         # =True to display the reference frame
        self.display_frames=True            # =True to display the sliced frames
        self.display_correlation=True       # =True to display the cross-correlation
        self.display_drift_curves=True      # =True to display the raw and interpolated drift curves
        
        # Output parameters
        self.interpolation_step=1e7         # Interpolation time step (in us) for the raw drift curve interpolation
        self.export_results=True            # =True to export the results
         
        #-----------------------------------------------------------------------------
        # END OF USER DEFINED PARAMETERS
        #-----------------------------------------------------------------------------
        
    def print_progress(self,current_index,max_index,period,indent):
        if current_index%period==0:
            print('    '*indent,current_index,'/',max_index,'(',"%0.1f"%(100.*current_index/max_index),'% )')

    # Frames generation
    def generate_single_frame(self,localizations):
        localizations=np.floor(localizations).astype(int)
        frame=np.zeros((self.frame_size[0],self.frame_size[1]))
        for k in np.arange(np.shape(localizations)[0]):
            frame[localizations[k,0],localizations[k,1]]+=1
        frame=scipy.ndimage.gaussian_filter(frame,self.smoothing_frame)
        return frame
    
    # Gaussian fit
    def GaussianFit2D(self,correlation):

        yv=np.arange(np.shape(correlation)[0])
        xv=np.arange(np.shape(correlation)[1])
        X,Y=np.meshgrid(yv,xv)
        position_max_init=[np.unravel_index(correlation.argmax(), correlation.shape)[0],np.unravel_index(correlation.argmax(), correlation.shape)[1]]
        
        initial_guess=np.zeros(5)
        self.fixed_param_fit=np.zeros(2)
        initial_guess[0]=position_max_init[0]
        initial_guess[1]=position_max_init[1]
        initial_guess[2]=2.*50./self.correction_pixel
        self.fixed_param_fit[0]=np.max(correlation)-np.min(correlation)
        self.fixed_param_fit[1]=np.min(correlation)
        
        ind_lim=[position_max_init[0]-self.width_ROI_gaussian,position_max_init[0]+self.width_ROI_gaussian+1,position_max_init[1]-self.width_ROI_gaussian,position_max_init[1]+self.width_ROI_gaussian+1]
        ind_lim=[np.max([0,ind_lim[0]]),np.min([np.shape(correlation)[0],ind_lim[1]]),np.max([0,ind_lim[2]]),np.min([np.shape(correlation)[1],ind_lim[3]])]
        X=X[ind_lim[0]:ind_lim[1],ind_lim[2]:ind_lim[3]]
        Y=Y[ind_lim[0]:ind_lim[1],ind_lim[2]:ind_lim[3]]
        correlation=correlation[ind_lim[0]:ind_lim[1],ind_lim[2]:ind_lim[3]]
        
        ErrorFunc = lambda param, y,x, data: (self.Gaussmodel2D(y,x, param).ravel() - data.ravel())
        [opt,norm] = leastsq(ErrorFunc, initial_guess, args=(Y,X, correlation), maxfev=500)
        position_max=[opt[0],opt[1]]
        
        if (position_max[0]-position_max_init[0])**2.+(position_max[1]-position_max_init[1])**2.>=0.5*np.sqrt(2.):position_max=position_max_init

        return position_max[0],position_max[1]
    
    def Gaussmodel2D(self,Y,X,a):
        F=self.fixed_param_fit[0]*np.exp(-((Y-a[0])**2)/(a[2]**2)-((X-a[1])**2)/(a[2]**2)) + self.fixed_param_fit[1]
        return F.ravel()
        
    # Maxima detection
    def MaximaDetection2D(self,correlation):
        ymax=np.unravel_index(correlation.argmax(), correlation.shape)[0]
        xmax=np.unravel_index(correlation.argmax(), correlation.shape)[1]
        return ymax,xmax
    
    def format_data(self,data_in):
        data_out=np.zeros((len(data_in),3))
        data_out[:,0]=data_in['y']
        data_out[:,1]=data_in['x']
        data_out[:,2]=data_in['t']
        return data_out
    
    def Main(self):
        
        # Data loading
        print('Loading data...')
        self.data_init=np.load(self.filepath)
        self.data=self.format_data(self.data_init)
        print('Data loaded')
        print('')
        
        # ROI filtering
        msk=(self.data[:,0]>=self.ROI[1])*(self.data[:,0]<self.ROI[3])*(self.data[:,1]>=self.ROI[0])*(self.data[:,1]<self.ROI[2])
        self.data=self.data[msk,:]
        
        # Initializations
        self.magnification_factor=self.initial_pixel/self.correction_pixel
        self.data[:,[0,1]]*=self.magnification_factor
        self.data[:,0]-=np.min(self.data[:,0])
        self.data[:,1]-=np.min(self.data[:,1])
        self.frame_size=np.array([np.max(self.data[:,0]),np.max(self.data[:,1])])
        self.frame_size=np.ceil(self.frame_size).astype(int)
        self.max_shift=int(np.ceil(self.max_shift/self.correction_pixel))
        self.nb_slices=int(np.ceil(np.max(self.data[:,2])/self.slice_size_t))
        print('Maximum shift allowed in correction pixels:',self.max_shift)
        print('Number of slices:',self.nb_slices)
        self.width_ROI_gaussian=int(np.ceil(self.width_ROI_gaussian/self.correction_pixel))
        self.xrange=np.arange(-self.max_shift,self.max_shift)
        
        # Generating initial frames
        msk=(self.data[:,2]>=self.frame_ref_t[0])*(self.data[:,2]<self.frame_ref_t[1])
        self.frame_ref=self.generate_single_frame(self.data[msk,:2])
        self.frames=np.zeros((self.frame_size[0],self.frame_size[1],self.nb_slices))
        for k in np.arange(self.nb_slices):
            msk=(self.data[:,2]>=k*self.slice_size_t)*(self.data[:,2]<(k+1)*self.slice_size_t)
            self.frames[:,:,k]=self.generate_single_frame(self.data[msk,:2])
        if self.display_ref_frame:
            plt.figure('Reference frame')
            figManager=plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.imshow(self.frame_ref,cmap='jet')
                    
        # Auto-correlation calculation
        autocorr=np.zeros((np.shape(self.xrange)[0],np.shape(self.xrange)[0]))
        for y in np.arange(np.shape(self.xrange)[0]):
            dy=self.xrange[y]
            for x in np.arange(np.shape(self.xrange)[0]):
                dx=self.xrange[x]
                translated_frame=np.roll(self.frame_ref,dy,axis=0)
                translated_frame=np.roll(translated_frame,dx,axis=1)
                autocorr[y,x]=np.sum(translated_frame*self.frame_ref)
        if self.smoothing_correlation>0:autocorr=scipy.ndimage.gaussian_filter(autocorr,self.smoothing_correlation)
        if self.display_correlation:
            plt.figure('Auto-correlation ( frame '+str(self.frame_ref_t[0])+' )')
            figManager=plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.imshow(autocorr,cmap='jet')
        if self.maximum_detection_method=='Gaussian':
            y0,x0=self.GaussianFit2D(autocorr)
        elif self.maximum_detection_method=='Max':
            y0,x0=self.MaximaDetection2D(autocorr)
        
        # Cross-correlation calculation
        print('')
        print('Measuring drift...')
        crosscorr=np.zeros((np.shape(self.xrange)[0],np.shape(self.xrange)[0],self.nb_slices))
        for k in np.arange(self.nb_slices):
            # print('Slice '+str(k+1)+' / '+str(self.nb_slices)+' ( '+str(100.*k/self.nb_slices)+' % )')
            print('    '+'Slice '+str(k+1)+' / '+str(self.nb_slices))
            self.frame_ref_crosscorr=self.frame_ref
            for y in np.arange(np.shape(self.xrange)[0]):
                self.print_progress(k*np.shape(self.xrange)[0]+y,self.nb_slices*np.shape(self.xrange)[0],1,2)
                dy=self.xrange[y]
                for x in np.arange(np.shape(self.xrange)[0]):
                    dx=self.xrange[x]
                    translated_frame=np.roll(self.frames[:,:,k],dy,axis=0)
                    translated_frame=np.roll(translated_frame,dx,axis=1)
                    crosscorr[y,x,k]=np.sum(translated_frame*self.frame_ref_crosscorr)
            if self.smoothing_correlation>0:crosscorr[:,:,k]=scipy.ndimage.gaussian_filter(crosscorr[:,:,k],self.smoothing_correlation)
            if self.display_correlation:
                plt.figure('Cross-correlation - Slice '+str(k)+' ( frame '+str(k*self.slice_size_t)+' )')
                figManager=plt.get_current_fig_manager()
                figManager.window.showMaximized()
                plt.imshow(crosscorr[:,:,k],cmap='jet')
            if self.display_frames:
                plt.figure('Frame '+str(k)+' ( '+str(k*self.slice_size_t)+' )')
                figManager=plt.get_current_fig_manager()
                figManager.window.showMaximized()
                plt.imshow(self.frames[:,:,k],cmap='jet')
                
        # Fitting
        drift_raw=np.zeros((self.nb_slices,3))
        for k in np.arange(self.nb_slices):
            drift_raw[k,2]=(k+0.5)*self.slice_size_t
            if self.maximum_detection_method=='Gaussian':
                y,x=self.GaussianFit2D(crosscorr[:,:,k])
            elif self.maximum_detection_method=='Max':
                y,x=self.MaximaDetection2D(crosscorr[:,:,k])
            drift_raw[k,0]=y-y0
            drift_raw[k,1]=x-x0
            if np.max(crosscorr[:,:,k])==0:
                drift_raw[k,:]=drift_raw[k-1,:]
        drift_interp=np.zeros((int(np.ceil(np.max(self.data[:,2]/self.interpolation_step)))+1,np.shape(drift_raw)[1]))
        drift_interp[:,2]=np.arange(np.shape(drift_interp)[0])*self.interpolation_step
        for k in np.arange(np.shape(drift_interp)[0]):
            if drift_interp[k,2]<=drift_raw[0,2]:
                drift_interp[k,0]=drift_raw[0,0]
                drift_interp[k,1]=drift_raw[0,1]
            elif drift_interp[k,2]>=drift_raw[-1,2]:
                drift_interp[k,0]=drift_raw[-1,0]
                drift_interp[k,1]=drift_raw[-1,1]
            else:
                ind=np.argmax((drift_raw[:,2]-drift_interp[k,2]>0))-1
                slope=(drift_raw[ind+1,0]-drift_raw[ind,0])/(drift_raw[ind+1,2]-drift_raw[ind,2])
                drift_interp[k,0]=drift_raw[ind,0]+slope*(drift_interp[k,2]-drift_raw[ind,2])
                slope=(drift_raw[ind+1,1]-drift_raw[ind,1])/(drift_raw[ind+1,2]-drift_raw[ind,2])
                drift_interp[k,1]=drift_raw[ind,1]+slope*(drift_interp[k,2]-drift_raw[ind,2])
        if self.display_drift_curves:
            plt.figure('Drift curves')
            figManager=plt.get_current_fig_manager()
            figManager.window.showMaximized()
            plt.plot(drift_raw[:,2],drift_raw[:,0]/self.magnification_factor,'xr',label='y (raw)')
            plt.plot(drift_raw[:,2],drift_raw[:,1]/self.magnification_factor,'xg',label='x (raw)')
            plt.plot(drift_interp[:,2],drift_interp[:,0]/self.magnification_factor,'--r',label='y (interpolated)')
            plt.plot(drift_interp[:,2],drift_interp[:,1]/self.magnification_factor,'--g',label='x (interpolated)')
            plt.xlabel('Time or frame')
            plt.ylabel('Drift value (correction pixels)')
            plt.legend()
        if self.display_correlation:
            for k in np.arange(self.nb_slices):
                plt.figure('Cross-correlation - Slice '+str(k)+' ( frame '+str(k*self.slice_size_t)+' )')
                plt.plot(drift_raw[k,1]+x0,drift_raw[k,0]+y0,'x')
            
        # Apply drift correction
        print('')
        print('Correcting drift...')
        drift_interp[:,0]/=self.magnification_factor
        drift_interp[:,1]/=self.magnification_factor
            
        self.data=self.format_data(self.data_init)
        for k in np.arange(np.shape(self.data_init)[0]):
            self.print_progress(k,np.shape(self.data_init)[0],50000,2)
            ind=np.argmax((drift_interp[:,2]-self.data[k,2]>0))-1
            self.data[k,[0,1]]+=drift_interp[ind,[0,1]]
        self.data_init['y']=self.data[:,0]
        self.data_init['x']=self.data[:,1]
            
        # Results export
        if self.export_results:
            np.save(self.filepath[:-4]+' - xydrift.npy',drift_interp)
            np.save(self.filepath[:-4]+'_xycorrected.npy',self.data_init)
            
            
        print('')
        print('Done')
        print('')
            
            
        
DCC_workflow=DCC_workflow()
DCC_workflow.Main()
