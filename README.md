# Evb-SMLM

Note: this is still work in progress currently. In particular, all the codes may not be available, or some documentation may be missing. You may want to check this repository regularly to have the last updates. If you are having trouble or would like to share feedback, do not hesitate to contact us.

--- General remarks

This set of codes is provided as part of the event-based Single Molecule Localization Microscopy (SMLM) project developed at Institut Langevin - ESPCI/CNRS, Paris, France by the team of Ignacio Izeddin and Clément Cabriel. See the published article [Cabriel, C., Monfort, T., Specht, C.G. et al. Event-based vision sensor for fast and dense single-molecule localization microscopy. Nat. Photon. 17, 1105–1113 (2023), https://doi.org/10.1038/s41566-023-01308-8] or the preprint [https://doi.org/10.1101/2022.07.22.501162] for more information.

Please use the following contact inforamtion if you would like to comment or contribute to our work, or if you have problems running the codes or questions about the technique or the data and processing steps:

* Clément Cabriel: clement.cabriel@espci.fr, cabriel.clement@gmail.com
* Ignacio Izeddin: ignacio.izeddin@espci.fr
    
If you use our code or data in your work and wish to publish, please make sure that you cite our work. If you would like to use them for commercial purposes, please contact us.
<!-- These codes and datasets is shared under the following licence [ADD LICENCE], all further work using them have to comply with the licencing terms. -->

--- Content of the repository and general comments

This repository contains some Python codes (.py) and some datasets containing either raw (.raw) or filtered or processed (.npy) datasets corresponding to event-based SMLM acquisitions. Note: some datasets are compressed sue to file size limitations. They must be decompressed before they can be used.

Acquisitions were performed with an event-based sensor manufactured by Prophesee (EVK V2 Gen4.1, Prophesee). The pixel size is 67nm, no PSF shaping is used. The acquisitions are performed with a dSTORM buffer under adequate illumination conditions to achieve blinking of the molecules (see the publication for more information).

All codes were developed and tested on a computer with Windows 10, Python version 3.8.8 and equipped with 128 GB of memory. Multithread options are available for some codes and were successfully tested. No GPU is required. The environment file "Eve_SMLM_environment.yml" provided in the repository can be used to run the codes.

Datasets containing initial data may be provided either in a .raw or in a .npy format. They are organized as arrays with four fields: 'x' (pixels), 'y' (pixels), 't' (us), 'p' (no unit). x and y are the lateral positions, t is the timestamp and p is the polarity (positive or negative).

Important: Some of the codes use the raw data reader provided by Metavision. Therefore, they require installing the Metavision Essentials package available here: [ https://www.prophesee.ai/metavision-intelligence-essentials-download/ ] or [ https://files.prophesee.ai/share/dists/public/windows/baiTh5si/ ]. We used Metavision Essentials 2.3.0 or 4.3.0. Note that these require specific versions of Python (3.8 or 3.9 for Metavision 4.3.0), and that the codes were tested only with Python 3.8.8. Other versions of Python will likely not work. It is important to follow closely the installation guidelines [ https://docs.prophesee.ai/stable/installation/windows.html ] until the "Installation" section included and to run the executable file with admin rights. Note that the functions provided by Metavision are only required to read the .raw data files, so .npy data files can be read without installing Metavision Essentials (if this installation step is skipped, please make sure that the line "from metavision_core.event_io.raw_reader import RawReader" is removed from the code).

See the section "Known issues" for information about fixes and workarounds if you encounter errors when running the codes.

--- Specific remarks about each dataset

* 'AF647_coverslip.raw': raw dataset containing 10 seconds of an acquisition performed on Alexa Fluor 647 deposited on a coverslip. Can be used with the read data code or with the localization code. Pixel size: 67 nm.
* 'AF647_tubulin.zip': raw dataset cotaining a part of an acquisition performed on Alexa Fluor 647 targetting the alpha-tubulin in fixed COS-7 cells. Can be used with the read data code or with the localization code. Pixel size: 67 nm.
* 'AF647_coverslip_sensitivity.zip': raw datasets containing 5 seconds of acquisitions performed on Alexa Fluor 647 deposited on a coverslip  with dSTORM buffer for different sensitivities. Acquisitions performed with the same conditions aside from the sensitivity. Contains 'AF647_coverslip-low_sentitivity', 'AF647_coverslip-standard_sensitivity', 'AF647_coverslip-high_sentitivity'. Can be used with the read data code or with the localization code. Pixel size: 67 nm.
* 'Beads40nm_coverslip_sensitivity.zip': raw datasets containing 8 seconds of acquisitions performed on 40-nm red beads deposited on a coverslip. Acquisitions were performed with the same conditions using a 10-Hz square modulated excitation. Contains 'Beads40nm_coverslip_low_sensitivity', 'Beads40nm_coverslip_standard_sensitivity', 'Beads40nm_coverslip_high_sensitivity'. Can be used with the read data code or with the localization code. Pixel size: 67 nm.
* 'Tubulin_DCC_test.zip': localized dataset containing the localizations of only the positive events from an acquisition of alpha-tubulin labelled with DNA-PAINT Atto655 imagers. The localization file is cropped and some fields are removed due to file size limitations. Intended to be used as a test sample for the 'Drift_correction.py' code. Pixel size: 80 nm.

--- Specific remarks about each code

* 'Read_data_evb.py': code to read raw datasets and perform some data cropping and optional basic data filtering and export the results as .npy datasets. Can be used to generate time frames and export them as .tif stack files.
* 'Localization_eventbased_blinking.py': code to perform the localization of molecules in raw blinking datasets
* 'Drift_correction.py': code to correct the drift after the localization of the molecules

--- Code user guide for 'Read_data_evb.py'

* Usage: download the repository and unzip it, unzip the datasets, set the parameters and run the code
* Recommended test datasets: 'AF647_coverslip_sensitivity.zip','AF647_coverslip.raw','AF647_tubulin.zip'
* Input: Events .raw file or events .npy file: (columns) 'x': x position (in pix), 'y': y position (in pix), 'p': event polarity, 't': timestamp (in us)
* Outputs:
    - Events .npy file: (columns) 'x': x position (in pix), 'y': y position (in pix), 'p': event polarity/algebric intensity, 't': timestamp (in us)
    - Frames .tif image stack
* Description of the parameters (more details are available directly in the code file):
    - 'filepath': complete path to the data file (e.g. 'C:/path/to/data.raw'). The data file can be .raw or .npy
    - 'buffer_size': memory allocated to the file import. Large files require higher values, wuich require the computer to have more memory. 1e9 works fine on a 128 GB memory computer, 1e7 should work on any computer and be sufficient for the small datasets provided
    - 'time_limits': time limits of the data pre-filtering
    - 'xy_limits': ROI for the data prefiltering
    - 'discard_up', 'discard_down': can be used to restrict the data to negative-only or positive-only events
    - 'skip_filtering': use to skip the spatio-temporal filtering step. This step is usually long and unnecessary for localization, but can be useful for data display. In the article, all the data used was unfiltered
    - 'space_gating', 'time_gating', 'detection_threshold': parameters of the spation-temporal filtering step. For each event, all the neighbors in a [-space_gating,+space_gating] pixels and [-time_gating,+time_gating] milliseconds neighborhood. The event is kept only if the number of neighbors is above or equal to detectionthreshold
    - 'export_event_data': used to export the filtered event list
    - 'export_frames': used to export the generated time frame sequence
    - 'xy_bin': xy binning to downsample the frames. Set =1 to skip the binning
    - 'time_bin': time bin of the frame generation
    - 'integrate_sighal': use to accumulate the events over time in order to reconstruct a camera-like image. This functionality is provided only for visualization convenience, and is not meant to be precise in terms of photometry
    - 'sign_display': event display type for the frames. Can be used to reconstruct the frames of the positive-only or negative-only events, or both regardless of their sign, or both respective of their sign
    - 'sum_all_frames': used to sum all the events of the filtered dataset in only one frame, typically to reconstruct a wide-field-like diffraction-limited image of the sample

--- Code user guide for 'Localization_eventbased_blinking.py':

* Usage: download the repository and unzip it, unzip the datasets, set the parameters and run the code
* Recommended test dataset: 'AF647_tubulin.zip'
* Input: Events .raw file or events .npy file: (columns) 'x': x position (in pix), 'y': y position (in pix), 'p': event polarity, 't': timestamp (in us)
* Outputs:
    - Events .npy file: (columns) ['id','t','y','x','N_up','N_down','t_up','t_down'] with all times in us and all distances in object plane pixels (more information in the code file)
    - Log .txt file
* Description of the parameters (more details are available directly in the code file):
    - 'filepath': complete path to the data file (e.g. 'C:/path/to/data.raw'). The data file can be .raw or .npy
    - 'buffer_size': memory allocated to import the data (in arbitrary units). Increase the buffer size if the input file is too large to be loaded. Decrease the buffer size if the PC is out of memory. Tested successfully with buffer size = 4e9 for 128 Gb of memory
    - 'pixel_size': pixel size in the object plane (in nm). =67 for the datasets provided
    - 'multithread': =True to use multithread parallelization; =False not to. Note: multithread parallelization increases calculation speed
    - 'time_limits': Time pre-filtering of the data. Only events in the range will be kept. =[t_min,t_max] (in seconds)
    - 'xy_limits': ROI selection (in pixels) [x_min,y_min,x_max,y_max]. =None not to filter. ='auto' to automatically detect the limits
    - 'threshold_detection': Wavelet detection threshold. Decreasing the value makes the detection more sensitive (and potentially more noisy)
    - 'exclusion_radius': Radius of the exclusion area (if two or more PSFs are closer than twice the value, they will all be discarded) (in pixels)
    - 'min_diameter': Minimum radius of the thresholded area (in pixels)
    - 'max_diameter': Maximum radius of the thresholded area (in pixels)
    - 'time_bin_frames': Time bin (in ms) of the frames. Used only for the PSF detection (the localization is performed on the eevent list and not on the frames)
    - 'area_radius': Radius of the fitting (or center of mass calculation) area (in pixels)
    - 'time_area_limits': Relative time limits (in ms) of the events to consider within an ROI. =[t_minus,t_plus]. The events considered for the localization are those within [t0-t_minus,t0+t_plus] where t0 is the detected time of the ROI
    - 'localization_method': Method for the PSF localization. ='COM' for center of mass calculation, ='Gaussian' for Gaussian fitting
    - 'display_frames': =True to display the frames, =False not to. Frames should be displayed to check the PSF detection parameters in particular, but the display should be disabled before running the complete localization code once the parameters are set.
    - 'number_frames_display': Number of frames to display
    - 'export_localized_results': =True to save the localized results; =False not to

--- Code user guide for 'Drift_correction.py':

* Usage: download the repository and unzip it, unzip the datasets, set the parameters and run the code. 
* Recommended test dataset: 'Tubulin_DCC_test.zip'
* Input: Localization .npy file. Required fields: 'y','x','t' (x,y in pixels, t in us). Note: t should start at 0. The localization file can be generated usiong the 'Localization_eventbased_blinking.py' code
* Outputs:
    - Drift-corrected .npy file: same format as input, same units as input
    - Drift .npy file: (columns) y,x,t with y,x in input pixels and t us
* Description of the parameters (more details are available directly in the code file):
    - 'filepath': complete path to the localized file (e.g. 'C:/path/to/data.npy'). The data file must be .npy
    - 'initial_pixel': Localization pixel size (in nm). Set =1 if the coordinates are already in nm. Value for the test dataset: =80
    - 'smoothing_frame': Sigma (in zoomed correction_pixel coordinates) of the Gaussian filtering of the frames. Default: =1
    - 'smoothing_correlation': Sigma (in zoomed correction_pixel coordinates) of the Gaussian filtering of the correlation functions (recommanded with maxima detection). Default: =1
    - 'frame_ref_t': Reference slice t limits (in us). =[t_min,t_max]. Suggested value: =[166e6,333e6]
    - 'slice_size_t': Slice width in t (in us). Suggested value: =150e6
    - 'max_shift': Max value of the allowed shift relative to the reference image (in nm). Default value: =1000
    - 'ROI': ROI limits (in initial localization pixels). =[x_min,y_min,x_max,y_max]. Default value: =[0,0,1e9,1e9] (no cropping)
    - 'maximum_detection_method': maximum_detection_method for the peak detection. ='Gaussian' or ='Max' (maximum detection, i.e. taking the highest pixel). Gaussian yields more precise results, but the fitting might fail in some cases. Default value: ='Gaussian'
    - 'width_ROI_gaussian': Half window of the fitting around the maximum value (in nm). Default value: =100.
    - 'interpolation_step': Interpolation time step (in us) for the raw drift curve interpolation. Default value: 1e7
    - 'export_results': set =True to export the results

--- Known issues:

Several errors returned by Python when running the code have been documented:

* "No module named 'metavision_core' ": this happens when Python cannot locate the installation folder of the Metavision SDK. Follow carefully the instructions given in the installation manual [ https://docs.prophesee.ai/stable/installation/windows.html ] until the "Installation" section included. In particular, make sure to execute the installation file with admin rights, follow the requirements with regards to the Python version, and make sure to update the PATH as described in the "Installing Dependencies" section. Alternatively, if you have .npy data files rather than .raw files, you can simply remove the line "from metavision_core.event_io.raw_reader import RawReader" from the code.
* "RawReader buffer size too small. Please increase max_events": this happens when the buffer allocated to import a .raw data file is too small. Increase the value of the variable "buffer_size" (a good value is usually between 1e8 and 4e9, depending on the available RAM in the computer).
* "Unable to allocate [A] GiB for an array with shape ([B],) and data type {'names':['x','y','p','t'], 'formats':['<u2','<u2','<i2','<i8'], 'offsets':[0,2,4,8], 'itemsize':16}" (where [A] and [B] may have various values): this happens when the buffer allocated to import a .raw data file is too large, with two possible reasons: either the requested buffer exceeds the RAM of the computer, or the memory was not released by Python after a previous code execution. To troubleshoot this problem, follow these steps: free as much RAM as possible by ending unnecessary processes and applications in the background, restart Python, and run the code. If the error persists, then it is due to the fact that the computer does not have enough RAM, so the solution is to decrease the value of the variable "buffer_size" (a good value is usually between 1e8 and 4e9, depending on the available RAM in the computer). If, on the contrary, the code can be run one or several times before the error reappears, then this is due to Python not releasing the memory after a run. In that case, it is required to reset the kernel (if running the code in the curent console or in a dedicated console) or close the active terminal(s) (if running the code in an external system terminal) before starting a new run. Note that it is currently not possible to load .raw files larger than approximately 25% of the RAM of the computer, however it is possible to slice a .raw file into several smaller files using the functions provided by Metavision [ https://docs.prophesee.ai/stable/samples/modules/driver/file_cutter.html?highlight=metavision%20raw%20cutter ] before loading them one by one in Python.
