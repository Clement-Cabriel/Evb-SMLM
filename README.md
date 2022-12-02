# Evb-SMLM

Note: this is still work in progress currently. In particular, all the codes may not be available, or some documentation might be missing. You may want to check this repository regularly to have the last updates. If you are having trouble, do not hesitate to contact us.

--- General remarks

This set of codes is provided as part of the event-based Single Molecule Localization Microscopy (SMLM) project developed at Institut Langevin - ESPCI/CNRS, Paris, France by the team of Ignacio Izeddin and Clément Cabriel. See the preprint [bioRxiv 2022.07.22.501162; https://doi.org/10.1101/2022.07.22.501162] for more information.
Please use the following contact inforamtion if you would like to comment or contribute to our work, or if you have problems running the codes or questions about the technique or the data and processing steps:
    Clément Cabriel: clement.cabriel@espci.fr, cabriel.clement@gmail.com
    Ignacio Izeddin: ignacio.izeddin@espci.fr
If you use our code or data in your work and wish to publish, please make sure that you cite our work. If you would like to use them for commercial purposes, please contact us.
<!-- These codes and datasets is shared under the following licence [ADD LICENCE], all further work using them have to comply with the licencing terms. -->

--- Content of the repository and general comments

This repository contains some Python codes (.py) and some datasets containing either raw (.raw) or filtered or processed (.npy) datasets corresponding to event-based SMLM acquisitions. Note: some datasets are compressed sue to file size limitations. They must be decompressed before they can be used.
Acquisitions were performed with an event-based sensor manufactured by Prophesee (EVK V2 Gen4.1, Prophesee). The pixel size is 67nm, no PSF shaping is used. The acquisitions are performed with a dSTORM buffer under adequate illumination conditions to achieve blinking of the molecules (see the publication for more information).
All codes were developed and tested on a computer with Windows 10, Python version 3.8.8 and equipped with 128 GB of memory. Multithread options are available for some codes and were successfully tested.
Files containing initial data may be provided either in a .raw or in a .npy format. They are organized as arrays with four fields: 'x' (pixels), 'y' (pixels), 't' (us), 'p' (no unit). x and y are the lateral positions, t is the timestamp and p is the polarity (positive or negative).

--- Specific remarks about each dataset

* 'AF647_coverslip.raw': raw dataset containing 10 seconds of an acquisition performed on Alexa Fluor 647 deposited on a coverslip. Can be used with the read data code or with the localization code.
* 'AF647_tubulin.zip': raw dataset cotaining a part of an acquiition performed on Alexa Fluor 647 targetting the alpha-tubulin in fixed COS-7 cells. Can be used with the read data code or with the localization code.
* 'AF647_coverslip_sensitivity.zip': raw datasets containing 5 seconds of acquisitions performed on Alexa Fluor 647 deposited on a coverslip  with dSTORM buffer for different sensitivities. Acquisitions performed with the same conditions aside from the sensitivity. Contains 'AF647_coverslip-low_sentitivity', 'AF647_coverslip-standard_sensitivity', 'AF647_coverslip-high_sentitivity'. Can be used with the read data code or with the localization code.
* 'Beads40nm_coverslip_sensitivity.zip': raw datasets containing 8 seconds of acquisitions performed on 40-nm red beads deposited on a coverslip. Acquisitions were performed with the same conditions using a 10-Hz square modulated excitation. Contains 'Beads40nm_coverslip_low_sensitivity', 'Beads40nm_coverslip_standard_sensitivity', 'Beads40nm_coverslip_high_sensitivity'. Can be used with the read data code or with the localization code.

--- Specific remarks about each code

* 'Read_data_evb.py': code to read raw datasets and perform some data cropping and optional basic data filtering and export the results as .npy datasets. Can be used to generate time frames and export them as .tif stack files.
* 'Localization_eventbased_blinking.py': [coming soon] code to perform the localization of molecules on blinking movies

--- Code user guide for 'Read_data_evb.py'

* Important: This code uses the raw data reader provided by Metavision. Therefore, it requires installing the Metavision Essentials package available here: https://www.prophesee.ai/metavision-intelligence-essentials-download/
  We used Metavision Essentials 2.3.0, which runs only with Python 3.7 or higher. Newer versions of Metavision Essentials may work.
* Usage: download the repository and unzip it, unzip the datasets, set the parameters and run the code
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

* Important: This code may use the raw data reader provided by Metavision (for .raw files only, not necessary for .npy files). Therefore, it requires installing the Metavision Essentials package available here: https://www.prophesee.ai/metavision-intelligence-essentials-download/
  We used Metavision Essentials 2.3.0, which runs only with Python 3.7 or higher. Newer versions of Metavision Essentials may work.
* Usage: download the repository and unzip it, unzip the datasets, set the parameters and run the code
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
    - 'display_frames': =True to display the frames, =False not to. Frames should be displayed to check the PSF detection parameters in particular.
    - 'number_frames_display': Number of frames to display
    - 'export_localized_results': =True to save the localized results; =False not to

