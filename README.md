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

This repository contains some Python codes (.py) and some datasets containing either raw (.raw) or filtered or processed (.npy) datasets corresponding to event-based SMLM acquisitions.
Acquisitions were performed with an event-based sensor manufactured by Prophesee (EVK V2 Gen4.1, Prophesee). The pixel size is 67nm, no PSF shaping is used. The acquisitions are performed with a dSTORM buffer under adequate illumination conditions to achieve blinking of the molecules (see the publication for more information).
All codes were developed and tested on a computer with Windows 10, Python version 3.8.8 and equipped with 128 GB of memory. Multithread options are available for some codes and were successfully tested.
Files containing initial data may be provided either in a .raw or in a .npy format. They are organized as arrays with four fields: 'x' (pixels), 'y' (pixels), 't' (us), 'p' (no unit). x and y are the lateral positions, t is the timestamp and p is the polarity (positive or negative).

--- Specific remarks about each dataset and code

* 'AF647_coverslip.raw': raw dataset containing 10 seconds of an acquisition performed on Alexa Fluor 647 deposited on a coverslip. Can be used with the read data code or with the localization code.
* 'AF647_tubulin.zip': [coming soon] raw dataset cotaining a part of an acquiition performed on Alexa Fluor 647 targetting the alpha-tubulin in fixed COS-7 cells. Can be used with the read data code or with the localization code.
* 'AF647_tubulin.zip': raw dataset cotaining a part of an acquiition performed on Alexa Fluor 647 targetting the alpha-tubulin in fixed COS-7 cells. Can be used with the read data code or with the localization code.
* 'AF647_coverslip_sensitivity.zip': raw datasets containing 5 seconds of acquisitions performed on Alexa Fluor 647 deposited on a coverslip  with dSTORM buffer for different sensitivities. Acquisitions performed with the same conditions aside from the sensitivity. Contains 'AF647_coverslip-low_sentitivity', 'AF647_coverslip-standard_sensitivity', 'AF647_coverslip-high_sentitivity'. Can be used with the read data code or with the localization code.
* 'Read_data_evb.py': code to read raw datasets and perform some data cropping and optional basic data filtering and export the results as .npy datasets. Can be used to generate time frames and export them as .tif stack files.
* 'Localization_eventbased_blinking.py': code to perform the localization of molecules on blinking movies

--- Code user guide for 'Read_data_evb.py'

* Important: This code uses the raw data reader provided by Metavision. Therefore, it requires installing the Metavision Essentials package available here: https://www.prophesee.ai/metavision-intelligence-essentials-download/
  We used Metavision Essentials 2.3.0, which runs only with Python 3.7 or higher. Newer versions of Metavision Essentials may work.
* Description of the parameters: [coming soon]

--- Code user guide for 'Localization_eventbased_blinking.py':
* [coming soon]
