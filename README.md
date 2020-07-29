# STM-C
A remote sensing driven permafrost carbon model

STM-C is a remote sensing driven permafrost carbon modelling framework, which mainly contains two modules: a 1-D soil thermal model and a terrestrial carbon flux model. This modelling framework is designed to simulate soil temperature and carbon decomposition dynamics in the northern permafrost region at intermediate scales. It is mainly driven by remote sensing data. STM-C is partly built on the original pan-Arctic permafrost hydrology model (PWBM). 

 1. Compiling the model in Linux
   gfortran -g *.f -o *.exe -mcmodel=medium
   
 2. Running the model in Linux
   *.exe *.init
   The initialization file provides the path and filename that were needed to run the model. An example was provided in ./test/
 
 3. Input ancillary files and datasets
   Two types of files are needed to run the model. Ancillay files including soil texture, peat fraction, annual mean temperature for the domain are needed. Model drivers 
   include land surface temperature, soil moisture, snow depth and density, downward solar radiation. An example of the input files will be uploaded to figshare. 
   
 4. References
 
   Yi, Y., Kimball, J. S., Chen, R. H., Moghaddam, M., Reichle, R. H., Mishra, U., Zona, D., and Oechel, W. C.: Characterizing permafrost active layer dynamics and sensitivity to landscape spatial heterogeneity in Alaska, The Cryosphere, 12, 145–161, https://doi.org/10.5194/tc-12-145-2018, 2018. 
   Yi, Y., Kimball, J. S., Watts, J. D., Natali, S. M., Zona, D., Liu, J., Ueyama, M., Kobayashi, H., Oechel, W., and Miller, C. E.: Investigating the sensitivity of soil respiration to recent snow cover changes in Alaska using a satellite-based permafrost carbon model, Biogeosciences Discuss., https://doi.org/10.5194/bg-2020-182, in review, 2020. 
   Rawlins, M. A., Nicolsky, D. J., McDonald, K. C. and Romanovsky, V. E.: Simulating soil freeze/thaw dynamics with an improved pan-Arctic water balance model: SOIL FREEZE/THAW MODELING, J. Adv. Model. Earth Syst., 5(4), 659–675, doi:10.1002/jame.20045, 2013.
 
