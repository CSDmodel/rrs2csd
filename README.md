## rrs2csd.m Documentation

This function calclulates CSD slopes from remote sensing reflectance by using the CSD model. Remote sensing reflectance at 8 MODIS bands (412, 443, 469, 488, 531, 547, 555,and 678 nm) are required for this function. A large CSD slope indicates a large proportion of smaller phytoplankton, whereas a small CSD slope suggests that larger phytoplankton dominate. To avoid unrealistic values of the CSD slope, we defined the upper and lower limits of the CSD slope as 3.0 and 0.0, respectively.

Remote sensing refrectance can be downloaded from Ocean Color Website (https://oceancolor.gsfc.nasa.gov/).

For a full description of the CSD model, see following papers:

- Waga, H., Hirawake, T., Fujiwara, A., Kikuchi, T., Nishino, S., Suzuki, K., et al. (2017). Differences in Rate and Direction of Shifts between Phytoplankton Size Structure and Sea Surface Temperature. Remote Sens. 9, 222. doi:10.3390/rs9030222.

- Waga, H., Hirawake, T., and Ueno, H. (in revision). Spatial and temporal variations in impacts of mesoscale eddies on phytoplankton size structure. Geophys. Res. Lett.
