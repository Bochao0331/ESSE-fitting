"ESSE-fitting: echo-shifted spin-echo complex non-li" <br />
Example Raw data locate in: **/mnt/sdata_new/Bochao/0528AppleESSE/RAW/ on MREL server**; <br />
Example MAT files locate in: **/mnt/sdata_new/Bochao/0528AppleESSE//MAT/ on MREL server** <br />
Two spaital-resolution (1x1 mm^2 and 1.5x1.5 mm^2) 2D Echo-shift Single-echo Spin Echo data are included.<br />

This script does reading, and seperat exponential fitting of 2D TSE data with different TE or GRE timing.<br />
To fit R2* (R2* = R2+R2'), need the zero-shift (zero_shift_ind) to the last postive echo-shift<br />
To fit R2up (R2up = R2-R2'), need the first negative echo-shift ot the zero-shift
R2 and R2' can be calculated by the above equations. <br />
**Need to be more improved by considering off-resoance and non-linear inversion**
