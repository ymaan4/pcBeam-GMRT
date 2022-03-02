# pcBeam-GMRT
### A package to form post-correlation beam using the phased-array (PA/CDPA) and incoherent-array (IA) data from GMRT.

#### Usage Information

$ bin/pcBeam_gmrt -h

pcBeam_gmrt - Construct post-correlation beam using GMRT's IA and PA/CDPA beam data

usage: pcBeam_gmrt -{options} 

options:

 -f1 <file-name> : GMRT data file1 (PA/CDPA)
 -f2 <file-name> : GMRT data file2 (IA)
 
 -nch1 <nchan>   : No. of channels in file1 
 -ts1  <tsamp>   : Sampling time (in ms) in file1 
 -nch2 <nchan>   : No. of channels in file2 
 -ts2  <tsamp>   : Sampling time (in ms) in file2 
 -scale <fact>   : scale input-2 by 'fact' before subtraction 
-o filename - specify output filename (def=stdout)


#
Yogesh Maan  <ymaan4[@]gmail.com>
