# Up-and-Down-states-detection

* Scripts to detect Up and Down states during SWS

Mukovski:          https://academic.oup.com/cercor/article/17/2/400/31809

Saleem:            https://link.springer.com/article/10.1007%2Fs10827-010-0228-5


These scripts are based on two papers mentioned above, Mukovski and Saleem. Up and Down states are detected based on the phase of LFP at low frequency bands <10 Hz
in general. Here I explain in brief about all the scripts being used to detect Up and Down states.

<vikash_intra_UP_DOWN.m>

\\\Input argument/s:\\\             *.abf file containing Intracellualr recording on one channel and LFP on another channel. 


The purpose of this script is to extract useful bandwidths and theta values specific to the brain region of recording or anesthesia. This scripts starts with the processing of
intracellular recording.First thing it does to the intracellular data is to median filter it with a user defined median filter window. The second step is to
use a bandpass filter with the bandwidth [0.1 20], 0.1 to remove the drift in the data. Then it plots the intracellular recording which is first figure
and holds it for some time. After that it does the histogram plot of this intracellular recording (upper 1 percentile excluded) with user defined bin size. In the next step,
it tries to fit histograms with gaussian mixture distribution. After this gaussian mixture distribution, it calculates mean and standard deviation and by using that it calculates
the upper and lower threshold for detecting UP and DOWN states. With this threshold it calculates the periods for UP and DOWN states. Then corrects for the
boundaries and discards those which are lesser than a pre defined period for a state (UP or DOWN). Then it plots these UP and DOWN states on the figure 1 which 
was on hold. 

In the second part of the script, it operates on LFP. It starts with passing the frequency below 500 Hz or 200 Hz, to remove the MUA activity. Then it decimates/downsamples the
LFP to 1 KHz. Next, in figure 3, it plots the LFP and then plots the UP and DOWN states obtained in the previous step from intracelluar recordings. In this figure, you can
