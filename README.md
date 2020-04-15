# SpectralFitofHRF
Estimation of the haemodynamic response function from resting-state fMRI data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% README: Matlab code for Spectral fit of HRF 			          %
%                                                                         %
% Author: Carolin Lennartz                                                %
% Email: carolin.lennartz@uniklinik-freiburg.de                           %
% Date: April 15th, 2020                                              	  %
%                                                                         %
% CITATION: Lennartz et. al, Estimation of the hemodynamic response func- %
% tion from resting state fMRI spectra, submitted to Human Brain Mapping	%
%									                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Given a data set, the code estimates the HRF from the time series using a HRF dictionary (HRF_dictionary)

For estimation use 'fitHRFtosprectrumwithDict.m' with 'Dict_big.HRF' from HRF_dictionary. 'balloon.m' contains 
code to model Balloon model HRFs.
'correctHRFlatency.m' allows to correct the time series with both time to peak shift and with Wiener filter,
where 'deconv_wiener.m' contains the function to apply the Wiener filter.
With 'fitHRFtospectrumwithDict_4D.m' a complete 4D dataset can be used. The dataset needs to be brain extracted
and all voxels outside of the brain need to be set to zero.

For fMRI data with higher TR than 0.1s, the same Dictionary can be used. The TR only needs to be indicated.
