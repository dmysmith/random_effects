%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO USING showSurf.m WITH FEMA OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is designed to provide an end-to-end demo of how to
% visualize VERTEXWISE results generated from running FEMA

% `showSurf.m` is a tool to plot surface based statistics onto a pial or
% inflated cortical surface based on FreeSurfer segmentation.  This tool
% produces MATLAB figures across different views of the brain (lateral and
% medial; left and right hemispheres) that can be used for publications.

% For publication high resolution quality images we recommend running
% vertexwise analyses at icosahedral order (ico) 7. This will take longer to run, so we
% recommend running analyses first at ico 4 or 5, and then re-running all
% analyses at ico 7 when ready to make figures.

% Code written by Anders Dale, John Iversen, Clare E Palmer and Diliana Pecheva, 2021
%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN showSurf.m
%
% CLONE GitHub directories for cmig_utils, showSurf
% ADD ALL ABCD CMIG tools directories to MATLAB path:

% e.g. if cloned into ~/github: 
% addpath(genpath('~/github/cmig_utils'))
% addpath(genpath('~/github/showSurf'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISULATION OF VERTEXWISE FEMA OUTPUT USING 'showSurf'

% Inputs to showsurf_ds 
dirname_out = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-01-30/designMat1_dmri_AgeSexScanSoft/FATSE';
fstem_imaging = 'area_ic5_sm1000'; % thickness_ic5_sm1000 area_ic5_sm1000 sulc_ic5_sm1000
dataRelease = '4.0'; % ABCD data release
ico = 5; % ico number
ncoeff = 1:5;
RandomEffects = {'F' 'A' 'T' 'S' 'E'};
savepath = '/home/d9smith/tmp/2023-02-01';


showsurf_ds(dirname_out, fstem_imaging, dataRelease, ico, RandomEffects, savepath);
