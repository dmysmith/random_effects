%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO USING showVol.m WITH FEMA OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is designed to provide an end-to-end demo of how to
% visualize VOXELWISE results generated from running FEMA.

% `showVol.m` is a tool to plot whole brain voxelwise statistics. It uses
% an interactive GUI that allows the user to move around the brain, explore
% effect from different planes and overlay regions of interest (ROIs).
% For the ABCD Study, subcortical parcellations from freesurfer and
% other external atlases have been registered to the ABCD atlas space. For
% these to provide accurate spatial labelling, statistical volumes also need
% to be within the ABCD atlas space. The function `convertFEMAVols` ensures
% volume data are aligned with these ROIs. This function also allows the user to
% set color overlays of statistical effects on top of a background image
% (default: T1w image).  The below demo shows the user how to
% visualize ACBD analysis results that used data registered the ABCD atlas space.

% Code written by Anders Dale, John Iversen, Clare E Palmer and Diliana Pecheva, 2021
%
% This software is Copyright (c) 2021 The Regents of the University of California. All Rights Reserved.
% See LICENSE.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO SCRIPT

% Run `FEMA_showVol_demo.m` in MATLAB as a script for a full end-to-end demo.
% The demo will prompt you to give user specific inputs in the command window.

% If you would like to edit the script manually to change some settings then SAVE A LOCAL
% COPY of this script before running.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN showVol.m
%
% CLONE GitHub directories for cmig_utils, showSurf
% ADD ALL ABCD CMIG tools directories to MATLAB path:

% e.g. if cloned into ~/github: 
addpath(genpath('/usr/pubsw/packages/matlab/R2022a/toolbox/shared/io/general'));
addpath(genpath('~/github/cmig_tools_internal'))
% addpath(genpath('/home/d9smith/github/cmig_tools_internal/showVol/utils/'))
% addpath(genpath('~/github/showVol'))

%git_dir=input('Please specify the path to your github directories: ','s');

%fprintf('\nAdding github directories (FEMA, cmig_utils, showSurf and showVol) to MATLAB path... \n')
%addpath(genpath(sprintf('%s/cmig_utils',git_dir)))
%addpath(genpath(sprintf('%s/showVol',git_dir)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISULATION OF VOXELWISE FEMA OUTPUT USING 'showVol'

% testing rgb argument
rgb = {'A' 'FT' 'S'};

% 1) Specify where FEMA_wrapper output is saved and load into MATLAB workspace
datestring = 'results_2023-02-17';
dirname_out=strcat('/space/syn50/1/data/ABCD/d9smith/random_effects','/',datestring); % directory of where FEMA output saved
designMat='designMat01_dmri_AgeSexScanSoft'; %  designMat02_t1w_AgeSexScanSoft designMat11_dmri_AgeSexScanSoftGest.txt
RandomEffects = {'F' 'A' 'T' 'S' 'E'};
fstem_imaging='RNI'; % MD FA RNI RND FNI RNT HNT JA (NOTE: use t1w for JA)

% 2) Specify ABCD release version as the atlas used for voxelwise registration is different from 3.0 to 4.0
dataRelease='4.0';

% 3) Specify path where you would like to save images
savepath = '/home/d9smith/projects/random_effects/plots/rgb';
savepath = sprintf('%s/%s/%s/%s/%s',savepath, datestring, designMat, strcat(RandomEffects{:}), fstem_imaging);

% call showvol_randomeffects
dirname_out = sprintf('%s/%s/%s/%s', dirname_out, designMat, strcat(RandomEffects{:}),dataRelease);
showvol_randomeffects(dirname_out, fstem_imaging, dataRelease, RandomEffects, rgb);

% create directory (if it doesn't exist) and navigate there
if ~exist(savepath,'dir'), mkdir(savepath); end
cd(savepath);

% While running, use these keyboard shortcuts to save images
% -- 'v': cycle to the next volume
% -- 'o': cycle orientations
% -- '!': save screenshot of main axis