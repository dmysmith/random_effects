function [fh] = showvol_randomeffects(dirname_out, fstem_imaging, dataRelease, RandomEffects)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEMO USING showVol.m WITH FEMA OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is designed to generate voxelwise random effects results from
% FEMA.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN showVol.m
%
% CLONE GitHub directories for cmig_utils, showSurf
% ADD ALL ABCD CMIG tools directories to MATLAB path:

% e.g. if cloned into ~/github: 
% addpath(genpath('~/github/cmig_utils'))
% addpath(genpath('~/github/showVol'))

%% Inputs to showvol_randomeffects:
% dirname_out    : path where FEMA results were saved
% fstem_imaging  : imaging phenotype e.g. 'area_ic5_sm1000'
% dataRelease    : ABCD data release (e.g.: '4.0')
% RandomEffects  : list of random effects e.g. {'F' 'A' 'T' 'S' 'E'}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZATION OF VOXELWISE FEMA RANDOM EFFECTS USING 'showVol'


% 1) Specify where FEMA_wrapper output is saved and load into MATLAB workspace
fname_results = sprintf('%s/FEMA_wrapper_output_voxel_%s.mat',dirname_out,fstem_imaging);
load(fname_results); % load FEMA output

% 2) Specify ABCD release version as the atlas used for voxelwise registration is different from 3.0 to 4.0
switch dataRelease
case '3.0'
  atlasVersion = 'ABCD1_cor10';
case '4.0'
  atlasVersion = 'ABCD2_cor10';
end
  
loadPrerenderedIfNeeded(atlasVersion) % loads atlas specific data
global PRI % saves atlas data as a global variable so does not need to be loaded every time you run showVol

%specify output to plot e.g. vol_z, vol_beta_hat etc
for i=1:length(RandomEffects)
    vol_stat{i} = fullvol(sig2mat(i,:),mask);
end

nvox=156662;
thresh_bonfcorr=abs(norminv(0.05/nvox));

% 3) Use a helper function to convert FEMA output into volumes for showVol. 

%  for each IV, creates a raw grayscale map and colorized map overlaid on T1 
limits = [-1 1];  % Limits for colormap: (1) min (2) max (3) threshold (leave third value empty to have unthresholded maps)
index = [1];                    % Indicates which IVs to plot: colnames_model(index)
showPval = false;                   % If showPval==true, creates a map of log10 p-values
cmap = redblackblue_alpha;                      % Colormap
cmap(:,4) = 1;
bg_vol = [];                        % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                       % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];                  % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution

for i=1:length(RandomEffects)
    vol{i} = convertFEMAVols(vol_stat{i}, fstem_imaging, RandomEffects{i}, colnames_model, index, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
end

% 4) Visualize using showVol
switch atlasVersion
      case {'ABCD1', 'ABCD1_cor10'}
            %for backwards compatibility, PRI naming is different for ABCD1
             showVol(vol{:}, PRI.vol_T1,PRI.vol_T2,PRI.vol_aseg_rgb,PRI.vol_fiber_rgb,PRI.vol_CO,PRI.vol_FOD_1mm)
      
      case {'ABCD2', 'ABCD2_cor10'}
            coords = [-7 -13 -4]; % Opens showVol to a specific location e.g. R NAcc
%             showVol(vol_F, vol_A, vol_S, vol_E, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)
            showVol(vol{:}, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)
end
