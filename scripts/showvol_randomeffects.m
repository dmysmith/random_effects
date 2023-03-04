function [fh] = showvol_randomeffects(dirname_out, fstem_imaging, dataRelease, RandomEffects, rgb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN showVol.m WITH FEMA RANDOM EFFECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is designed to call showVol.m with the necessary arguments to
% display multiple random effects.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs to showvol_randomeffects:
% dirname_out    : path where FEMA results were saved
% fstem_imaging  : imaging phenotype e.g. 'area_ic5_sm1000'
% dataRelease    : ABCD data release (e.g.: '4.0')
% RandomEffects  : list of random effects e.g. {'F' 'A' 'T' 'S' 'E'}
% rgb            : if provided, list of random effects to use for RGB map.

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

nvox=156662;
thresh_bonfcorr=abs(norminv(0.05/nvox));



%  for each IV, creates a raw grayscale map and colorized map overlaid on T1 
limits = [-1 1];  % Limits for colormap: (1) min (2) max (3) threshold (leave third value empty to have unthresholded maps)
index = [];                    % Indicates which IVs to plot: colnames_model(index)
showPval = false;                   % If showPval==true, creates a map of log10 p-values
bg_vol = [];                        % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                       % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];                  % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution

% Colormap
cmap = redblackblue_alpha;
cmap(:,4) = 1;

if exist('rgb')==1
    assert(length(rgb)==3,'Argument rgb must have length 3.');

    redi = find(ismember(RandomEffects, rgb{1}(:)));
    grni = find(ismember(RandomEffects, rgb{2}(:)));
    blui = find(ismember(RandomEffects, rgb{3}(:)));
    
    redvals = fullvol(sum(sig2mat(redi,:),1), mask);
    grnvals = fullvol(sum(sig2mat(grni,:),1), mask);
    bluvals = fullvol(sum(sig2mat(blui,:),1), mask);

    rgbvals = cat(4, redvals, grnvals, bluvals);

    % Use a helper function to convert FEMA output into volumes for showVol.
    rgbvols = convertFEMAVols(rgbvals, fstem_imaging, strcat(rgb{:}), colnames_model, index, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, rgb);

    switch atlasVersion
        case {'ABCD1', 'ABCD1_cor10'}
            % for backwards compatibility, PRI naming is different for ABCD1
            showVol(rgbvols, PRI.vol_T1,PRI.vol_T2,PRI.vol_aseg_rgb,PRI.vol_fiber_rgb,PRI.vol_CO,PRI.vol_FOD_1mm)
          
        case {'ABCD2', 'ABCD2_cor10'}
            coords = [-7 -13 -4]; % Opens showVol to a specific location e.g. R NAcc
            % showVol(vol_F, vol_A, vol_S, vol_E, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)
            % showVol(rgbvols, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)
            showVol(rgbvols, PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)

    end



else

    for i=1:length(RandomEffects)
        % specify output to plot e.g. vol_z, vol_beta_hat etc
        vol_stat{i} = fullvol(sig2mat(i,:),mask);
        % Use a helper function to convert FEMA output into volumes for showVol. 
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

end