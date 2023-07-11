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
addpath(genpath('~/github/cmig_tools_internal'))
% addpath(genpath('~/github/showSurf'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISULATION OF VERTEXWISE FEMA OUTPUT USING 'showSurf'

% Inputs to showsurf_ds 
dirname_out = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03/designMat02_t1w_AgeSexScanSoft/FASE';
% dirname_out = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-02-17/designMat02_t1w_AgeSexScanSoft/FATSE';
modality = 'smri'; % dmri or smri
dataRelease = '4.0'; % ABCD data release
ico = 3; % ico number
RandomEffects = {'F' 'A' 'S' 'E'};
rgb = {'A' 'F' 'S'};
savepath = '/home/d9smith/projects/random_effects/plots/results_2023-03-03/FASE/horiz';

% Specify visual preferences for plotting
legendPosition = 'South'; % 'South' or 'East'
title = 1; % whether to include title at top of plot
polarity = 1;
curvcontrast = [0.2 0.2]; % contrast of gyri/sulci
bgcol = [0 0 0];

if modality == 'smri'
    pheno_list = {'area_ic5_sm1000' 'thickness_ic5_sm1000' 'sulc_ic5_sm1000'};
    
elseif modality == 'dmri'
    pheno_list = {'FA-gm' 'FA-wm' 'FNI-gm' 'FNI-wm' 'HNT-gm' 'HNT-wm' 'LD-gm' 'LD-wm' 'MD-gm' 'MD-wm' 'RND-gm' 'RND-wm' 'RNI-gm' 'RNI-wm' 'RNT-gm' 'RNT-wm' 'TD-gm' 'TD-wm'}; % name of imaging phenotype - data already saved as ico=5
    pheno_list = strcat(pheno_list, '_ic5_sm1000'); % add '_ic5_sm1000'
end

load SurfView_surfs.mat % load surface templates
icnum=ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order

for i=1:length(pheno_list)
    fstem_imaging = pheno_list{i};

    fname_results = sprintf('%s/FEMA_wrapper_output_vertex_%s.mat',dirname_out,fstem_imaging); % FEMA output filename
    load(fname_results); % load FEMA output

    % plot total residual variance
    statname = 'total residual variance';
    vertvals = sig2tvec; % specify statistics to plot
    vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
    vertvals_rh = vertvals(icnvert+[1:icnvert]);
    
    % specify limits for plot based on vertvals
    fmax = min(300,max(abs(vertvals))); % max limit for plotting purposes
    fmin = 0.0; % min limit for plotting purposes
    fmid = fmax/2; % middle value for plotting purposes
    fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
    
    % set colorbar limits - usually [fmin fmax] or [-fmax fmax]
    clim = [fmin fmax]; 
    
    % Create figure
    fh = figure('Units', 'centimeters', 'Position', [10 10 16 5], 'Color', bgcol, 'InvertHardcopy', 'off');
    
    % Define spacing for axes
    hvgap = [0.02 0.02];
    if strcmpi(legendPosition, 'south')
        lrgap = [0.02 0.02];
        btgap = [0.2 0.01];
    else
        if strcmpi(legendPosition, 'east')
          % lrgap = [0.02 0.138];
          lrgap = [0.02 0.17];
          btgap = [0.018 0.018];
        end
    end
    
    % Create axes
    allH = tight_subplot(1,4, hvgap, btgap, lrgap);
    hold(allH(:), 'on');

    cm = fire;
    
    axes(allH(1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    axes(allH(2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    axes(allH(4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    
    
    % Set colorbar
    colormap(cm);
    cb                    = colorbar('color', 'w');
    cb.FontSize           = 10;
    cb.Label.Interpreter  = 'latex';
    cb.Label.String       = '$\sigma^2_{total}$';
    cb.Label.FontSize     = 12;
    cb.Label.FontWeight   = 'bold';   
    cb.Box                = 'off';
    
    if strcmpi(legendPosition, 'south')
        cb.Location = 'south';
        cb.Position(1)      = allH(1).Position(1);
        cb.Position(2)      = cb.Position(2) - btgap(1);
        cb.Position(3)      = 1- allH(1).Position(1) - hvgap(1);
    else
        if strcmpi(legendPosition, 'east')
          cb.Location = 'eastoutside';
          cb.Position(1)      = allH(4).Position(1) + allH(4).Position(3) + 0.16;
          cb.Position(2)      = allH(3).Position(2) + .035;
          % cb.Position(4)      = allH(1).Position(4)*2 + hvgap(1);
          cb.Position(4)      = allH(1).Position(4) - .05;
        end
    end
    caxis(clim);

    % save
    if ~exist(savepath,'dir'), mkdir(savepath); end
    saveas(fh, sprintf('%s/%s_%s_%s.png',savepath, fstem_imaging, strrep(dirname_out(68:end),'/','_'),'sig2tvec'));

    % create plots for individual random effects
    fh = showsurf_randomeffects(dirname_out, fstem_imaging, dataRelease, ico, RandomEffects, savepath, 'legendPosition', legendPosition, 'polarity', polarity);

    % put vertices into bins for ROI-wise analysis
    lookup = generateLookup_vertices(sig2mat, [], [], [1 2 3 4], 'none', []);

    for randeff=1:length(lookup)
        writetable(lookup{randeff}, sprintf('%s/%s_table_%s.csv', dirname_out, pheno_list{i}, RandomEffects{randeff}), 'Delimiter', ',', 'QuoteStrings', 'all');
    end

end

% create plots for RGB values
for i=1:length(pheno_list)
    fstem_imaging = pheno_list{i};
    showsurf_randomeffects(dirname_out, fstem_imaging, dataRelease, ico, RandomEffects, savepath, 'rgb', rgb, 'legendPosition', legendPosition, 'polarity', polarity);
end

