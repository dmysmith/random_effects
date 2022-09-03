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

dirname_out='/users/dsmith/OneDrive - UC San Diego/PhD/ABCD Projects/random_effects/imaging/results/results_20220825/designMat1_allcovs/FASE'; % directory of where FEMA output saved
fstem_imaging='sulc_ic5_sm1000'; % thickness_ic5_sm1000 area_ic5_sm1000 sulc_ic5_sm1000

outpath = 'Downloads';

dataRelease='4.0'; % ABCD data release

fname_results = sprintf('%s/FEMA_wrapper_output_vertex_%s.mat',dirname_out,fstem_imaging); % FEMA output filename
load(fname_results); % load FEMA output

% In order to plot statistics on a cortical surface, we require a template
% of the cortical surface to project the results onto.  Below we load these
% surface templates.  `SurfView_surfs.mat` contains surfaces for all icos.

load SurfView_surfs.mat % load surface templates
ico=5; % ico number
icnum=ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order
  
% The demo below will produce figures for the IVs (columns of X) from
% the FEMA analysis specfiied by `ncoeff`

%% DIANA TODO: make matrix of vertvals so that you can plot them all together!

statname = 'sig2mat'; % specify name of statistic plotting for figure label
vertvals = sig2mat; % specify statistics to plot
vertvals_lh = vertvals(:,1:icnvert); % divide statistics by hemisphere for plotting
vertvals_rh = vertvals(:,icnvert+[1:icnvert]);

% specify limits for plot based on vertvals
fmax = 1; % max limit for colorbar
fmin = 0; % min limit of colorbar
fmid = fmax/2; % middle of colorbar
fvals = [fmin fmid fmax];
clim = [fmin fmax]; % set colorbar limits

cm = parula; % set colormap
% cm = redblackblue;

curvcontrast = [0.2 0.2]; % contrast of gyri/sulci

bgcol = [0 0 0]; % change to [1 1 1] for white background

fh = figure(101); clf; % number matlab figure window
set(fh,'Color',bgcol); fh.InvertHardcopy = 'off';

ncoeff=size(sig2mat,1);


for coeffi=1:ncoeff
      rowi = coeffi-1;
      subplot(ncoeff,ncoeff,(rowi*ncoeff+1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh(coeffi,:),vertvals_rh(coeffi,:),fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      subplot(ncoeff,ncoeff,(rowi*ncoeff+2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh(coeffi,:),vertvals_rh(coeffi,:),fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      subplot(ncoeff,ncoeff,(rowi*ncoeff+3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh(coeffi,:),vertvals_rh(coeffi,:),fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
      subplot(ncoeff,ncoeff,(rowi*ncoeff+4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh(coeffi,:),vertvals_rh(coeffi,:),fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;

end

titleAx = axes;
set(titleAx,'position',[0 0 1 1],'units','normalized');axis off;
%text(titleAx, 0.5,1,sprintf('%s ~ %s [%s]', fstem_imaging, colnames_model{coeffnum}, statname),'color','w','fontweight','bold','interpreter','none','verticalalignment','top','horizontalalignment','center','fontsize',14)

% Set colorbar
colormap(cm);
cb = colorbar('color', 'w');
% cb.Label.Interpreter = 'latex';
cb.Label.String = strcat('sig2mat');
cb.Label.FontSize = 10;
cb.Box = 'off';
cb.Position = [.92 .08 .02 .8150];
caxis(clim);
saveas(fh, strcat(outpath,"/",fstem_imaging,".png"))
