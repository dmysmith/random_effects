function [fh] = showsurf_randomeffects(dirname_out, fstem_imaging, dataRelease, ico, RandomEffects, savepath, rgb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% REQUIREMENTS TO RUN showSurf.m
%
% CLONE GitHub directories for cmig_utils, showSurf
% ADD ALL ABCD CMIG tools directories to MATLAB path:

% e.g. if cloned into ~/github: 
% addpath(genpath('~/github/cmig_utils'))
% addpath(genpath('~/github/showSurf'))

%% Inputs to showsurf_ds:
% dirname_out    : path where FEMA results were saved
% fstem_imaging  : imaging phenotype e.g. 'area_ic5_sm1000'
% dataRelease    : ABCD data release (e.g.: '4.0')
% ico            : ico number
% RandomEffects  : list of random effects e.g. {'F' 'A' 'T' 'S' 'E'}
% savepath       : path to folder to save figures

%% Optional arguments:
% rgb            : if provided, list of random effects to use for RGB map.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISULATION OF VERTEXWISE FEMA RANDOM EFFECTS USING 'showSurf'

fname_results = sprintf('%s/FEMA_wrapper_output_vertex_%s.mat',dirname_out,fstem_imaging); % FEMA output filename
load(fname_results); % load FEMA output

% In order to plot statistics on a cortical surface, we require a template
% of the cortical surface to project the results onto.  Below we load these
% surface templates.  `SurfView_surfs.mat` contains surfaces for all icos.

load SurfView_surfs.mat % load surface templates
% ico=5; % ico number
icnum=ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order
  
% The below code will produce figures for the random effects (sig2mat)

ncoeff=1:length(RandomEffects);

% specify limits for plot based on vertvals
fmax = 1; % max limit for colorbar
fmin = 0; % min limit of colorbar
fmid = fmax/2; % middle of colorbar
fvals = [fmin fmid fmax];
clim = [fmin fmax]; % set colorbar limits

cm = hot; % set colormap

curvcontrast = [0.2 0.2]; % contrast of gyri/sulci

bgcol = [0 0 0]; % change to [1 1 1] for white background

if exist('rgb')
    assert(length(rgb)==3,'Argument rgb must have length 3.');
    
    redi = find(ismember(RandomEffects, rgb{1}(:)));
    grni = find(ismember(RandomEffects, rgb{2}(:)));
    blui = find(ismember(RandomEffects, rgb{3}(:)));
    
    redvals = sum(sig2mat(redi,:),1);
    grnvals = sum(sig2mat(grni,:),1);
    bluvals = sum(sig2mat(blui,:),1);
    
    rgbvals = cat(1, redvals, grnvals, bluvals);
    rgbvals_lh = rgbvals(:,1:icnvert,:); % divide statistics by hemisphere for plotting
    rgbvals_rh = rgbvals(:,icnvert+[1:icnvert],:);
    
    fh = figure; clf; % number matlab figure window
    set(fh,'Color',bgcol); fh.InvertHardcopy = 'off';
    
    % subplot(2,2,1); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;

    subplot(2,2,1); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    subplot(2,2,2); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    subplot(2,2,3); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    subplot(2,2,4); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    titleAx = axes;
    set(titleAx,'position',[0 0 1 1],'units','normalized');axis off;
    %text(titleAx, 0.5,1,sprintf('%s ~ %s [%s]', fstem_imaging, colnames_model{coeffnum}, statname),'color','w','fontweight','bold','interpreter','none','verticalalignment','top','horizontalalignment','center','fontsize',14)

    % TODO: Figure out color bar
    
    % save
    if ~exist(savepath,'dir'), mkdir(savepath); end
    saveas(fh, sprintf('%s/%s_%s_%s.png',savepath, fstem_imaging, strrep(dirname_out(68:end),'/','_'),strjoin(rgb,'_')));
    
else

    for coeffnum = ncoeff

          statname = 'sig2mat'; % specify name of statistic plotting for figure label
          vertvals = sig2mat(coeffnum,:); % specify statistics to plot
          vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
          vertvals_rh = vertvals(icnvert+[1:icnvert]);



          fh = figure(coeffnum + 100*(str2num(dataRelease(1))-3)); clf; % number matlab figure window
          set(fh,'Color',bgcol); fh.InvertHardcopy = 'off';
          subplot(2,2,1); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          subplot(2,2,2); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          subplot(2,2,3); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          subplot(2,2,4); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
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

          % save
          if ~exist(savepath,'dir'), mkdir(savepath); end
          saveas(fh, sprintf('%s/%s_%s_%s.png',savepath, fstem_imaging, strrep(dirname_out(68:end),'/','_'),RandomEffects{coeffnum}));
    end
end