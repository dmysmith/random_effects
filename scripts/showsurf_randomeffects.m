function [fh] = showsurf_randomeffects(dirname_out, fstem_imaging, dataRelease, ico, RandomEffects, savepath, varargin)

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

inputs = inputParser;
addParamValue(inputs, 'rgb', []);
addParamValue(inputs, 'legendPosition', 0);
addParamValue(inputs, 'title', 0)
addParamValue(inputs, 'polarity', 2)
addParamValue(inputs, 'colorbar', []);

parse(inputs, varargin{:});
legendPosition = inputs.Results.legendPosition;
title = inputs.Results.title;
rgb = inputs.Results.rgb;
polarity = inputs.Results.polarity;

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
fmid = 0.5; % middle of colorbar
fvals = [fmin fmid fmax];
clim = [fmin fmax]; % set colorbar limits

% cm = hot; % set colormap
% cm = blueblackred(); % set colormap
% cm = cm(51:end,:);

cm = fire;

curvcontrast = [0.2 0.2]; % contrast of gyri/sulci

bgcol = [0 0 0]; % change to [1 1 1] for white background

fh = figure('Units', 'centimeters', 'Position', [10 10 16 10], 'Color', bgcol, 'InvertHardcopy', 'off');

% Define spacing for axes
% hvgap controls the horizontal and vertical spaces between the axes
% btgap controls the space from the bottom of the figure and the top
% of the figure respectively
% lrgap controls the space from the left of the figure and the right
% of the figure respectively
hvgap = [0.02 0.02];
if strcmpi(legendPosition, 'south')
  lrgap = [0.02 0.02];
  if title
      btgap = [0.12 0.08];
  else
      btgap = [0.12 0.01];
  end
else
  if strcmpi(legendPosition, 'east')
      lrgap = [0.02 0.138];
      if title
          btgap = [0.018 0.08];
      else
          btgap = [0.018 0.018];
      end
  end
end

if ~isempty(rgb)
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

    btgap = [0.02 0.02];
    
    allH = tight_subplot(2, 2, hvgap, btgap, lrgap);
    hold(allH(:), 'on');
      
    axes(allH(1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    axes(allH(2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
    axes(allH(4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,rgbvals_lh,rgbvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
 
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



          allH = tight_subplot(1, 4, hvgap, btgap, lrgap);
          hold(allH(:), 'on');
    
          axes(allH(1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          axes(allH(2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          % axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          axes(allH(4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},polarity,curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
          
          if(title)
              titleAx = axes;
              set(titleAx,'position',[0 0 1 1],'units','normalized');axis off;
              text(titleAx, 0.5,1,sprintf('%s ~ %s [%s]', fstem_imaging, colnames_model{coeffnum}, statname),'color','w','fontweight','bold','interpreter','none','verticalalignment','top','horizontalalignment','center','fontsize',14)
          end
          
          if ~isempty(colorbar)
              % Set colorbar
    
              colormap(cm);
              cb                    = colorbar('color', 'w');
              cb.FontSize           = 10;
              % cb.Label.String       = strcat('z-score');
              cb.Label.FontSize     = 12;
              cb.Label.FontWeight   = 'bold';   
              cb.Box                = 'off';
        
              if strcmpi(legendPosition, 'south')
                  cb.Location = 'south';
                  if title
                      cb.Position(1)      = allH(1).Position(1);
                      cb.Position(2)      = cb.Position(2) - hvgap(1);
                      cb.Position(3)      = allH(4).Position(3)*4 + hvgap(1);
                  else
                      cb.Position(1)      = allH(1).Position(1);
                      cb.Position(2)      = cb.Position(2) - btgap(1);
                      cb.Position(3)      = 1- allH(1).Position(1) - hvgap(1);
                  end
              else
                  if strcmpi(legendPosition, 'east')

                      cb.Location = 'westoutside';
                      if title
                          cb.Position(1)      = allH(4).Position(1) + allH(4).Position(3) + 0.01;
                          cb.Position(2)      = allH(3).Position(2);
                          cb.Position(4)      = allH(1).Position(4)*2 + hvgap(1);
                      else                          
                          
                          cb.Position(1)      = allH(4).Position(1) + 0.16;
                          cb.Position(2)      = allH(4).Position(2) + hvgap(1);
                          cb.Position(4)      = allH(1).Position(4)*4;% + hvgap(1);
                          
                      end
                  end
              end
          end

          caxis(clim);

          % save
          if ~exist(savepath,'dir'), mkdir(savepath); end
          saveas(fh, sprintf('%s/%s_%s_%s.png',savepath, fstem_imaging, strrep(dirname_out(68:end),'/','_'),RandomEffects{coeffnum}));
    end
end