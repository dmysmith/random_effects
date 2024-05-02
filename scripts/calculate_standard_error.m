% Calculate standard error of estimates for random effects paper
% Diana Smith
% May 2024

%% Calculate standard error based on 100 MoM permutations
% load results file
results_dir = '/space/ceph/1/ABCD/users/d9smith/random_effects/results_mom_100perms_2024-05-01/designMat02_t1w_AgeSexScanSoft/FASE/nonnullWB_100perms/';

% path to save plots
savedir = '/home/d9smith/projects/random_effects/plots/std_err';

RandomEffects = {'F' 'A' 'S' 'E'};
modality = {'thickness' 'area' 'sulc'};

legendPosition = 'south';
ico = 3;

% loop through results files
for m=1:length(modality)
    clear sig2mat_perm;
    results_file = sprintf('FEMA_wrapper_output_vertex_%s_ic5_sm1000.mat', modality{m});
    load(fullfile(results_dir,results_file));

    % The standard error is defined as the standard deviation / sqrt(n)
    std_err = std(sig2mat_perm(:,:,2:end),[],3)/sqrt(size(sig2mat_perm(:,:,2:end),3));  

    for re=1:length(RandomEffects)
        % save map of standard error
        vertvals = std_err(re,:);
        statname = sprintf('%s ~ %s (standard error)', modality{m}, RandomEffects{re});
        savefile = sprintf('%s_std_err_%s.png', modality{m}, RandomEffects{re});
        
        if ~exist(savedir, 'dir')
            mkdir(savedir)
        end

        % specify limits for plot based on vertvals
        fmax = 0.0222; % max limit for colorbar - max for thickness 0.0169; area .0222, sulc is .0208
        fmin = 0; % min limit of colorbar
        fmid = fmax/2; % middle of colorbar
        fvals = [fmin fmid fmax];
        clim = [fmin fmax]; % set colorbar limits
        
        FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
        'colormap', fire, 'legendPosition', legendPosition, 'ico', ico, 'savepath', fullfile(savedir, savefile));
    
        % save thresholded map based on 95% CI
        ci = prctile(sig2mat_perm(re,:,2:end),[2.5 97.5], 3);
        mask = ci(:,:,1)>0;

        vertvals = sig2mat(re,:) .* mask;
        statname = sprintf('%s ~ %s (thresholded)', modality{m}, RandomEffects{re});
        savefile = sprintf('%s_thresholded_%s.png', modality{m}, RandomEffects{re});

        % specify limits for plot based on vertvals
        fmax = 1; % max limit for colorbar - max for thickness is 0.0169
        fmin = 0; % min limit of colorbar
        fmid = 0.5; % middle of colorbar
        fvals = [fmin fmid fmax];
        clim = [fmin fmax]; % set colorbar limits

        FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
        'colormap', fire, 'legendPosition', legendPosition, 'ico', ico, 'savepath', fullfile(savedir, savefile));
    
    
    end

    

end