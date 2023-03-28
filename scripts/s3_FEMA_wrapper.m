%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run FEMA for random effects imaging analysis 
%% Diana Smith
%% Created March 2022
%% Last Updated August 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify where to store results
outpath = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03';

if ~exist(outpath, 'dir')
      mkdir(outpath)
end

% start diary
diary_path=strcat(outpath,'/','diary_', datestr(now, 'yyyy-mm-dd_HHMM'));
diary(diary_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD ALL ABCD CMIG tools directories to MATLAB path:
addpath(genpath('/home/d9smith/github/cmig_tools_internal'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify data release

dataRelease = '4.0'; %'3.0' or '4.0'

% run abcdConfig
cfg = abcdConfig('FEMA');
abcd_sync_path=cfg.data.abcd_sync;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify which imaging analyses to demo
doVertexwiseSmri = 1; % run vertexwise smri analysis (datatype = 'vertex')
doVertexwiseDmri = 0; % run vertexwise dmri analysis
doVoxelwiseSmri = 0; % run voxelwise smri analysis (datatype = 'voxel')
doVoxelwiseDmri = 0; % run voxelwise dmri analysis (datatype = 'voxel')
doMOSTest = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUTS TO FEMA_wrapper.m

dirname_tabulated = fullfile(abcd_sync_path,dataRelease,'tabulated/released'); % directory to tabulated imaging data on abcd-sync 
atlasVersion = 'ABCD2_cor10';
dirname_tabulated = fullfile(abcd_sync_path,'4.0','tabulated/released'); %KNOWN ISSUE: breaks when using txt files following NDA release --> must use pre-release csv files
fname_pihat = fullfile('/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat'); 

%To run multiple deisgn matrices with same imaging data populate each row with path to each design matrix
designmat_dir = '/space/syn50/1/data/ABCD/d9smith/random_effects/designMat';
designmat_file = dir(sprintf('%s/designMat*.txt', designmat_dir));
% designmat_file = {designmat_file.name}';
% fname_design = strcat(designmat_dir, '/', designmat_file);

% only running one design mat at a time
fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/designMat/designMat02_t1w_AgeSexScanSoft.txt';
designmat_file = 'designMat02_t1w_AgeSexScanSoft.txt';

outdir_file = strrep(designmat_file, '.txt', '');
outdir_path=strcat(outpath,'/',outdir_file);

% for running just one designmat
% fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/designMat/designMat11_dmri_AgeSexScanSoftGest.txt';
% outdir_path = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-02-07/designMat11_dmri_AgeSexScanSoftGest.txt';

% Optional inputs for `FEMA_wrapper.m` depending on analysis
contrasts=[]; % Contrasts relate to columns in design matrix e.g. [1 -1] will take the difference between cols 1 and 2 in your design matrix (X)
ranknorm = 0; % Rank normalizes dependent variables (Y) (default = 0)
nperms = 0; % Number of permutations - if wanting to use resampling methods nperms>0
mediation = 0; % If wanting to use outputs for a mediation analysis set mediation=1 - ensures same resampling scheme used for each model in fname_design
PermType = 'wildbootstrap'; %Default resampling method is null wild-bootstrap - to run mediation analysis need to use non-null wild-bootstrap ('wildboostrap-nn')
tfce = 0; % If wanting to run threshold free cluster enhancement (TFCE) set tfce=1 (default = 0)
RandomEstType = 'ML'; % specify random effects estimator (default is MoM)
Hessflag=0;
logLikflag=0;
ciflag=0;

colsinterest=[1]; % Only used if nperms>0. Indicates which IVs (columns of X) the permuted null distribution and TFCE statistics will be saved for (default 1, i.e. column 1)
niter=0; % decrease number of iterations -- change when you want to run for real!

fname_pregnancyID = fullfile('/home/sabad/requests/pregnancy_ID_01172023.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION 1: F, A, T, S, E

% specify random effects - updating to accept a cell array
% RandomEffects = {{'F','A','T','S','E'}};
RandomEffects = {{'F','A','T','S','E'}; {'F','A','S','E'}};

for r=1:length(RandomEffects)

      % specify where to store results
      outdir_label = strcat(RandomEffects{r}{:});
      outDir = strcat(outdir_path, '/', outdir_label);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% VERTEXWISE ANALYSIS
      if doVertexwiseDmri
            datatype='vertex'; % imaging modality selected
            modality='dmri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)
  
            % Uses path structure in abcd-sync to automatically find data
            dirname_imaging = fullfile(abcd_sync_path, dataRelease, 'imaging_concat/vertexwise/', modality); % filepath to imaging data
            dirname_out = fullfile(outDir); % filepath to save FEMA output
  
            switch dataRelease
                  case '3.0'
                        %fstem_imaging = 'area-sm256';
                        fstems = 'thickness-sm256'; % name of imaging phenotype
                        %fstem_imaging = 'sulc-sm256';
                  case '4.0'
                        fstems ={'FA-gm' 'FA-wm' 'FNI-gm' 'FNI-wm' 'HNT-gm' 'HNT-wm' 'LD-gm' 'LD-wm' 'MD-gm' 'MD-wm' 'RND-gm' 'RND-wm' 'RNI-gm' 'RNI-wm' 'RNT-gm' 'RNT-wm' 'TD-gm' 'TD-wm'}; % name of imaging phenotype - data already saved as ico=5
                        fstems = strcat(fstems, '_ic5_sm1000'); % add '_ic5_sm1000'
            end
  
            ico = 5; % icosahedral number
  
            % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line
            for m=1:length(fstems)
                  fstem_imaging = fstems{m};
                  % RUN FEMA
                  % RUN FEMA
                  [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
                  'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects{r}, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'colsinterest',colsinterest);
            end
  
            %%
            if doMOSTest && nperms>0
                  % Compute the MOSTest statistics with the permutation scheme
                  if strcmp(modality,'smri')
                        % Set regularization factor
                        % has only been computed for smri thickness & surface area for now
                        if contains(fstem_imaging, 'thickness')
                              k=47;
                        elseif contains(fstem_imaging, 'area')
                              k=79;
                        else
                              warning('Regularization has not been optimized for this set of parameters. MOSTest results are probably off.')
                        end
                        alpha = 0.05;
                        col_interest = 1;
                        [vertex_MOSTest_perm_pval, vertex_MOSTest_extrap_pval, vertex_MOSTest_mostvec, ...
                              vertex_MOSTest_cthresh] = FEMA_MOSTest(zmat_perm, col_interest, alpha, k);
                        sprintf('MOSTest results: p-value=%.3f, extrapolated p-value=%.3f.', perm_pval, extrap_pval)
                  end
            end 
      end

      if doVertexwiseSmri
          
          
          datatype='vertex'; % imaging modality selected
          modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

          % Uses path structure in abcd-sync to automatically find data
          dirname_imaging = fullfile(abcd_sync_path, dataRelease, 'imaging_concat/vertexwise/', modality); % filepath to imaging data
          dirname_out = fullfile(outDir); % filepath to save FEMA output

          switch dataRelease
                case '3.0'
                      %fstem_imaging = 'area-sm256';
                      fstems = 'thickness-sm256'; % name of imaging phenotype
                      %fstem_imaging = 'sulc-sm256';
                case '4.0'
                      % fstems ={'area_ic5_sm1000' 'sulc_ic5_sm1000' 'thickness_ic5_sm1000'}; % name of imaging phenotype - data already saved as ico=5
                      fstems = {'thickness_ic5_sm1000' 'area_ic5_sm1000' 'sulc_ic5_sm1000'};

          end

          ico = 3; % icosahedral number

          % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line

          for m=1:length(fstems)
                fstem_imaging = fstems{m};
                
                % RUN FEMA
                [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
                'ico', ico, 'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects{r}, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'colsinterest',colsinterest,...
                'Hessflag',Hessflag,'ciflag',ciflag,'logLikflag',logLikflag,'RandomEstType',RandomEstType);
          end

          %%
          if doMOSTest && nperms>0
                % Compute the MOSTest statistics with the permutation scheme
                if strcmp(modality,'smri')
                      % Set regularization factor
                      % has only been computed for smri thickness & surface area for now
                      if contains(fstem_imaging, 'thickness')
                            k=47;
                      elseif contains(fstem_imaging, 'area')
                            k=79;
                      else
                            warning('Regularization has not been optimized for this set of parameters. MOSTest results are probably off.')
                      end
                      alpha = 0.05;
                      col_interest = 1;
                      [vertex_MOSTest_perm_pval, vertex_MOSTest_extrap_pval, vertex_MOSTest_mostvec, ...
                            vertex_MOSTest_cthresh] = FEMA_MOSTest(zmat_perm, col_interest, alpha, k);
                      sprintf('MOSTest results: p-value=%.3f, extrapolated p-value=%.3f.', perm_pval, extrap_pval)
                end
          end
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% VOXELWISE ANALYSES

      if doVoxelwiseDmri

          datatype = 'voxel'; % imaging modality selected
          modality='dmri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)

          % uses path structure in abcd-sync to automatically find data
          dirname_imaging = fullfile(abcd_sync_path, dataRelease, '/imaging_concat/voxelwise/', atlasVersion, modality); % filepath to imaging data
          dirname_out = fullfile(outDir,dataRelease); % filepath to save FEMA output

          modality = {'RNT' 'RNI' 'RND' 'RIF' 'RDF' 'HNT' 'HNI' 'HND' 'HIF' 'HDF' 'FNI' 'FA' 'MD' 'JA'};
          % modality = {'RNI' 'RND' 'FNI'}; % for running just one modality

          % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line
          for m=1:length(modality)
                fstem_imaging=modality{m};

                % Run FEMA
                [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
                'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects{r}, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'colsinterest',colsinterest);
          end    
      end

      if doVoxelwiseSmri

            datatype = 'voxel'; % imaging modality selected
            modality='smri'; % concatenated imaging data stored in directories based on modality (smri, dmri, tfmri, roi)
  
            % uses path structure in abcd-sync to automatically find data
            dirname_imaging = fullfile(abcd_sync_path, dataRelease, '/imaging_concat/voxelwise/', atlasVersion, modality); % filepath to imaging data
            dirname_out = fullfile(outDir,dataRelease); % filepath to save FEMA output
  
            % modality = {'RNT' 'RNI' 'RND' 'RIF' 'RDF' 'HNT' 'HNI' 'HND' 'HIF' 'HDF' 'FNI' 'FA' 'MD' 'JA'};
            modality = {'nu'}; % for running just one modality
  
            % Once all filepaths and inputs have been specified FEMA_wrapper.m can be run in one line
            for m=1:length(modality)
                  fstem_imaging=modality{m};
                  keyboard;
  
                  % Run FEMA
                  [fpaths_out beta_hat beta_se zmat logpmat sig2tvec sig2mat beta_hat_perm beta_se_perm zmat_perm sig2tvec_perm sig2mat_perm inputs mask tfce_perm analysis_params] = FEMA_wrapper(fstem_imaging, fname_design, dirname_out, dirname_tabulated, dirname_imaging, datatype,...
                  'ranknorm', ranknorm, 'contrasts', contrasts, 'RandomEffects', RandomEffects{r}, 'pihat_file', fname_pihat, 'nperms', nperms, 'mediation',mediation,'PermType',PermType,'tfce',tfce,'preg_file',fname_pregnancyID,'colsinterest',colsinterest);
            end    
        end
end

diary off
