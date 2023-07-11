% visualize relatedness of analytic sample
% Diana Smith
% July 2023

% specify paths to design matrix, results, and GRM
fname_design = '/space/syn50/1/data/ABCD/d9smith/random_effects/designMat/designMat02_t1w_AgeSexScanSoft.txt';
% results_file = '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03/designMat02_t1w_AgeSexScanSoft/FASE/FEMA_wrapper_output_vertex_area_ic5_sm1000.mat';
GRM_file = '/space/amdale/1/tmp/ABCD_cache/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat';

% load results file (need this for list of participant IDs)
designmat = readtable(fname_design);
obs = designmat(:, 1:3);
obs.Properties.VariableNames = ["iid","eid","fid"];

% load GRM file
load(GRM_file);

% get list of IDs
rows = strcmp(obs.eid,'baseline_year_1_arm_1');
baseline_ids = obs(rows,:).iid;

% index GRM by IDs
idx = ismember(iid_list,baseline_ids);
GRM_baseline = GRM(idx,idx);

GRM_baseline_tri = triu(GRM_baseline, 1);
GRM_baseline_tri(GRM_baseline_tri==0)=NaN;

% Default figure settings
width = 6;     % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% Part 1: histogram of relatedness using all participants at baseline
figure(1);clf;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); % Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); % Set properties
histogram(GRM_baseline_tri, 300);
set(gca,'YScale', 'log');
xlabel('Genetic Relatedness');
title('Genetic Relatedness (Full Sample, Baseline Visit)');

exportgraphics(gcf, '/home/d9smith/projects/random_effects/plots/hist_fullsample.png', 'Resolution',300);

% Part 2: histogram of relatedness using only participants in same family
% at baseline

% initialize empty matrix
GRM_related = NaN(size(GRM));

% find duplicate fids within analytic_sample_ids
fid_dups = obs(rows,:).fid(find(hist(obs(rows,:).fid,length(unique(obs(rows,:).fid)))>1));

% from FEMA_parse_family
for fi = 1:length(fid_dups) % for each fid in fid_dups:

    fi_ids = obs(obs.fid==fid_dups(fi)&strcmp(obs.eid,'baseline_year_1_arm_1'),:).iid; % Identify all iids for a given family at baseline

    % get relatedness values for pts in this family
    idx = ismember(iid_list,fi_ids);
    fi_GRM = GRM(idx,idx);

    % Impute missing values within families (assume all proper siblings => phat = 0.5)
    if any((~isfinite(fi_GRM)),'all')
        [row,col] = find(~isfinite(fi_GRM));
        fi_GRM(row,col) = 0.5;
    end

    % populate GRM_related with values from fi_GRM
    GRM_related(idx,idx) = fi_GRM;

end

GRM_related_tri = triu(GRM_related, 1);
GRM_related_tri(GRM_related_tri==0)=NaN;

% histogram of upper triangle of new matrix
figure(2);clf;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); % Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); % Set properties
histogram(GRM_related_tri, 300);
% set(gca,'YScale', 'log');
xlabel('Genetic Relatedness');
title('Genetic Relatedness (Related Sample, Baseline Visit)');

exportgraphics(gcf, '/home/d9smith/projects/random_effects/plots/hist_relsample.png', 'Resolution',300);

% Part 3: Number of families and related pairs - for table
% baseline, pairs that share fid
length(unique(obs(strcmp(obs.eid,'baseline_year_1_arm_1'),:).fid)) % 9529 families
sum(GRM_related_tri(:)>0.25&GRM_related_tri(:)<0.75) % 685 pairs between 0.25-0.75
sum(GRM_related_tri(:)>0.75) % 155 pairs >0.75

% checking full baseline sample - why so much more?
sum(GRM_baseline_tri(:)>0.25&GRM_baseline_tri(:)<0.75) % 1387 pairs between 0.25-0.75
sum(GRM_baseline_tri(:)>0.75) % 375 pairs >0.75

% year 2 follow up
y2 = strcmp(obs.eid,'2_year_follow_up_y_arm_1');
length(unique(obs(y2,:).fid)) % 6524 families

y2_ids = obs(y2,:).iid;

% index GRM by IDs
idx = ismember(iid_list,y2_ids);
GRM_y2 = GRM(idx,idx);

GRM_y2_tri = triu(GRM_y2, 1);
GRM_y2_tri(GRM_y2_tri==0)=NaN;

% full y2 sample
sum(GRM_y2_tri(:)>0.25&GRM_y2_tri(:)<0.75) % 845 pairs between 0.25-0.75
sum(GRM_y2_tri(:)>0.75) % 262 pairs >0.75

%y2, pairs that share fid

% initialize empty matrix
GRM_related_y2 = NaN(size(GRM));

% find duplicate fids within analytic_sample_ids
fid_dups = obs(y2,:).fid(find(hist(obs(y2,:).fid,length(unique(obs(y2,:).fid)))>1));

% from FEMA_parse_family
for fi = 1:length(fid_dups) % for each fid in fid_dups:

    fi_ids = obs(obs.fid==fid_dups(fi)&strcmp(obs.eid,'2_year_follow_up_y_arm_1'),:).iid; % Identify all iids for a given family at y2

    % get relatedness values for pts in this family
    idx = ismember(iid_list,fi_ids);
    fi_GRM = GRM(idx,idx);

    % Impute missing values within families (assume all proper siblings => phat = 0.5)
    if any((~isfinite(fi_GRM)),'all')
        [row,col] = find(~isfinite(fi_GRM));
        fi_GRM(row,col) = 0.5;
    end

    % populate GRM_related with values from fi_GRM
    GRM_related_y2(idx,idx) = fi_GRM;

end

GRM_related_y2_tri = triu(GRM_related_y2, 1);
GRM_related_y2_tri(GRM_related_y2_tri==0)=NaN;

sum(GRM_related_y2_tri(:)>0.25&GRM_related_y2_tri(:)<0.75) % 452 pairs between 0.25-0.75
sum(GRM_related_y2_tri(:)>0.75) % 151 pairs >0.75