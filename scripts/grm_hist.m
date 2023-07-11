load('/space/syn65/1/data/abcd-sync/4.0/genomics/ABCD_rel4.0_grm.mat'); % This does not include fids
iid = iid_list;

% initialize family IDs from design matrix
tbl_birthid = readtable('/space/syn50/1/data/ABCD/d9smith/random_effects/designMat/designMat02_t1w_AgeSexScanSoft.txt');
[dummy IA IB] = intersect(iid,tbl_birthid.src_subject_id,'stable');

fid = NaN(size(iid));
fid(IA) = tbl_birthid.rel_family_id(IB);

GRM_tmp = GRM(IA,IA);
fid_tmp = tbl_birthid.rel_family_id(IB);
FRM = fid_tmp==fid_tmp'; % Family Relatedness Matrix

indvec = find(tril(FRM<1));
grmvec = GRM_tmp(indvec);
[ivec jvec] = ind2sub(size(FRM),indvec);

[sv si] = sort(GRM_tmp(indvec),'descend'); si = si(isfinite(sv)); sv = sv(isfinite(sv));

for i = rowvec(si(sv>0.2))
  fprintf(1,'id1=%5d id2=%5d GRM=%0.2f\n',ivec(i),jvec(i),grmvec(i));
end

trilmat = tril(true(size(GRM_tmp)),-1);
trilind = find(trilmat);

if 0 % old code
    figure; hist(GRM_tmp(FRM(trilind)<1),100);
    
    figure; hist(GRM(FRM(trilind)>0.5),100);
    figure; hist(GRM_tmp(FRM(trilind)>0.5),100);
end

maskmat_fam = (FRM>0.5).*trilmat;
maskmat_notfam = (FRM<0.5).*trilmat;

[hc hv] = hist(GRM_tmp(find(maskmat_fam)),100); chc = cumsum(hc); % Inside family
figure(1); 
subplot(2,1,1); plot(hv,hc,'LineWidth',2); h=title('Histogram of pairwise GRM within families'); set(h,'FontSize',14,'FontWeight','bold'); axis tight;
subplot(2,1,2); semilogy(hv,sum(hc)-cumsum(hc),'LineWidth',2); h=[title('Cumulative histogram of pairwise GRM within families') xlabel('GRM') ylabel('Count')]; set(h,'FontSize',14,'FontWeight','bold'); axis tight;

[hc hv] = hist(GRM_tmp(find(maskmat_notfam)),100); chc = cumsum(hc); % outside family
figure(2); 
subplot(2,1,1); plot(hv,hc,'LineWidth',2); h=[title('Histogram of pairwise GRM between families') xlabel('GRM') ylabel('Count')]; set(h,'FontSize',14,'FontWeight','bold'); axis tight;
subplot(2,1,2); semilogy(hv,sum(hc)-cumsum(hc),'LineWidth',2); h=[title('Cumulative histogram of pairwise GRM between families') xlabel('GRM') ylabel('Count')]; set(h,'FontSize',14,'FontWeight','bold'); axis tight;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures and Info for paper

% Default figure settings
width = 6;     % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% Histogram of relatedness in full sample
figure(3);clf;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); % Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); % Set properties
histogram(GRM_tmp(trilind), 100);
set(gca,'YScale', 'log');
xlabel('Genetic Relatedness');
title('Genetic Relatedness (Full Sample)');

if input('\nSave figure? (1 or 0) ')
    exportgraphics(gcf, '/home/d9smith/projects/random_effects/plots/hist_fullsample.png', 'Resolution',300);
end

% Histogram of relatedness within families
figure(4);clf;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); % Set size
set(gca, 'FontSize', fsz, 'LineWidth', alw); % Set properties
histogram(GRM_tmp(find(maskmat_fam)), 100);
% set(gca,'YScale', 'log');
xlabel('Genetic Relatedness');
title('Genetic Relatedness (Related Sample)');

if input('\nSave figure? (1 or 0) ')
    exportgraphics(gcf, '/home/d9smith/projects/random_effects/plots/hist_relsample.png', 'Resolution',300);
end

% Number of families at baseline and y2
% baseline
numfams_bl = length(unique(tbl_birthid(strcmp(tbl_birthid.eventname,'baseline_year_1_arm_1'),:).rel_family_id));
fprintf('%i unique families at baseline\n',numfams_bl);

% y2
numfams_y2 = length(unique(tbl_birthid(strcmp(tbl_birthid.eventname,'2_year_follow_up_y_arm_1'),:).rel_family_id));
fprintf('%i unique families at Year 2\n',numfams_y2);

% number of related pairs (total)
GRM_unrelpairs = sum(GRM_tmp(trilind)<0.25,'all');
GRM_sibpairs = sum((GRM_tmp(trilind)>0.25&GRM_tmp(trilind)<0.75),'all');
GRM_twinpairs = sum(GRM_tmp(trilind)>0.75,'all');

fprintf('%d total pairs with relatedness < 0.25 across timepoints\n',GRM_unrelpairs);
fprintf('%d total pairs with relatedness 0.25 - 0.75 across timepoints\n',GRM_sibpairs);
fprintf('%d total pairs with relatedness > 0.75 across timepoints\n',GRM_twinpairs);

% Number of within-family related pairs
fprintf('%i within-family pairs across timepoints\n',sum(maskmat_fam,'all'));
fprintf('%i within-family pairs with relatedness data across timepoints\n',sum(isfinite(GRM_tmp(find(maskmat_fam))),'all'));

unrelpairs = sum(GRM_tmp(find(maskmat_fam))<0.25,'all');
sibpairs = sum((GRM_tmp(find(maskmat_fam))>0.25&GRM_tmp(find(maskmat_fam))<0.75),'all');
twinpairs = sum(GRM_tmp(find(maskmat_fam))>0.75,'all');

fprintf('%d within-family pairs with relatedness < 0.25 across timepoints\n',unrelpairs);
fprintf('%d within-family pairs with relatedness 0.25 - 0.75 across timepoints\n',sibpairs);
fprintf('%d within-family pairs with relatedness > 0.75 across timepoints\n',twinpairs);


