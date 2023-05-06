load('/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03/designMat02_t1w_AgeSexScanSoft/FASE/FEMA_wrapper_output_vertex_sulc_ic5_sm1000.mat');
obs = table(iid,eid);
write.table(obs, '/space/syn50/1/data/ABCD/d9smith/random_effects/results_2023-03-03/designMat02_t1w_AgeSexScanSoft/FASE/obs.txt');   