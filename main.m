clear;
addpath(genpath("../FID-A-master"));


cfg = struct();
cfg.labels_path = 'data/skeleton/labels.nii';
cfg.metabs_list = ["NAA","NAAG","Cr","PCr","GPC","PCh", ...
        "Glu","Gln","GABA","GSH","Asp","Tau","PE", "Lac"];
cfg.metabs_path = 'data/metabolites/metab_df.csv';
phantom = Phantom(cfg);


sm = SignalModel();

seq_func = @sim_press;
seq_params = struct();
seq_params.n = 1024;
seq_params.sw = 2000;
seq_params.B0 = 3;
seq_params.lw = 1/(150E-3);
seq_params.sys = 0; 
TE   = 30e-3;  
seq_params.tau1 = TE/4;             
seq_params.tau2 = TE/4;
voxrange = [24, 1, 1];
signal = sm.extract_signal(phantom,@sim_press,seq_params, voxrange);