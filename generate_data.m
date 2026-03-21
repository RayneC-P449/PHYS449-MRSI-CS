clear;
addpath(genpath("../FID-A-master"));
cfg = struct(); 
cfg.metabs = readtable('data/metabolites/metab_df.csv');
cfg.metabs_list = ["Asp","NAAG","Cr","PCr","GPC","PCh", ...
    "Glu","Gln","GABA","GSH","NAA","Tau","PE", "Lac", "Ins"];

cfg.labels = niftiread('data/skeleton/labels.nii');
cfg.labels = cfg.labels(41:140,59:158,41:140);
pd_data = niftiread('data\skeleton\pd.nii');
label_data = cfg.labels;
pd_data = pd_data(41:140,59:158,41:140);
mm_json = jsondecode(fileread('data/macromolecules/mm_json.json'));
cfg.mm_list = mm_json.names;
cfg.mm = mm_json;
cfg.water = zeros(2,2);
cfg.water(1,:) = [36000, 52.6E-3];
cfg.water(2,:) = [43300, 63.4E-3];
cfg.pd = pd_data;
phantom = Phantom(cfg);
seqparams = struct;
seq_params.sw = 2000;
seq_params.B0 = 3;
seq_params.lw = 2;
TE = 30E-3;
seq_params.tau1 = TE/4;
seq_params.tau2 = TE/4;
seq_params.n = 1024;
load spinSystems
for l=1:numel(phantom.metabs_list)
   metab = phantom.metabs_list(l);
   seq_params.sys = eval(['sys' char(metab)]);
   basis = sim_press(seq_params.n, seq_params.sw, seq_params.B0, seq_params.lw, seq_params.sys, seq_params.tau1, seq_params.tau2);
   basis.fids = basis.fids;
   metab_bases{l} = basis;
end
orientations = {[2,3,1],[2,3,1],[2,3,1]};
nd_labels = ndims(phantom.labels);
nd_metab  = ndims(phantom.metab_data);
nd_mm     = ndims(phantom.mm_data);
nd_water  = ndims(phantom.water_data);
sm = SignalModel;
for o = 1:3
    orient = orientations{o};
    phantom.labels      = permute(phantom.labels, [orient, 4:nd_labels]);
    phantom.metab_data  = permute(phantom.metab_data, [orient, 4:nd_metab]);
    phantom.mm_data     = permute(phantom.mm_data, [orient, 4:nd_mm]);
    phantom.water_data  = permute(phantom.water_data, [orient, 4:nd_water]);
    pd_data = permute(pd_data, orient);
    label_data = permute(label_data, orient);
    slice_thickness = 10;
    for start_idx = 1:slice_thickness:size(phantom.labels,3)-slice_thickness+1
        end_idx = start_idx + slice_thickness - 1;
        path = fullfile(pwd, 'simulated', num2str(o), num2str(start_idx));
        [~,~] = mkdir(path);
        label_img = label_data(:,:,round((start_idx+end_idx)/2));
        vals = unique(label_img);
        n = numel(vals);
        % Map values to 1..n
        [~, img_idx] = ismember(label_img, vals);
        fig = figure('Visible','off');
        imagesc(img_idx);
        axis image off
        colormap(lines(n));
        caxis([0.5 n+0.5]);
        cb = colorbar;
        cb.Ticks = 1:n;
        cb.TickLabels = vals;   
        label_img_path = fullfile(path, 'labels.png');
        saveas(fig,label_img_path);
        close(fig);
        pd_img = pd_data(:,:,round((start_idx+end_idx)/2));
        fig = figure('Visible', 'off');
        pd_img = im2gray(pd_img);
        imshow(pd_img);
        pd_img_path = fullfile(path, 'MRI.png');
        saveas(fig,pd_img_path);
        close(fig);
        voxrange = [12,12]; sigma = 0; beta = 0;
        slice = start_idx:end_idx;
        data = struct;
        [data.kspace, data.t, data.ppm] = sm.extract_kspace(metab_bases, seq_params, phantom, slice, voxrange, sigma, beta); 
        k_path = fullfile(path, 'kspace.mat'); 
        data.MRI = pd_img;
        data.labels = label_img;
        save(k_path, 'data');
        disp(k_path);
    end
end
folder = fullfile(pwd,'simulated', 'basis');
[~,~] = mkdir(folder);
for m = 1:numel(phantom.metabs_list)
    path = fullfile(folder, sprintf('%s.RAW',phantom.metabs_list{m}));
    io_writelcmraw(metab_bases{m}, path, phantom.metabs_list(m));
end