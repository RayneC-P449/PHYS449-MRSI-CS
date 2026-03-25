clear;
addpath(genpath("../FID-A-master"));

function phantom = load_phantom()
    cfg = struct(); 
    cfg.metabs = readtable('data/metabolites/metab_df.csv');
    cfg.metabs_list = ["Asp","NAAG","Cr","PCr","GPC","PCh", ...
        "Glu","Gln","GABA","GSH","NAA","Tau","PE", "Lac", "Ins"];
    cfg.labels = niftiread('data/skeleton/labels.nii');
    slices = {50:130,78:138, 50:130};
    cfg.labels = cfg.labels(slices{:});
    cfg.pd = niftiread('data\skeleton\pd.nii');
    cfg.pd = cfg.pd(slices{:});
    mm_json = jsondecode(fileread('data/macromolecules/mm_json.json'));
    cfg.mm_list = mm_json.names;
    cfg.mm = mm_json;
    cfg.water = zeros(2,2);
    cfg.water(1,:) = [36000, 52.6E-3];
    cfg.water(2,:) = [43300, 63.4E-3];    
    phantom = Phantom(cfg);
end


function permuted_phantom = permute_phantom(phantom, permutation)
    permuted_phantom = phantom;
    permuted_phantom.labels = permute(phantom.labels, permutation);
    permuted_phantom.pd = permute(phantom.pd, permutation);
    permuted_phantom.metab_data  = permute(phantom.metab_data, [permutation, 4:ndims(phantom.metab_data)]);
    permuted_phantom.mm_data     = permute(phantom.mm_data, [permutation, 4:ndims(phantom.mm_data)]);
    permuted_phantom.water_data  = permute(phantom.water_data, [permutation, 4:ndims(phantom.water_data)]);
end

function simulate(setup_params_list, orientations, permutations)
    simulation_path = fullfile(pwd, sprintf('%s', 'simulations'));
    [~,~] = mkdir(simulation_path);
    for i = 1:numel(setup_params_list)
        setup_params = setup_params_list{i};
        setup_path = fullfile(simulation_path, sprintf('%s_%d', 'setup', setup_params.num));
        [~,~] = mkdir(setup_path);
        phantom = load_phantom();
        seq_params = setup_params.seq_params;
        mbasis = load_mbasis(phantom.metabs_list, seq_params, seq_params.sequence);
        basis_path = fullfile(setup_path, sprintf('%s', 'basis'));
        [~,~] = mkdir(basis_path);
        for m = 1:numel(mbasis)
            path = fullfile(basis_path, sprintf('%s.RAW', phantom.metabs_list{m}));
            io_writelcmraw(mbasis{m}, path, phantom.metabs_list(m));
        end
        lcm = LCMHelper();
        lcm.make_basis(basis_path);
        voxels = setup_params.voxels; snr = setup_params.snr; 
        beta = setup_params.beta; B0_map = setup_params.B0_map;
        slice_thickness = setup_params.slice_thickness;
        sm = SignalModel();
        for o = 1:numel(orientations)
            orientation = orientations{o};
            orientation_path = fullfile(setup_path, orientation);
            [~,~] = mkdir(orientation_path);
            permutation = permutations{o};
            permuted_phantom = permute_phantom(phantom, permutation);
            for start_idx = 1 : slice_thickness : size(permuted_phantom.labels,3) - slice_thickness + 1
                end_idx = start_idx + slice_thickness - 1;
                slice_path = fullfile(orientation_path, sprintf('%s_%s_%d-%d',orientation,'slice',start_idx,end_idx));
                [~,~] = mkdir(slice_path);
                slice = start_idx:end_idx;
                [kspace, t, ppm] = sm.extract_kspace(mbasis, seq_params, permuted_phantom, slice, voxels, snr, beta, B0_map);
                data = struct();
                data.K = kspace; data.t = t; data.ppm = ppm; data.X = fftshift(fftn(kspace));
                data.setup_params = setup_params;
                data_path = fullfile(slice_path, 'data.mat');
                pd_slice_idx = round((start_idx + end_idx) / 2);
                pd_slice = permuted_phantom.pd(:,:, pd_slice_idx);
                data.pd = pd_slice;
                save(data_path, 'data');
                fig = figure('Visible','off');
                pd_img = im2gray(pd_slice);
                imshow(pd_img);
                pd_img_path = fullfile(slice_path, 'MRI.png');
                saveas(fig, pd_img_path);
                close(fig);
                labels_slice_idx = round((start_idx + end_idx) / 2);
                labels_slice = permuted_phantom.labels(:,:, labels_slice_idx);
                unique_labels = unique(labels_slice);
                [~, labels_img] = ismember(labels_slice, unique_labels);
                fig = figure('Visible','off');
                imagesc(labels_img);
                axis image off
                n = numel(unique_labels);
                colormap(lines(n));
                clim([0.5, n+0.5]);
                cb = colorbar;
                cb.Ticks = 1:n;
                cb.TickLabels = unique_labels;   
                labels_img_path = fullfile(slice_path, 'labels.png');
                saveas(fig, labels_img_path);
                close(fig);
                fprintf('Simulated %s\n', slice_path);
                vis = Visualizer();
                mrsi_path = fullfile(slice_path, 'MRSI.png');
                vis.visualize(flip(data.ppm), data.X, 'reverse', [0, max(abs(data.X),[],'all')], data.pd, mrsi_path, false);
                fig = figure('Visible', 'off');
                plot(flip(data.ppm),squeeze(abs(data.X(8,8,:))));
                spec_path = fullfile(slice_path, 'MRS.png');
                saveas(fig, spec_path);
            end
        end
    end
end

seq_params = struct();
seq_params.sw = 2000;
seq_params.B0 = 3;
seq_params.lw = 2;
seq_params.te = 30E-3;
seq_params.tau1 = seq_params.te/4;
seq_params.tau2 = seq_params.te/4;
seq_params.n = 1024;
seq_params.sequence = 'press';
setup_params_list = {struct('voxels', [16, 16], 'snr', Inf, 'beta', 0, 'B0_map', ones([16,16]), 'slice_thickness', 10, 'seq_params', seq_params, 'num', 2)};
orientations = {'Axial', 'Coronal', 'Sagittal'};
permutations = {[1,2,3], [1,3,2], [2,3,1]};

simulate(setup_params_list, orientations, permutations);




