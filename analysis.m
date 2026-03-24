clear;
addpath(genpath('../FID-A-master'));
accel = 2;
rc_params = struct();
rc_params.rho = [1];
rc_params.ptol = [1e-3];
rc_params.dtol = [1e-3];
rc_params.mu = [10];
rc_params.gu = [2];
rc_params.gl = [2];
rc_params.verbose = false;
data_paths = {...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Axial\Axial_slice_51-60', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Coronal\Coronal_slice_21-30', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Sagittal\Sagittal_slice_41-50', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Axial\Axial_slice_21-30', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Coronal\Coronal_slice_1-10', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Sagittal\Sagittal_slice_1-10', ...
};
dwtmode('per','nodisp');
dwtmode;
method = @dwt1;
method_name = 'W_f';
analyses_path = fullfile(pwd,'analyses');
setup = 'setup_1';
basis_path = fullfile(pwd,'simulations',setup,'basis','basis.basis');
setup_path = fullfile(analyses_path, setup);
method_path = fullfile(setup_path, method_name);
if strcmp(method_name, 'FULL')
    accel_path = method_path;
    [~,~] = mkdir(accel_path);
elseif strcmp(method_name, 'ZF')
    accel_path = fullfile(method_path, sprintf('accel_%d', accel));
else
    accel_path = fullfile(method_path, sprintf('accel_%d', accel));
    tune_path = fullfile(accel_path, 'tune.mat');
    tune = load(tune_path).tune;
    optimal_lambda = tune.optimal_lambda;
    disp(optimal_lambda);
    rc_params.lambda = optimal_lambda;
end
for idx = 1:numel(data_paths)
    data_path = data_paths{idx};
    [~, name, ~] = fileparts(data_path);
    folder_path = fullfile(accel_path, name);
    [~,~] = mkdir(folder_path);
    data_path = fullfile(data_path, 'data.mat');
    data = load(data_path).data;
    K = data.K;
    us = Undersampler();
    [Ku, U] = us.undersample(K, accel, 'poisson_disc');
    if strcmp(method_name, 'FULL')
        Krc = K;
    elseif strcmp(method_name, 'ZF')
        Krc = Ku;
    else
        Krc = method(Ku, U, rc_params);
    end
    Xrc = fftshift(fftn(Krc));
    Xref = data.X;
    ascale = (Xrc(:)' * Xref(:)) / (Xrc(:)' * Xrc(:));
    res = Xref - ascale * Xrc;
    nrmse = norm(res(:)) / norm(Xref(:));
    result = struct();
    result.nrmse = nrmse;
    disp(result.nrmse);
    [Nx,Ny,Ns] = size(Xrc);
    mfields = cellstr(["Asp","NAAG","Cr","PCr","GPC","PCh", ...
        "Glu","Gln","GABA","GSH","NAA","Tau","PE", "Lac", "Ins"]);
    dmfields = strcat('d_', mfields);
    result.mfields = mfields;
    result.snrs = zeros(Nx,Ny);
    result.concs = zeros(Nx,Ny,numel(mfields));
    result.dconcs = zeros(Nx,Ny,numel(dmfields));
    fids = fftshift(fftshift(fft(fft(Krc, [], 1), [], 2), 1), 2);
    for nx = 1:Nx
        for ny = 1:Ny
            lcm = LCMHelper();
            lcm_path = fullfile(folder_path, 'lcm');
            [~,~] = mkdir(lcm_path);
            filbas = lcm.wsl_path(basis_path);
            filraw = fullfile(lcm_path,sprintf("%d_%d.RAW",nx,ny));
            fid = squeeze(fids(nx,ny,:));
            io_writelcm(fid,filraw,data.setup_params.seq_params.te*1E3);
            filraw = lcm.wsl_path(filraw);
            filps = lcm.wsl_path(fullfile(lcm_path, sprintf("%d_%d.ps",nx,ny)));
            filcon = fullfile(lcm_path, sprintf("%d_%d.control",nx,ny));
            lcm.write_control(filcon, filraw, filbas, filps);
            filcon = lcm.wsl_path(filcon);
            [status, cmdout] = system(['wsl -d Ubuntu-24.04 bash -c "export PATH=\$PATH:/home/rayne/.lcmodel/bin' ...
                        '    && cd /mnt ' ...
                        '    && lcmodel < ' char(filcon) '"']);
            % [~,~] = mkdir(pdf_path);
            % eps2pdf(char(fullfile(lcm_path, sprintf("%d_%d.ps",nx,ny))), char(fullfile(pdf_path,  sprintf("%d_%d.pdf",nx,ny))));
            out = io_readlcmtab(fullfile(lcm_path, sprintf("%d_%d.table",nx,ny)));
            result.snrs(nx,ny) = out.SNR;
            result.concs(nx,ny,:) = cellfun(@(a) out.(a), cellstr(mfields));
            result.dconcs(nx,ny,:) = cellfun(@(a) out.(a), cellstr(dmfields));
            disp(result);
        end
    end
    result.K = Krc;
    result.X = Xrc;
    result_path = fullfile(folder_path, 'result.mat');
    save(result_path, 'result');
end