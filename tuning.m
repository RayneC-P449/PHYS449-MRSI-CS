clear;
accel = 4;
setup ='setup_1';
data_paths = {...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Axial\Axial_slice_1-10', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Coronal\Coronal_slice_51-60', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Sagittal\Sagittal_slice_31-40', ...
};
dwtmode('per','nodisp');
dwtmode;
method = @hkm1;
method_name = 'HKM';
analyses_path = fullfile(pwd,'analyses');
setup_path = fullfile(analyses_path, setup);
rc_params = struct();
rc_params.rho = [1,1];
rc_params.ptol = [1e-2,1e-2];
rc_params.dtol = [1e-2,1e-2];
rc_params.mu = [10,10];
rc_params.gu = [2,2];
rc_params.gl = [2,2];
rc_params.verbose = true;
[~,~] = mkdir(analyses_path);
[~,~] = mkdir(setup_path);
method_path = fullfile(setup_path, method_name);
[~,~] = mkdir(method_path);
accel_path = fullfile(method_path, sprintf('accel_%d', accel));
[~,~] = mkdir(accel_path);
tune_path = fullfile(accel_path, 'tune.mat');
lambdas = num2cell(logspace(-9,-2,20));
% a = logspace(-7, -4, 3);
% b = logspace(-7, -4, 3);
% [A, B] = ndgrid(a, b);
% lambdas = num2cell([A(:), B(:)], 2);
% lambdas = {[1e-6, 1e-6]};
tune = struct();
tune.lambdas = {};    
tune.nrmse = [];
for i = 1:numel(lambdas)
    lambda = lambdas{i};
    disp(lambda);
    rc_params.lambda = lambda;
    tune.lambdas{end+1} = lambda;
    all_nrmse = [];
    for idx = 1:numel(data_paths)
        data_path = fullfile(data_paths{idx});
        parts = strsplit(data_path, filesep);
        folder_name = parts{end};
        data_path = fullfile(data_path, 'data.mat');
        data = load(data_path).data;
        K = data.K;
        us = Undersampler();
        [Ku, U] = us.undersample(K, accel, 'poisson_disc');
        % Krc = method(Ku, U, rc_params);
        % disp(numel(U)/nnz(U));
        % Krc = Ku;
        Xrc = fftshift(fftn(Krc));
        Xref = data.X;
        ascale = (Xrc(:)' * Xref(:)) / (Xrc(:)' * Xrc(:));
        res = Xref - Xrc;
        nrmse = norm(res(:)) / norm(Xref(:));
        disp(nrmse);
        all_nrmse(end+1) = nrmse;
    end
    avg_nrmse = mean(all_nrmse);
    tune.nrmse(end+1) = avg_nrmse;
    disp(tune.nrmse);
end
[~,idx] = min(tune.nrmse);
tune.optimal_lambda = tune.lambdas{idx};
tune.minimal_nrmse = tune.nrmse(idx);
save(tune_path, 'tune');



