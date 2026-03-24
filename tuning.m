clear;
accel = 2;
data_paths = {...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Axial\Axial_slice_1-10', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Coronal\Coronal_slice_51-60', ...
    'C:\Users\rayne\PHYS449\Thesis\simulations\setup_1\Sagittal\Sagittal_slice_31-40', ...
};
dwtmode('per','nodisp');
dwtmode;
method = @dwt1;
method_name = 'W_f';
setup ='setup_1';
analyses_path = fullfile(pwd,'analyses');
setup_path = fullfile(analyses_path, setup);
rc_params = struct();
rc_params.rho = [1];
rc_params.ptol = [1e-2];
rc_params.dtol = [1e-2];
rc_params.mu = [10];
rc_params.gu = [2];
rc_params.gl = [2];
rc_params.verbose = false;
[~,~] = mkdir(analyses_path);
[~,~] = mkdir(setup_path);
method_path = fullfile(setup_path, method_name);
[~,~] = mkdir(method_path);
accel_path = fullfile(method_path, sprintf('accel_%d', accel));
[~,~] = mkdir(accel_path);
tune_path = fullfile(accel_path, 'tune.mat');
lambdas = num2cell(logspace(-9,-2,16));
tune = struct();
tune.lambdas = {};    
tune.nrmse = [];
% for i = 1:numel(lambdas)
%     lambda = lambdas{i};
%     disp(lambda);
%     rc_params.lambda = lambda;
%     tune.lambdas{end+1} = lambda;
%     all_nrmse = [];
%     for idx = 1:numel(data_paths)
%         data_path = fullfile(data_paths{idx});
%         parts = strsplit(data_path, filesep);
%         folder_name = parts{end};
%         data_path = fullfile(data_path, 'data.mat');
%         data = load(data_path).data;
%         K = data.K;
%         us = Undersampler();
%         [Ku, U] = us.undersample(K, accel, 'poisson_disc');
%         Krc = method(Ku, U, rc_params);
%         Xrc = fftshift(fftn(Krc));
%         Xref = data.X;
%         ascale = (Xrc(:)' * Xref(:)) / (Xrc(:)' * Xrc(:));
%         res = Xref - ascale * Xrc;
%         nrmse = norm(res(:)) / norm(Xref(:));
%         disp(nrmse);
%         all_nrmse(end+1) = nrmse;
%     end
%     avg_nrmse = mean(all_nrmse);
%     tune.nrmse(end+1) = avg_nrmse;
%     disp(tune.nrmse);
% end
% [~,idx] = min(tune.nrmse);
% tune.optimal_lambda = tune.lambdas{idx};
% tune.minimal_nrmse = tune.nrmse(idx);
tune = struct();
tune.optimal_lambda = 4.6416e-05;
disp(tune.optimal_lambda);
save(tune_path, 'tune');



