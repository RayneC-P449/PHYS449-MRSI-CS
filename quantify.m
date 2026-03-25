clear;
addpath(genpath('../FID-A-master'));
accel = 2;
method_name = 'HKM';
analyses_path = fullfile(pwd,'analyses');
setup = 'setup_1';
method_path = fullfile(analyses_path, setup, method_name, sprintf('accel_%d',accel));
full_path = fullfile(analyses_path,setup, 'FULL');

data_folders = dir(method_path);
% Loop through folders
for k = 1:numel(data_folders)
    if data_folders(k).isdir && ~strcmp(data_folders(k).name, '.') && ~strcmp(data_folders(k).name, '..')
        method_data = fullfile(method_path, data_folders(k).name);
        name = data_folders(k).name;
        full_data = fullfile(full_path, name);
        method_result = load(fullfile(method_data,'result.mat')).result;
        full_result = load(fullfile(full_data,'result.mat')).result;
        NAA_idx = 11;
        nrmse = method_result.nrmse;
        
        full_conc = full_result.concs / sum(full_result.concs(:,:,NAA_idx),'all');
        method_conc = method_result.concs / sum(method_result.concs(:,:,NAA_idx),'all');
        

        
        conc_err = abs(full_conc - method_conc);
        conc_err = mean(conc_err,'all');
        % disp(conc_err);
        disp(method_result.nrmse);
        
    end
end