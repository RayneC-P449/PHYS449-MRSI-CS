clear;
addpath(genpath('../FID-A-master'));
accel = 2;
full_name = "full";
method_name = 'W_f';
analyses_path = fullfile(pwd,'analyses');
setup = 'setup_1';
method_path = fullfile(analyses_path, setup, method_name, sprintf('accel_%d',accel));
full_path = fullfile(analyses_path,setup, 'FULL');

data_folders = dir(method_path);
% Loop through folders
for k = 1:length(data_folders)
    if data_folders(k).isdir && ~strcmp(data_folders(k).name, '.') && ~strcmp(data_folders(k).name, '..')
        method_data = fullfile(method_path, data_folders(k).name);
        name = data_folders(k).name;
        full_data = fullfile(full_path, name);
        method_result = load(fullfile(method_data,'result.mat')).result;
        full_result = load(fullfile(full_data,'result.mat')).result;

        disp(method_result);
        disp(full_result);
    end
end