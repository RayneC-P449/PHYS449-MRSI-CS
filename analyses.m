% path = fullfile(pwd, 'simulated','1','81');
% bg_img =  mat2gray(imread(fullfile(path,'MRI.png')));
% data = load(fullfile(path,'kspace.mat')).data;
% vis = Visualizer;
% imspace = fftshift(fftn(data.kspace));
% vis.visualize(imspace, flip(data.ppm), 1, bg_img, max(abs(imspace(:))), 'reverse', [0 1]);
% s = abs(squeeze(imspace(8,8,:)));
% plot(flip(data.ppm), s);

lcm = LCMHelper();
% lcm.make_basis();




spath = fullfile(pwd, 'simulated');
rcpath = fullfile(pwd, 'reconstructed', 'full');
for o = 1:3
    folders = fullfile(spath, num2str(o));
    folders = dir(folders); 
    for f = 1:length(folders)
        if folders(f).isdir && ~strcmp(folders(f).name, '.') && ~strcmp(folders(f).name, '..')
            folder = folders(f);
            kpath = fullfile(folder.folder, folder.name);
            data = load(fullfile(kpath, 'kspace.mat')).data;
            data.te = 30E-3;
            te = data.te * 1E3;
            kspace = data.kspace;
            fids = fftshift(fftshift(fft(fft(kspace,[], 1),[],2), 1), 2);
            for vx = 1:size(kspace, 1)
                for vy = 1:size(kspace, 2)
                    vpath = fullfile(rcpath, num2str(o), folder.name, num2str(vx) + "-" + num2str(vy));
                    [~,~] = mkdir(vpath);
                    temp = squeeze(fids(vx,vy,:));
                    temppath = fullfile('reconstructed', 'full', num2str(o), folder.name, num2str(vx) + "-" + num2str(vy));
                    filraw = fullfile(vpath,'fid.RAW');
                    io_writelcm(temp,filraw,te*1E3);
                    filraw = lcm.wsl_path(filraw);
                    filps = lcm.wsl_path(fullfile(vpath, 'out.ps'));
                    filbas = lcm.wsl_path(fullfile(spath, 'basis', 'basis.basis'));
                    filcon = fullfile(vpath, 'lcm.control');
                    lcm.write_control(filcon, filraw, filbas, filps)
                    filcon = lcm.wsl_path(filcon);
                    [status, cmdout] = system(['wsl -d Ubuntu-24.04 bash -c "export PATH=\$PATH:/home/rayne/.lcmodel/bin' ...
                        '    && cd /mnt ' ...
                        '    && lcmodel < ' char(filcon) '"']);
                    disp(cmdout);
                end
            end
        end

    end
end


