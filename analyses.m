% path = fullfile(pwd, 'simulated','1','81');
% bg_img =  mat2gray(imread(fullfile(path,'MRI.png')));
% data = load(fullfile(path,'kspace.mat')).data;
% vis = Visualizer;
% imspace = fftshift(fftn(data.kspace));
% vis.visualize(imspace, flip(data.ppm), 1, bg_img, max(abs(imspace(:))), 'reverse', [0 1]);
% s = abs(squeeze(imspace(8,8,:)));
% plot(flip(data.ppm), s);

clear;
lcm = LCMHelper();
% lcm.make_basis();
vis = Visualizer();


fields = cellstr(["Asp","NAAG","Cr","PCr","GPC","PCh", ...
    "Glu","Gln","GABA","GSH","NAA","Tau","PE", "Lac", "Ins"]);
dfields = strcat('d_', fields);

addpath(genpath('../FiD-A-master'));

function u = sample(Rx, Ry, sp, nx, ny)
    c = 2;
    lx = ceil((nx+1)/2-c):floor((nx+1)/2+c);
    ly = ceil((ny+1)/2-c):floor((ny+1)/2+c);
    u = false(nx,ny);
    u(lx,ly) = 1;
    pxs = zeros(nx*ny, 1);
    pys = zeros(nx*ny, 1);
    pxs(1) = randi(nx);
    pys(1) = randi(ny);
    num_actives = 1;
    while num_actives >= 1
        i = randi(num_actives);
        px = pxs(i);
        py = pys(i);
        rx = Rx(px,py);
        ry = Ry(px,py);
        done = false;
        k = 1;
        while ~done && k <= sp
            v = sqrt(rand(1) * 3 + 1);
            t = 2 * pi * rand(1);
            qx = px + v * rx * cos(t);
            qy = py + v * ry * sin(t);
            if qx>=1 && qx < nx + 1 && qy >=1 && qy < ny + 1
                startx = max(fix(qx - rx), 1);
                endx = min(fix(qx+rx+1),nx);
                starty = max(fix(qy - ry), 1);
                endy = min(fix(qy+ry+1),ny);
                done = true;
                for x = startx:endx
                    for y = starty:endy
                        if (u(x,y) == 1 && (((qx-x)/Rx(x,y))^2 + ((qy-y)/Ry(x,y))^2 < 1))                            
                            done = false;
                            break;
                        end
                    end
                end
            end
            k = k+1;
        end
        if done
            num_actives = num_actives + 1;
            pxs(num_actives) = fix(qx);
            pys(num_actives) = fix(qy);
            u(fix(qx), fix(qy)) = 1;
        else
            pxs(i) = pxs(num_actives);
            pys(i) = pys(num_actives);
            num_actives = num_actives - 1;
        end
    end
end
function u = poisson_disc(nx, ny, accel, sp, tol)
    x = -(nx-1)/2:(nx-1)/2;
    y = -(ny-1)/2:(ny-1)/2;
    [x,y] = ndgrid(x,y);
    r = sqrt(x.^2 + y.^2);
    smin = 0;
    smax = max(nx,ny);
    while smin < smax
        s = (smax + smin) / 2;
        Rx = (1 + s*r) * nx / max(nx, ny);
        Ry = (1 + s*r) * ny / max(nx, ny);
        u = sample(Rx, Ry, sp, nx, ny);
        actual_accel = (nx * ny) / nnz(u);
        if abs(actual_accel - accel) < tol
            break;
        end
        if actual_accel < accel
            smin = s;
        else
            smax = s;
        end
    end
end
function [Ku,U] = undersample(kspace, accel)
    [nx,ny,nt] = size(kspace);
    U = false(size(kspace));
    for kt = 1:nt
        U(:,:,kt) = poisson_disc(nx, ny, accel, 20, 0.05);
    end
    Ku = kspace;
    Ku(~U) = 0;
    clf;
    nslices = 3;
    indices = round(linspace(1,size(U,3),nslices));
    slice(double(permute(U,[2 3 1])), indices,[], []);
    drawnow;
end

function Krc = reconstruct(K, scheme, params, accel)
    switch scheme
        case 'full'
            Krc = K;
        case 'dwt1'
            [Ku, U] = undersample(K, accel);
            Krc = dwt1(Ku, U, params);
        case 'hkm1'
            [Ku, U] = undersample(K, accel);
            Krc = hkm1(Ku, U, params);
        case 'zfll'
            [Ku, U] = undersample(K, accel);
            Krc = Ku;
    end
end
spath = fullfile(pwd, 'simulated');
rcscheme = 'full';
accel = 3;
rcpath = fullfile(pwd, 'reconstructed', 'full');
params = struct;
params.lambda = [1e-4];
% 1e-5, 1e-3 wow
% 1e-4 might be good
params.rho = [1];
params.mu = 10;
params.gu = 2;
params.gl = 2;
params.ptol = [1e-3];
params.dtol = [1e-3];
params.label = '1';
cspath = fullfile(pwd, 'reconstructed', rcscheme, params.label, num2str(accel));
file_per_orient = 1;
for o = 3:3
    folders = fullfile(spath, num2str(o));
    folders = dir(folders); 
    fcount = 0;
    for f = 1:length(folders)
        if folders(f).isdir && ~strcmp(folders(f).name, '.') && ~strcmp(folders(f).name, '..')
            fcount = fcount + 1;
            folder = folders(f);
            kpath = fullfile(folder.folder, folder.name);
            data = load(fullfile(kpath, 'kspace.mat')).data;
            data.te = 30E-3;
            te = data.te * 1E3;
            kspace = data.kspace;
            kspace = kspace / norm(kspace, 'fro');
            Xref = fftshift(fftn(kspace));
            kspace = reconstruct(kspace, rcscheme, params, accel);
            % kspace = kspace / norm(kspace, 'fro');
            Xrc = fftshift(fftn(kspace));
            temp = fullfile(cspath, num2str(o), folder.name);
            [~,~] = mkdir(temp);
            temp = fullfile(temp, 'MRSI.png');
            vis.visualize(flip(data.ppm),Xrc,'reverse',[0,max(abs(Xrc),[],'all')],data.MRI, temp);

            fids = fftshift(fftshift(fft(fft(kspace,[], 1),[],2), 1), 2);
            imageData = struct;
            imageData.rmse = rms(Xrc - Xref, 'all');
            disp(imageData.rmse);
            imageData.specs = fftshift(fftn(kspace));
            imageData.snrs = zeros(size(kspace, [1, 2]));
            imageData.concs = zeros(size(kspace,1), size(kspace,2), numel(fields));
            imageData.dconcs = zeros(size(kspace,1), size(kspace,2), numel(dfields));
            for vx = 1:size(kspace, 1)
                for vy = 1:size(kspace, 2)
                    vpath = fullfile(cspath, num2str(o), folder.name, num2str(vx) + "-" + num2str(vy));
                    [~,~] = mkdir(vpath);
                    temp = squeeze(fids(vx,vy,:));
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
                    out = io_readlcmtab(fullfile(vpath, 'out.table'));
                    imageData.snrs(vx,vy) = out.SNR;
                    imageData.concs(vx,vy,:) = cellfun(@(a) out.(a), cellstr(fields));
                    imageData.dconcs(vx,vy,:) = cellfun(@(a) out.(a), cellstr(dfields));
                    eps2pdf(char(fullfile(vpath, 'out.ps')), char(fullfile(vpath, 'out.pdf')));
                    disp(vpath);
                end
            end
            if fcount == file_per_orient
                break;
            end
        end
    end
end
% 
% 
% 



