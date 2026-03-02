clear;
addpath(genpath("../FID-A-master"));
cfg = struct();
cfg.metabs = readtable('data/metabolites/metab_df.csv');
% cfg.metabs_list = ["NAA","NAAG","Cr","PCr"];
cfg.metabs_list = ["NAA","NAAG","Cr","PCr","GPC","PCh", ...
    "Glu","Gln","GABA","GSH","Asp","Tau","PE", "Lac"];
cfg.labels = niftiread('data/skeleton/labels.nii');
mm_json = jsondecode(fileread('data/macromolecules/mm_json.json'));
cfg.mm_list = mm_json.names;
cfg.mm = mm_json;
cfg.water = zeros(2,2);
cfg.water(1,:) = [36000, 63.4E-3];
cfg.water(2,:) = [43300, 52.6E-3];

phantom = Phantom(cfg);
sm = SignalModel();
seq_func = @sim_press;
seq_params = struct();
seq_params.n = 1024;
seq_params.sw = 2000;
seq_params.B0 = 3;
seq_params.lw = 4;
seq_params.sys = 0; 
TE   = 30e-3;  
seq_params.tau1 = TE/4;             
seq_params.tau2 = TE/4;
zslice = 90:90;
voxrange = [16,16];
sigma = 1000;
params_order = {'n','sw','B0','lw','sys','tau1','tau2'};
beta = 0.0003;
[kspace, t, imspace, ppm] = sm.extract_kspace(phantom,@sim_press,seq_params, params_order, zslice, voxrange, sigma, beta);
figure(1);
clf;
[Nx, Ny, Nz] = size(phantom.labels);
ax_bg = axes('Units','normalized','Position',[0 0 1 1]);
img = niftiread('data/skeleton/labels.nii');
img = mat2gray(img(:,:,zslice));
imshow(img, 'Parent', ax_bg);
axis(ax_bg, 'off');
set(ax_bg,'DataAspectRatioMode','auto')
w = 1/voxrange(1);
h = 1/voxrange(2);
max_val = max(abs(imspace(:)));
for r = 1:voxrange(1)
    for c = 1:voxrange(2)
        x = (c-1)/voxrange(1);
        y = 1-r/voxrange(2);
        ax = axes('Units','normalized',...
                      'Position',[x y w h]);
        plot(ppm, flip(abs(squeeze(imspace(r,c,:)) / max_val)), 'g');
        ylim([-0.1, 1])
        box on;
        ax.Color = 'None';
        ax.XTick = [];
        ax.YTick = [];
        set (ax,'XDir','reverse' )
    end
end
figure(2);
clf;
max_val = rms(abs(kspace(:)));
mask = zeros(voxrange(1), voxrange(2), seq_params.n);
for r = 1:voxrange(1)
    for c = 1:voxrange(2)
        x = (c-1)/voxrange(1);
        y = 1-r/voxrange(2);
        ax = axes('Units','normalized',...
                      'Position',[x y w h]);
        plot(t, abs(squeeze(kspace(r,c,:)) / max_val), 'g');
        ylim([-0.1, 2])
        xlim([min(t), max(t)])
        box on;
        ax.Color = 'None';
        ax.XTick = [];
        ax.YTick = [];
        mask(r,c,:) = rand(seq_params.n, 1) > 0.5;
        temp = logical(mask(r,c,:));
        temp = temp(:);        
    end
end

