clear;
addpath(genpath("../FID-A-master"));
cfg = struct();
cfg.metabs = readtable('data/metabolites/metab_df.csv');
cfg.metabs_list = ["NAA","NAAG","Cr","PCr"];
% cfg.metabs_list = ["NAA","NAAG","Cr","PCr","GPC","PCh", ...
%     "Glu","Gln","GABA","GSH","Asp","Tau","PE", "Lac"];
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
voxrange = [8,8];
sigma = 1000;
params_order = {'n','sw','B0','lw','sys','tau1','tau2'};
beta = 0.0003;
[kspace, t, imspace, ppm] = sm.extract_kspace(phantom,@sim_press,seq_params, params_order, zslice, voxrange, sigma, beta);
figure(1);
clf;
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
        xlim([min(ppm), max(ppm)])
        box on;
        ax.Color = 'None';
        ax.XTick = [];
        ax.YTick = [];
        set(ax,'XDir','reverse' )
    end
end
figure(2);
clf;
max_val = rms(abs(kspace(:)));
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
    end
end
p = 0.8;
mask = rand(voxrange(1), voxrange(2)) > (1-p);
mask = repmat(mask, 1, 1, seq_params.n);
Ku = mask .* kspace;
figure(3);
clf;
for r = 1:voxrange(1)
    for c = 1:voxrange(2)
        x = (c-1)/voxrange(1);
        y = 1-r/voxrange(2);
        ax = axes('Units','normalized',...
                      'Position',[x y w h]);
        plot(t, abs(squeeze(Ku(r,c,:)) / max_val), 'g');
        ylim([-0.1, 2])
        xlim([min(t), max(t)])
        box on;
        ax.Color = 'None';
        ax.XTick = [];
        ax.YTick = [];      
    end
end


function out = f(P, Ku, mask)
    temp = mask .* ifftn(P) - Ku;
    out = 0.5 * norm(temp(:), 2)^2;
end

function out = g(Q, lambda)
    out = lambda *  norm(Q(:), 1);
end

function [Q, B] = psi(P)
    Ns = size(P,3);
    Q = cell(1, Ns);
    B = cell(1, Ns);
    for ns = 1:Ns
        P_slice = P(:,:,ns);
        [Q_slice, B_slice] = wavedec2(P_slice, 3, 'db2');
        Q{ns} = Q_slice;
        B{ns} = B_slice;
    end
    Q = cat(3, Q{:});
    B = cat(3, B{:});
end

function P = psi_adj(Q, B)
    Ns = size(Q);
    Ns = Ns(end);
    P = cell(1, Ns);
    idx_Q = repmat({':'}, 1, ndims(Q));
    idx_B = repmat({':'}, 1, ndims(B));
    for ns = 1:Ns
        idx_Q{end} = ns;
        Q_slice = Q(idx_Q{:});
        idx_B{end} = ns;
        B_slice = B(idx_B{:});
        P{ns} = waverec2(Q_slice, B_slice, 'db2');
    end
    P = cat(3, P{:});
end

P = fftn(Ku);
max_val = max(abs(P(:)));
P = P / max_val;
Ku = Ku / max_val;
[Q,B] = psi(P);
W = zeros(size(Q));

function out = A_pcg(v, mask, rho)
    out = reshape(v, size(mask));
    [temp1, temp2] = psi(out);
    out = fftn(mask .* ifftn(out)) + rho * psi_adj(temp1, temp2);
    out = out(:);
end
rho = 0.01;
lambda = 1;
disp(f(P, Ku, mask));
disp(g(Q, lambda));
% for iter = 1:2
%     rhs = fftn(Ku) + rho * psi_adj(Q - W, B);  
%     rhs = rhs(:);
%     p0 = P(:);
%     A = @(v) A_pcg(v, mask, rho);
%     p = pcg(A, rhs, 0.01, 50, [], [], p0);
%     P = reshape(p, size(Ku));
%     temp = psi(P) + W;
%     q = wthresh(temp(:), 's', lambda/rho);
%     Q = reshape(q, size(Q));
% 
%     disp(f(P, Ku, mask));
%     disp(g(Q, lambda));
% 
%     W = W+psi(P)-Q;
% end



P = fftshift(P);
figure(4);
clf;
ax_bg = axes('Units','normalized','Position',[0 0 1 1]);
img = niftiread('data/skeleton/labels.nii');
img = mat2gray(img(:,:,zslice));
imshow(img, 'Parent', ax_bg);
axis(ax_bg, 'off');
set(ax_bg,'DataAspectRatioMode','auto')
w = 1/voxrange(1);
h = 1/voxrange(2);
max_val = max(abs(P(:)));
for r = 1:voxrange(1)
    for c = 1:voxrange(2)
        x = (c-1)/voxrange(1);
        y = 1-r/voxrange(2);
        ax = axes('Units','normalized',...
                      'Position',[x y w h]);
        plot(ppm, flip(abs(squeeze(P(r,c,:)) / max_val)), 'g');
        ylim([-0.1, 1])
        xlim([min(ppm), max(ppm)])
        box on;
        ax.Color = 'None';
        ax.XTick = [];
        ax.YTick = [];
        set(ax,'XDir','reverse' )
    end
end









