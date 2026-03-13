clear;
addpath(genpath("../FID-A-master")); addpath(genpath(".."));
cfg = struct(); cfg.metabs = readtable('data/metabolites/metab_df.csv');
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
zslice = 90;
voxrange = [8,8];
sigma = 0;
params_order = {'n','sw','B0','lw','sys','tau1','tau2'};
beta = 0.0003;
[kspace, t, imspace, ppm] = sm.extract_kspace(phantom,@sim_press,seq_params, params_order, zslice, voxrange, sigma, beta);
bg_img = niftiread('data/skeleton/labels.nii');
bg_img = mat2gray(bg_img(:,:,zslice));
vis = Visualizer();
vis.visualize(kspace, t, 1, bg_img, rms(abs(kspace(:))), 'normal', [0, 2]);
vis.visualize(imspace, flip(ppm), 2, bg_img, max(abs(imspace(:))), 'reverse', [0, 1]);
accel = 2;
[Nx,Ny,Nt] = size(kspace);
Nsamp = Nx*Ny/accel;
U_slice = rand(size(kspace,[1,2])) .* rms(kspace,3);
[~,idx] = sort(U_slice(:), 'descend');
idx = idx(1:Nsamp);
U_slice = false(size(U_slice));
U_slice(ind2sub(size(U_slice), idx)) = true;
U = repmat(U_slice, 1, 1, Nt);
% for nt = 1:Nt
%     U_slice = rand(size(kspace,[1,2])) .* abs(kspace(:,:,nt));
%     [~,idx] = sort(U_slice(:), 'descend');
%     Nsamp = round(Nx*Ny / accel);
%     idx = idx(1:Nsamp);
%     U_slice = false(size(U_slice));
%     U_slice(ind2sub(size(U_slice), idx)) = true;
%     U(:,:,nt) = U_slice;
% end
% U = rand(size(kspace)) .* abs(kspace);
% 
% 
% 
% 
% [~,idx] = sort(U(:), 'descend');
% idx = idx(1:Nsamp);
% U = false(size(U));
% U(ind2sub(size(U), idx)) = true;





function Y = phi(X,U)
    Y = U .* ifftn(X) * sqrt(numel(U));
end

function X = phiH(Y,U)
    X = fftn(U .* Y) / sqrt(numel(U));
end

function Z = psi(X)
    Z = cell(2,1);
    [temp,B] = wavedec2(real(X(:,:,1)),3,'db2');
    Ns = size(X,3);
    Z{1} = zeros(size(temp,1),size(temp,2),Ns);
    for ns = 1:Ns
        Z{1}(:,:,ns) = wavedec2(real(X(:,:,ns)),3,'db2') + 1i * wavedec2(imag(X(:,:,ns)),3,'db2');
    end
    Z{2} = B;
end

function X = psiH(Z)
    temp = waverec2(real(Z{1}(:,:,1)),Z{2},'db2');
    Ns = size(Z{1},3);
    X = zeros(size(temp,1),size(temp,2),Ns);
    for ns = 1:Ns
        X(:,:,ns) = waverec2(real(Z{1}(:,:,ns)),Z{2},'db2') + 1i * waverec2(imag(Z{1}(:,:,ns)),Z{2},'db2');
    end
end

function Z = psi_add(Z1,Z2)
    Z = Z1;
    Z{1} = Z1{1} + Z2{1};
end

function Z = psi_sub(Z1,Z2)
    Z = Z1;
    Z{1} = Z1{1} - Z2{1};
end

function Z = psi_init(Z0)
    Z = Z0;
    Z{1} = zeros(size(Z0{1}));
end

function Z = psi_sthresh(Z0,T)
    Z = Z0;
    temp = wthresh(Z0{1}(:), 's', T);
    Z{1} = reshape(temp, size(Z0{1}));
end

function xnew = A_pcg(x,U,rho)
    temp = reshape(x,size(U));
    Xtemp = phiH(phi(temp,U),U) + rho * psiH(psi(temp));
    xnew = Xtemp(:);
end

function Xnew = updateX(X,Z,W,U,Y,rho)
    rhs = phiH(Y,U) + rho*psiH(psi_sub(Z,W));
    rhs = rhs(:);
    A = @(x) A_pcg(x,U,rho);
    X0 = X;
    x0 = X0(:);
    xnew = pcg(A,rhs,1e-6,50,[],[],x0);
    Xnew = reshape(xnew,size(X));
end

function Znew = updateZ(X,Z,W,rho,lambda)
    temp = psi_add(psi(X), W);
    Znew = psi_sthresh(temp, lambda / rho);
end
function Wnew = updateW(X,Z,W)
    Wnew = psi_add(W, psi_sub(psi(X), Z));
end
function monitor(X,Z,W,U,Y)


end
Y = U .* kspace;
disp(rms(Y(:)) / rms(kspace(:)));
disp(nnz(U)/numel(U));

Y_norm = norm(Y(:));
Y = Y/Y_norm;
X0 = phiH(Y,U);
Z0 = psi(X0);
W0 = psi_init(Z0);
rho = 1;
lambda = 0.001;
max_iter = 5;
rc = Reconstruction();
[X,Z,W] = rc.reconstruct(X0, Z0, W0, U, Y, rho, lambda, @updateX, @updateZ, @updateW, max_iter, @monitor);
X = fftshift(X) * Y_norm * sqrt(numel(U));
vis.visualize(X, flip(ppm), 3, bg_img, max(abs(imspace(:))), 'reverse', [0, 1]);



