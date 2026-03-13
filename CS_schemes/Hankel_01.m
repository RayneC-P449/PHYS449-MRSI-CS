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
voxrange = [9,9];
sigma = 0;
params_order = {'n','sw','B0','lw','sys','tau1','tau2'};
beta = 0.0003;
[kspace, t, imspace, ppm] = sm.extract_kspace(phantom,@sim_press,seq_params, params_order, zslice, voxrange, sigma, beta);
bg_img = niftiread('data/skeleton/labels.nii');
bg_img = mat2gray(bg_img(:,:,zslice));
vis = Visualizer();
vis.visualize(kspace, t, 1, bg_img, rms(abs(kspace(:))), 'normal', [0, 2]);
vis.visualize(imspace, flip(ppm), 2, bg_img, max(abs(imspace(:))), 'reverse', [0, 1]);
accel = 3;
[Nx,Ny,Nt] = size(kspace);
U = zeros(size(kspace));
Nsamp = round(Nx*Ny/accel);
[kx,ky] = ndgrid(linspace(-1,1,Nx), linspace(-1,1,Ny));
for nt = 1:Nt
    U_slice = rand(size(kspace,[1,2])) .* (1 - 0.2*(1/sqrt(2))*sqrt(kx.^2 + ky.^2));
    [~,idx] = sort(U_slice(:), 'descend');
    idx = idx(1:Nsamp);
    U_slice = false(size(U_slice));
    U_slice(ind2sub(size(U_slice), idx)) = true;
    U(:,:,nt) = U_slice;
end


function Y = phi(X,U)
    Y = U .* X;
end

function X = phiH(Y,U)
    X = U .* Y;
end

function Z = psi(X)
    L = numel(X) / 2;
    Z = hankel(X(1:L), X(L:end));
end

function X = psiH(Z)
    [L, Lp1] = size(Z); 
    n = L + Lp1 - 1;    
    X = zeros(n,1);
    for i = 1:L
        for j = 1:Lp1
            k = i + j - 1;
            X(k) = X(k) + Z(i,j);
        end
    end
end

function xnew = A_pcg(x, U, rho)
    xnew = U .* x + rho * psiH(psi(x));
end
function Xnew = updateX(X,Z,W,U,Y,rho)
    rhs = Y + rho*psiH(Z - W);
    rhs = rhs(:);
    A = @(x) A_pcg(x, U, rho);
    X0 = X;
    x0 = X0(:);
    [xnew, flag] = pcg(A, rhs, 1e-6, 50, [], [], x0);
    Xnew = reshape(xnew, size(X));
end

function Znew = updateZ(X,Z,W,rho,lambda)
    [U,S,V] = svd(psi(X) + W,'econ');
    S = max(S - lambda/rho, 0);
    Znew = U * S * V';
end

function Wnew = updateW(X,Z,W)
    Wnew = W + psi(X) - Z;
end

function out = monitor(X,Z,W,U,Y)

end

rc = Reconstruction();
X_final = zeros(size(kspace));
for nx = 1:Nx
    for ny = 1:Ny
        U_line = U(nx,ny,:);
        Y = U_line .* kspace(nx,ny,:);
        Y_norm = norm(Y(:));
        Y = Y / Y_norm;
        Y = Y(:);
        U_line = U_line(:);
        X0 = Y;
        Z0 = psi(X0);
        W0 = zeros(size(Z0));
        rho = 1;
        lambda = 0.00001;
        max_iter = 7;
        [X,Z,W] = rc.reconstruct(X0, Z0, W0, U_line, Y, rho, lambda, @updateX, @updateZ, @updateW, max_iter, @monitor);
        disp([nx, ny]);
        X_final(nx,ny,:) = X * Y_norm;
    end
end

temp = U .* kspace;
disp(rms(temp(:)) / rms(kspace(:)));

X_img = fftshift(fftn(X_final));
vis.visualize(X_img, flip(ppm), 3, bg_img, max(abs(imspace(:))), 'reverse', [0, 1]);
disp(U(:,:,1));
disp(U(:,:,2));

disp(sum(U~=0,3));