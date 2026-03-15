





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
    Z = wavedec3(X,3,'db2');
end

function X = psiH(Z)
    X = waverec3(Z);
end

function Z = psi_add(Z1,Z2)
    Z = Z1;
    Z.dec = cellfun(@plus, Z1.dec, Z2.dec, 'UniformOutput', false);
end

function Z = psi_sub(Z1,Z2)
    Z = Z1;
    Z.dec = cellfun(@minus, Z1.dec, Z2.dec, 'UniformOutput', false);
end

function Z = psi_init(Z0)
    Z = Z0;
    Z.dec = cellfun(@(z) zeros(size(z)), Z0.dec, 'UniformOutput', false);
end

function Z = psi_sthresh(Z0,T)
    Z = Z0;
    Z.dec = cellfun(@(z) reshape(wthresh(z(:),'s',T),size(z)), Z0.dec, 'UniformOutput', false);
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
    [xnew, flag] = pcg(A,rhs,1e-6,50,[],[],x0);
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
    % temp = psi_sub(psi(X),Z);
    % temp = temp.dec;
    % out1 = sqrt(sum(cellfun(@(z) sum(abs(z(:)).^2), temp)));
    % out2 = sum(cellfun(@(z) sum(abs(z),'all'), Z.dec));
    % temp = phi(X,U) - Y;
    % out3 = 0.5 * norm(temp(:))^2;
end
% Y = U .* kspace;
% disp(rms(Y(:)) / rms(kspace(:)));
% disp(nnz(U)/numel(U));
% 
% Y_norm = norm(Y(:));
% Y = Y/Y_norm;
% X0 = phiH(Y,U);
% Z0 = psi(X0);
% W0 = psi_init(Z0);
% rho = 1;
% lambda = 0.00001;
% max_iter = 50;
% rc = Reconstruction();
% [X,Z,W] = rc.reconstruct(X0, Z0, W0, U, Y, rho, lambda, @updateX, @updateZ, @updateW, max_iter, @monitor);
% X = fftshift(X) * Y_norm * sqrt(numel(U));
% vis.visualize(X, flip(ppm), 3, bg_img, max(abs(imspace(:))), 'reverse', [0, 1]);



