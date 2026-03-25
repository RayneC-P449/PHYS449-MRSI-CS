function [Krc, fval, gval] = dwt1(Ku, U, params)
    rc = Reconstruction();
    Y = Ku;
    Y_norm = norm(Y(:));
    Y = Y / Y_norm;
    mem = struct('Y', Y, 'U', U);
    X = phi(Y,mem);
    [z,mem.B] = psi_test(X,mem);
    [Nx,Ny,Ns] = size(X);
    mem.X_dims = [Nx,Ny,Ns];
    mem.psi_dims = [Nx,Ny,numel(z)];
    Z = psi(X,mem);
    W1 = zeros(size(Z));

    prox1 = @(rc) updateX(rc, mem);
    prox2 = @(rc) updateZ(rc, mem);
    prox3 = @(rc) updateG(rc, mem);
    ascent1 = @(rc) updateW1(rc, mem);
    ascent2 = @(rc) updateW2(rc, mem);
    pres1 = @(rc) get_pres1(rc, mem);
    dres1 = @(rc) get_dres1(rc, mem);
    pres2 = @(rc) get_pres2(rc, mem);
    dres2 = @(rc) get_dres2(rc, mem);
    fX = @(rc) get_fX(rc, mem);
    fZ = @(rc) get_fZ(rc, mem);

    G = grad3(X,mem);
    W2 = zeros(size(G));
    rc.primal = {X,Z,G};
    rc.proximal = {prox1, prox2, prox3};
    rc.dual = {W1,W2};
    rc.ascent = {ascent1, ascent2};
    rc.pres = {pres1, pres2};
    rc.dres = {dres1, dres2};
    rc.obj = {fX, fZ};
    rc.rho = params.rho;
    rc.mu = params.mu;
    rc.gu = params.gu;
    rc.gl = params.gl;
    rc.lambda = params.lambda;
    rc.ptol = params.ptol;
    rc.dtol = params.dtol;
    rc.iter_max = 200;
    rc.admm(params.verbose);
    fval = get_fX(rc, mem);
    gval = get_fZ(rc, mem);
    Xrc = rc.primal{1} * Y_norm;
    mem.U = ones(size(mem.U));
    Krc = phi(Xrc, mem);
end

function out = get_fX(rc, mem)
    out = 0.5 * norm(mem.Y - phi(rc.primal{1},mem), 'fro')^2;
end

function out = get_fZ(rc, mem)
    temp = psi(rc.primal{1},mem);
    out = norm(temp(:),1);
end

function Y = phi(X,mem)
    Y = mem.U .* ifftn(ifftshift(X)) * sqrt(numel(mem.U));
end

function X = phiH(Y,mem)
    X = fftshift(fftn(mem.U .* Y)) / sqrt(numel(mem.U));
end

function Z = psi(X,mem)
    Z = zeros(mem.psi_dims);
    for nx = 1:mem.X_dims(1)
        for ny = 1:mem.X_dims(2)
            u = X(nx,ny,:);
            Z(nx,ny,:) = wavedec(u(:), 4, 'db4');
        end
    end
end

function [z,B] = psi_test(X,mem)
    u = X(1,1,:);
    [z,B] = wavedec(u(:), 4, 'db4');
end

function X = psiH(Z,mem)
    X = zeros(mem.X_dims);
    for nx = 1:mem.X_dims(1)
        for ny = 1:mem.X_dims(2)
            v = Z(nx,ny,:);
            X(nx,ny,:) = waverec(v(:), mem.B, 'db4');
        end
    end
end

function x = A_pcg(x,mem,rho1, rho2)
    temp = reshape(x,size(mem.U));
    X = phiH(phi(temp,mem),mem) + rho1 * temp - rho2*div3(grad3(temp));
    x = X(:);
end

function X = updateX(rc, mem)
    rhs = phiH(mem.Y, mem) + rc.rho(1) *psiH(rc.primal{2} - rc.dual{1}, mem) + rc.rho(2)*-div3(rc.primal{3} - rc.dual{2}, mem);
    rhs = rhs(:);
    A = @(x) A_pcg(x, mem, rc.rho(1),rc.rho(2));
    X = rc.primal{1};
    x = X(:);
    [x, flag] = pcg(A,rhs,1e-8,200,[],[],x);
    X = reshape(x, size(X));
end

function Z = updateZ(rc, mem)
    V = psi(rc.primal{1},mem) + rc.dual{1};
    Z = wthresh(V(:), 's', rc.lambda(1) / rc.rho(1));
    Z = reshape(Z, size(V));
end

function G = updateG(rc,mem)
    V = grad3(rc.primal{1},mem) + rc.dual{2};
    G = wthresh(V(:), 's', rc.lambda(2) / rc.rho(2));
    G = reshape(G,size(V));
end

function W = updateW1(rc, mem)
    W = rc.dual{1} + psi(rc.primal{1},mem) - rc.primal{2};
end

function W = updateW2(rc, mem)
    W = rc.dual{2} + grad3(rc.primal{1},mem) - rc.primal{3};
end


function [pres, prel] = get_pres1(rc,mem)
    pres = norm(rc.primal{2} - psi(rc.primal{1}, mem), 'fro');
    prel = max(norm(psi(rc.primal{1}, mem), 'fro'), norm(rc.primal{2}, 'fro'));
end

function [dres, drel] = get_dres1(rc,mem)
    dres = norm(psiH(rc.primal{2} - rc.primal_prev{2}, mem), 'fro');
    drel = norm(psiH(rc.dual{1}, mem), 'fro');
end

function [pres, prel] = get_pres2(rc,mem)
    pres = norm(rc.primal{3} - grad3(rc.primal{1}, mem), 'fro');
    prel = max(norm(grad3(rc.primal{1}, mem), 'fro'), norm(rc.primal{3}, 'fro'));
end

function [dres, drel] = get_dres2(rc,mem)
    dres = norm(-div3(rc.primal{3} - rc.primal_prev{3}, mem), 'fro');
    drel = norm(-div3(rc.dual{2}, mem), 'fro');
end


function V = grad3(X,mem)
    V = zeros(size(X,1), size(X,2), size(X,3), 3);
    V(:,:,:,1) = X([2:end 1],:,:) - X;
    V(:,:,:,2) = X(:,[2:end 1],:) - X;
    V(:,:,:,3) = X(:,:,[2:end 1],:) - X;
end

function X = div3(V,mem)
    X = zeros(size(V,1), size(V,2), size(V,3));
    X = X + V(:,:,:,1) - V([end 1:end-1],:,:,1);
    X = X + V(:,:,:,2) - V(:,[end 1:end-1],:,2);
    X = X + V(:,:,:,3) - V(:,:,[end 1:end-1],3);
end