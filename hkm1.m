function Krc = hkm1(Ku, U, params)
    [Nx, Ny, Nt] = size(Ku);
    Xrc = zeros(Nx*Ny,Nt);
    parfor_progress(Nx*Ny); 
    parfor idx = 1:(Nx*Ny)
        [nx, ny] = ind2sub([Nx, Ny], idx);
        rc = Reconstruction();
        Y = Ku(nx,ny,:);
        Y_norm = norm(Y(:)) + 1e-16;
        Y = Y / Y_norm;
        Y = Y(:);
        U_line = U(nx,ny,:);
        U_line = U_line(:);
        mem = struct("Y", Y, "U", U_line);
        L = 40;
        mem.L = L;
        K = Nt-L+1;
        [gi, gj] = ndgrid(1:L,1:K);
        mem.g_vect = gi + gj - 1;
        mem.g_vect = mem.g_vect(:);
        id = 1:Nt;
        mem.c = zeros(Nt,1);
        mem.c(1:L) = id(1:L);
        mem.c(L+1:K) = L;
        mem.c(K+1:Nt) = L+K-id(K+1:Nt);
        X = phiH(Y, mem);
        Z = psi(X,mem);
        W = zeros(size(Z));
        prox1 = @(rc) updateX(rc, mem);
        prox2 = @(rc) updateZ(rc, mem);
        ascent1 = @(rc) updateW(rc, mem);
        pres1 = @(rc) get_pres(rc, mem);
        dres1 = @(rc) get_dres(rc, mem);
        fX = @(rc) get_fX(rc, mem);
        fZ = @(rc) get_fZ(rc, mem);
        rc.primal = {X,Z};
        rc.proximal = {prox1, prox2};
        rc.dual = {W};
        rc.ascent = {ascent1};
        rc.pres = {pres1};
        rc.dres = {dres1};
        rc.obj = {fX, fZ};
        rc.rho = params.rho;
        rc.mu = params.mu;
        rc.gu = params.gu;
        rc.gl = params.gl;
        rc.lambda = params.lambda;
        rc.ptol = params.ptol;
        rc.dtol = params.dtol;
        rc.iter_max = 1000;
        rc.admm();
        Xrc(idx,:) = rc.primal{1} * Y_norm;
        parfor_progress; 
    end
    parfor_progress(0);
    Xrc = reshape(Xrc, [Nx, Ny, Nt]);
    Krc = Xrc;
end

function Y = phi(X,mem)
    Y = mem.U .* X;
end

function X = phiH(Y,mem)
    X = mem.U .* Y;
end

function Z = psi(X,mem)
    Z = hankel(X(1:mem.L), X(mem.L:end));
end

function X = psiH(Z,mem)
    X = accumarray(mem.g_vect, Z(:));
end

function x = A_pcg(x, mem, rho)
    temp = mem.c .* x;
    x = mem.U .* x + rho * temp;
end
function X = updateX(rc,mem)
    rhs = mem.Y + rc.rho(1)*psiH(rc.primal{2} - rc.dual{1},mem);
    rhs = rhs(:);
    A = @(x) A_pcg(x, mem, rc.rho(1));
    X = rc.primal{1};
    x = X(:);
    [x, flag] = pcg(A, rhs, 1e-6, 100, [], [], x);
    X = reshape(x, size(X));
end

function Z = updateZ(rc,mem)
    [U,S,V] = svd(psi(rc.primal{1},mem) + rc.dual{1},'econ');
    S = max(S - rc.lambda(1)/rc.rho(1), 0);
    Z = U * S * V';
end

function W = updateW(rc,mem)
    W = rc.dual{1} + psi(rc.primal{1},mem) - rc.primal{2};
end

function [pres, prel] = get_pres(rc,mem)
    pres = norm(rc.primal{2} - psi(rc.primal{1}, mem), 'fro');
    prel = max(norm(psi(rc.primal{1}, mem), 'fro'), norm(rc.primal{2}, 'fro'));
end

function [dres, drel] = get_dres(rc,mem)
    dres = norm(psiH(rc.primal{2} - rc.primal_prev{2}, mem), 'fro');
    drel = norm(psiH(rc.dual{1}, mem), 'fro');
end

function out = get_fX(rc, mem)
    out = 0.5 * norm(mem.Y - phi(rc.primal{1},mem), 'fro')^2;
end

function out = get_fZ(rc, mem)
    temp = psi(rc.primal{1},mem);
    [U,S,V] = svd(temp, 'econ');
    out = sum(S(:));
end