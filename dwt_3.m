function Krc = dwt1(Ku, U, params)
    rc = Reconstruction();
    Y = Ku;
    Y_norm = norm(Y(:));
    Y = Y / Y_norm;
    mem = struct('Y', Y, 'U', U);
    X = phi(Y,mem);
    [Nx,Ny,Ns] = size(X);
    mem.X_dims = [Nx,Ny,Ns];

    
    u = X(1,1,:);
    c = u;
    Q = X(:,:,1);
    d = wavedec2(Q,3,'db1');
    mem.psi_dims = [numel(c),numel(d)];
    [Z,mem.B] = psi(X,mem);
    
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
    rc.iter_max = 200;
    rc.admm(true);
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
    Y = mem.U .* ifftn(X) * sqrt(numel(mem.U));
end

function X = phiH(Y,mem)
    X = fftn(mem.U .* Y) / sqrt(numel(mem.U));
end

function [Z,B] = psi(X,mem)
    B = cell(2,1);
    V = zeros(size(X,1),size(X,2),mem.psi_dims(1));
    for nx = 1:mem.X_dims(1)
        for ny = 1:mem.X_dims(2)
            u = X(nx,ny,:);
            temp = u;
            V(nx,ny,:) = reshape(temp,1,1,[]);
        end
    end
    Z = zeros(mem.psi_dims);
    for n = 1:size(V,3)
        U = V(:,:,n);
        [temp,B{2}] = wavedec2(real(U),3,'db1');
        Z(n,:) = reshape(temp + 1i*wavedec2(imag(U),3,'db1'),1,[]);
    end
end

function X = psiH(Z,mem)
    V = zeros(mem.X_dims(1), mem.X_dims(2), mem.psi_dims(1));
    for n = 1:mem.psi_dims(1)
        V(:,:,n) = waverec2(real(Z(n,:)),mem.B{2},'db1');
        V(:,:,n) = V(:,:,n)+1i*waverec2(imag(Z(n,:)),mem.B{2},'db1');
    end
    X = zeros(mem.X_dims);
    for nx = 1:mem.X_dims(1)
        for ny = 1:mem.X_dims(2)
            v = V(nx,ny,:);
            X(nx,ny,:) = v;
        end
    end
end

function x = A_pcg(x,mem,rho)
    temp = reshape(x,size(mem.U));
    X = phiH(phi(temp,mem),mem) + rho * temp;
    x = X(:);
end

function X = updateX(rc, mem)
    rhs = phiH(mem.Y, mem) + rc.rho(1) * psiH(rc.primal{2} - rc.dual{1}, mem);
    rhs = rhs(:);
    A = @(x) A_pcg(x, mem, rc.rho(1));
    X = rc.primal{1};
    x = X(:);
    [x, flag] = pcg(A,rhs,1e-6,500,[],[],x);
    X = reshape(x, size(X));
end

function Z = updateZ(rc, mem)
    V = psi(rc.primal{1},mem) + rc.dual{1};
    Z = wthresh(V(:), 's', rc.lambda(1) / rc.rho(1));
    Z = reshape(Z, size(V));
end

function W = updateW(rc, mem)
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


