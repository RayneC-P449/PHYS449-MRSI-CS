function Krc = tgv(Ku, U, params)
    rc = Reconstruction();
    Y = Ku;
    Y_norm = norm(Y(:));
    Y = Y / Y_norm;
    mem = struct('Y', Y, 'U', U);
    X = phiH(Y,mem);
    V = grad3(X,mem);
    Z = grad3(X,mem) - V; % has dual
    E = eps3(V,mem); % has dual
    W1 = zeros(size(Z));
    W2 = zeros(size(E));

    rc.primal = {X,V,Z,E};
    rc.dual = {W1,W2};
    
    
    mem.alpha0 = 0.000005;
    mem.alpha1 = 0.0001;
    mem.V_size = size(V);

    rc.rho = params.rho;
    rc.lambda = params.lambda;
    rc.ptol = params.ptol;
    rc.dtol = params.dtol;
    rc.mu = params.mu;
    rc.gu = params.gu;
    rc.gl = params.gl;

    rc.primal{1} = updateX(rc,mem);
    rc.primal{2} = updateV(rc,mem);
    rc.primal{3} = updateZ(rc,mem);
    rc.primal{4} = updateE(rc,mem);
    rc.dual{1} = updateW1(rc,mem);
    rc.dual{2} = updateW2(rc,mem);

    prox1 = @(rc) updateX(rc, mem);
    prox2 = @(rc) updateV(rc, mem);
    prox3 = @(rc) updateZ(rc, mem);
    prox4 = @(rc) updateE(rc, mem);
    ascent1 = @(rc) updateW1(rc, mem);
    ascent2 = @(rc) updateW2(rc, mem);
    pres1 = @(rc) get_pres1(rc, mem);
    dres1 = @(rc) get_dres1(rc, mem);
    pres2 = @(rc) get_pres2(rc, mem);
    dres2 = @(rc) get_dres2(rc, mem);
    rc.primal = {X,V,Z,E};
    rc.proximal = {prox1, prox2, prox3, prox4};
    rc.ascent = {ascent1,ascent2};
    rc.pres = {pres1,pres2};
    rc.dres = {dres1,dres2};
    fX = @(rc) get_fX(rc,mem);
    rc.obj = {fX};

    
   

    rc.iter_max = 1000;
    disp("hey");
    rc.admm(true);

    Krc = phi(rc.primal{1},mem) * Y_norm;
end

function Y = phi(X,mem)
    Y = mem.U .* ifftn(ifftshift(X)) * sqrt(numel(mem.U));
end

function X = phiH(Y,mem)
    X = fftshift(fftn(mem.U .* Y)) / sqrt(numel(mem.U));
end

function V = grad3(X,mem)
    V = zeros(size(X,1), size(X,2), size(X,3), 2);
    V(:,:,:,1) = X([2:end 1],:,:) - X;
    V(:,:,:,2) = X(:,[2:end 1],:) - X;
end

function fval = get_fX(rc, mem)
    fval = 0.5 * norm(mem.Y - phi(rc.primal{1},mem), 'fro')^2;
end

function E = eps3(V,mem)
    E = zeros(size(V,1),size(V,2),size(V,3),2,2);
    v1 = V(:,:,:,1);
    v2 = V(:,:,:,2);
    dv1_dx = v1([2:end 1],:,:) - v1;
    dv1_dy = v1(:,[2:end 1],:) - v1;
    dv2_dx = v2([2:end 1],:,:) - v2;
    dv2_dy = v2(:,[2:end 1],:) - v2;
    E(:,:,:,1,1) = dv1_dx;
    E(:,:,:,2,2) = dv2_dy;
    E(:,:,:,1,2) = 0.5*(dv1_dy + dv2_dx);
    E(:,:,:,2,1) = E(:,:,:,1,2);
end



function W1 = updateW1(rc, mem)
    X = rc.primal{1}; V = rc.primal{2}; Z = rc.primal{3};
    W1 = rc.dual{1};
    W1 = W1 + grad3(X,mem) - V - Z;
end


function W2 = updateW2(rc, mem)
    V = rc.primal{2}; E = rc.primal{4};
    W2 = rc.dual{2};
    W2 = W2 + eps3(V) - E;
end

function x = A_pcg(x, mem, rho1)
    x = reshape(x,size(mem.U));
    x = phiH(phi(x,mem),mem) + rho1 * -div3(grad3(x));
    x = x(:);
end

function X = updateX(rc, mem)
    rho1 = rc.rho(1);
    X = rc.primal{1}; Y = mem.Y; b = rc.primal{2} + rc.primal{3} - rc.dual{1};
    rhs = phiH(Y,mem) + rho1 * -div3(b);
    rhs = rhs(:);
    A = @(x) A_pcg(x,mem,rho1);
    [x,~] = pcg(A,rhs,1e-5,200,[],[],X(:));
    X = reshape(x,size(X));
end

function v = B_pcg(v,mem,rho1,rho2)
    v = reshape(v, mem.V_size);
    v = rho1*v + rho2 * -ediv3(eps3(v));
    v = v(:);
end

function V = updateV(rc, mem)
    X = rc.primal{1}; V = rc.primal{2}; Z = rc.primal{3}; E = rc.primal{4};
    rho1 = rc.rho(1); rho2 = rc.rho(2);
    W1 = rc.dual{1}; W2 = rc.dual{2};
    Gx = grad3(X,mem);
    rhs = rho1*(Gx - Z + W1) + rho2*-ediv3(E-W2);
    rhs = rhs(:);
    B = @(v) B_pcg(v,mem,rho1,rho2);
    [v,~] = pcg(B,rhs,1e-5,200,[],[],V(:));
    V = reshape(v, size(V));
end


function X = div3(V)
    X = zeros(size(V,1), size(V,2), size(V,3));
    X = X + V(:,:,:,1) - V([end 1:end-1],:,:,1);
    X = X + V(:,:,:,2) - V(:,[end 1:end-1],:,2);
end



function V = ediv3(E)
    V = zeros(size(E,1),size(E,2),size(E,3),2);
    E11 = E(:,:,:,1,1);
    E22 = E(:,:,:,2,2);
    E12 = E(:,:,:,1,2);
    E21 = E(:,:,:,2,1);
    V(:,:,:,1) = (E11 - E11([end 1:end-1],:,:)) + 0.5*(E12 - E12(:,[end 1:end-1],:)) + 0.5*(E21 - E21([end 1:end-1],:,:));
    V(:,:,:,2) = 0.5*(E12 - E12([end 1:end-1],:,:)) + 0.5*(E21 - E21(:,[end 1:end-1],:)) + (E22 - E22(:,[end 1:end-1],:));
end

function Z = updateZ(rc,mem)
    X = rc.primal{1}; V = rc.primal{2}; W1 = rc.dual{1}; rho1 = rc.rho(1);
    temp = grad3(X) - V + W1;
    Z = wthresh(temp(:),'s', rc.lambda(1) * mem.alpha1 / rho1);
    Z = reshape(Z, size(W1));
end


function E = updateE(rc,mem)
    V = rc.primal{2}; W2 = rc.dual{2}; rho2 = rc.rho(2);
    temp = eps3(V) + W2;
    E = wthresh(temp(:), 's', rc.lambda(1) * mem.alpha0 / rho2);
    E = reshape(E, size(W2));
end

function [pres, prel] = get_pres1(rc,mem)
    G = grad3(rc.primal{1}); V = rc.primal{2}; Z = rc.primal{3};
    temp = G - V - Z;
    pres = norm(temp(:));
    prel = max([norm(G(:)),norm(V(:)),norm(Z(:))]) + 1e-9;
end

function [dres, drel] = get_dres1(rc,mem)
    dres = -div3(rc.primal{3} - rc.primal_prev{3});
    dres = norm(dres, 'fro');
    drel = -div3(rc.dual{1});
    drel = norm(drel, 'fro') + 1e-9;
end

function [pres, prel] = get_pres2(rc,mem)
    V = rc.primal{2}; E = rc.primal{4};
    epsV = eps3(V);
    temp = eps3(V) - E;
    pres = norm(temp(:));
    prel = max([norm(epsV(:)),norm(E(:))]) + 1e-9;
end

function [dres, drel] = get_dres2(rc,mem)
    dres = -ediv3(rc.primal{4} - rc.primal_prev{4});
    dres = norm(dres, 'fro');
    drel = -ediv3(rc.dual{2});
    drel = norm(drel, 'fro') + 1e-9;
end