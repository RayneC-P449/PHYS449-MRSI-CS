classdef Reconstruction
    properties    
      
    end
    methods
        function obj = Reconstruction()

        end
        % Sparse 2D spatial CS reconstruction from proposal
        % Implementation of ADMM based on 
        % https://www.stat.cmu.edu/~ryantibs/convexopt/lectures/admm.pdf
        % with constraint Ts = z (where x is s here and so A = T, B = -1)
        function S_cs = sparse_reconstruct_admm(~, ku, mask, T, T_adj, lambda, rho, niter)
            % TODO: write C and C_adj using wrappers for ffts that work on vectors 
            % (i.e. flattened k-space or metabolite mappings)
            C  = @(s) mask .* fft2(s);
            C_adj = @(k) ifft2(mask .* k);
            s = ifft2(ku);              
            z = T(s);
            u = zeros(size(z));
            cg_iters = 10;
            b = C_adj(ku);
            for n = 1:niter
                % Use conjugate gradient method to optimize Lagrangian
                rhs = b + rho * T_adj(z - u);
                s = pcg(@(v) normal_eq(v, C, C_adj, T, T_adj, rho), rhs, ...
                        1e-6, cg_iters, [], [], s);
                Ts = T(s);
                z = sign(Ts + u) .* max(abs(Ts + u) - lambda / rho, 0);
                u = u + rho*(Ts - z);
            end
        end
    end
end
