classdef Reconstruction < handle
    properties    
        primal
        proximal
        dual
        ascent
        pres
        dres
        obj
        rho
        mu
        gu
        gl
        iter_max
        primal_prev
        lambda
        ptol
        dtol
    end
    methods
        function rc = Reconstruction()

        end
        
        function admm(rc, verbose)
            arguments
                rc Reconstruction
                verbose = false
            end
            for iter = 1:rc.iter_max
                rc.primal_prev = rc.primal;
                for i = 1:numel(rc.primal)
                    rc.primal{i} = rc.proximal{i}(rc);
                  
                end
                for i = 1:numel(rc.dual)
                    
                    rc.dual{i} = rc.ascent{i}(rc);
                end
               
                converged = true;
                for i = 1:numel(rc.dual)
                    res_ok = true;
                    [pres,prel] = rc.pres{i}(rc);
                    [dres,drel] = rc.dres{i}(rc);
                    pres_rel = pres / prel;
                    dres_rel = dres / drel;
                    if pres >= rc.ptol(i) * prel || dres >= rc.dtol(i) * drel
                        converged = false;
                        res_ok = false;
                    end
                    if pres_rel > rc.mu(i) * dres_rel && ~res_ok
                        rc.rho(i) = rc.rho(i) * rc.gu(i);
                    elseif dres_rel > rc.mu(i) * pres_rel && ~res_ok
                        rc.rho(i) = rc.rho(i) / rc.gl(i);
                    end
                    if verbose
                        fprintf("f = %.5e\n", rc.obj{1}(rc))
                        fprintf('pres_rel %d = %2.5f\n', i, pres_rel);
                        fprintf('dres_rel %d = %2.5f\n', i, dres_rel);
                    end
                end

              

                if converged
                    break;
                end
            end
        end
    end
end
