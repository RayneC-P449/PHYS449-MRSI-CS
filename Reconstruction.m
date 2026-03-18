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
        function obj = Reconstruction()

        end
        
        function admm(rc, verbose)
            function say(msg)
                if verbose
                    disp(msg);
                end
            end
            for iter = 1:rc.iter_max
                rc.primal_prev = rc.primal;
                for i = 1:numel(rc.primal)
                    rc.primal{i} = rc.proximal{i}(rc);
                end
                for i = 1:numel(rc.dual)
                    rc.dual{i} = rc.ascent{i}(rc);
                end
                say(join(["iter = ", iter]));
                converged = true;
                for i = 1:numel(rc.dual)
                    [pres,prel] = rc.pres{i}(rc);
                    [dres,drel] = rc.dres{i}(rc);
                    pres_rel = pres / prel;
                    dres_rel = dres / drel;
                    say(join(["dual_id = ", i]));
                    say(join(["rho = ", rc.rho(i), "lambda = ", rc.lambda(i)]))
                    say(join(["pres_rel = ", pres_rel]));
                    say(join(["dres_rel = ", dres_rel]));
                    if pres >= rc.ptol(i) * prel || dres >= rc.dtol(i) * drel
                        converged = false;
                    end
                    if pres_rel > rc.mu(i) * dres_rel
                        rc.rho(i) = rc.rho(i) * rc.gu(i);
                        % rc.lambda(i) = rc.lambda(i) * rc.gu(i);
                    elseif dres_rel > rc.mu(i) * pres_rel
                        rc.rho(i) = rc.rho(i) / rc.gl(i);
                        % rc.lambda(i) = rc.lambda(i) / rc.gl(i);
                    end
                end
                for i = 1:numel(rc.primal)
                    say(join(["primal_id = ", i]));
                    say(join(["obj_val = ", rc.obj{i}(rc)]));
                end
                say("-------");
                if converged
                    % disp("Done");
                    break;
                end
            end
        end

        
    end
end
