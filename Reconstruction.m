classdef Reconstruction
    properties    
      
    end
    methods
        function obj = Reconstruction()

        end
        
        function [X, Z, W] = reconstruct(obj, X0, Z0, W0, U, Y, rho, lambda, updateX, updateZ, updateW, max_iter, monitor)
            X = X0;
            Z = Z0;
            W = W0;
            for iter = 1:max_iter
                X = updateX(X,Z,W,U,Y,rho);
                Z = updateZ(X,Z,W,rho,lambda);
                W = updateW(X,Z,W);
                % monitor(X,Z,W,U,Y);
            end
        end
    end
end
