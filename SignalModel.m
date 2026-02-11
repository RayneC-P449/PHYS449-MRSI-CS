classdef SignalModel
   properties    
      
   end
   methods
       function sm = SignalModel()

       end
       function [kspace, t, imspace, ppm] = extract_kspace(obj, phantom, seq_func, seq_params, zslice, voxrange)
           %   INPUTS:
           %       phantom - the brain phantom model
           %       seq_func - MRS function to simulate (e.g. PRESS)
           %       seq_params - arguments for seq_func
           %       zslice - range of z-indices of phantom to sum over for 
           %            2D MRSI signal over
           %       voxrange - vector indicating number of MRSI voxels in 
           %            x and y direction. Typically smaller
           %            than the phantom's spatial dimensions
           %   OUTPUTS:
           %       TBA
           kspace_dims = [voxrange, seq_params.n];
           kspace = zeros(kspace_dims);
           bases = {};
           load spinSystems
           % Acquire basis for each metabolite system in
           % phanton.metabs_list
           for l=1:numel(phantom.metabs_list)
               metab = phantom.metabs_list(l);
               seq_params.sys = eval(['sys' char(metab)]);
               args = struct2cell(seq_params);
               basis = seq_func(args{:});
               basis.fids = basis.fids;
               bases{l} = basis;
           end
           [Nx, Ny, Nz] = size(phantom.labels);
           x = linspace(-0.5, 0.5, Nx+1); kx = -voxrange(1)/2:1:voxrange(1)/2-1;
           x = x(1:end-1);
           y = linspace(-0.5, 0.5, Ny+1); ky = -voxrange(2)/2:1:voxrange(2)/2-1;
           y = y(1:end-1);
           % Simulate phase encoding steps with voxrange(1) steps in 
           % x-direction and voxrange(2) steps in y-direction
           for vx = 1:voxrange(1)
               for vy = 1:voxrange(2)
                   phase_x = squeeze(kx(vx) * x);
                   phase_y = squeeze(ky(vy) * y);
                   for l = 1:numel(phantom.metabs_list)
                       basis_fid = reshape(bases{l}.fids, 1, 1, []);
                       for sc = zslice
                           temp = phantom.data(:,:,sc,l,1);
                           temp = temp .* exp(1i*-2*pi*phase_x.');
                           temp = temp .* exp(1i*-2*pi*phase_y);
                           % TODO: incorporate variation in basis linewidth 
                           % based on randomness + tissue type     
                           s = temp .* basis_fid;
                           s = sum(sum(s,1),2);
                           s = squeeze(s);
                           s = reshape(s, 1, 1, []);
                           kspace(vx, vy, :) = kspace(vx, vy, :) + s;
                       end
                   end
               end
           end
           t = bases{1}.t;
           ppm = bases{1}.ppm;
           imspace = fftshift(ifftn(kspace));
       end
   end
end