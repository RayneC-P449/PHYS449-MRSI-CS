classdef SignalModel
   properties    
      
   end
   methods
       function sm = SignalModel()

       end
       function [kspace, t, imspace, ppm] = extract_kspace(obj, phantom, seq_func, seq_params, params_order, zslice, voxrange, sigma, beta)
           kspace_dims = [voxrange, seq_params.n];
           kspace = zeros(kspace_dims);
           metab_bases = {};
           load spinSystems
           for l=1:numel(phantom.metabs_list)
               metab = phantom.metabs_list(l);
               seq_params.sys = eval(['sys' char(metab)]);
               args = struct2cell(orderfields(seq_params, params_order));
               basis = seq_func(args{:});
               basis.fids = basis.fids;
               metab_bases{l} = basis;
           end
           t_vect = metab_bases{1}.t;
           for m=1:numel(phantom.mm_list)
               f_m = seq_params.B0 * phantom.gyro * phantom.mm_freqs(m);
               T_m = phantom.mm_T2(m) / 10E3;
               mm_bases{m} = exp(-t_vect/T_m).*exp(2*pi*1i*f_m*t_vect);
           end
           T2_prime = 1/(pi*seq_params.lw);
           water_bases = {};
           for w = 1:5
               f_water = 4.4 + 0.5 * rand(1);
               f_water = (phantom.ppm_ref - f_water) .* seq_params.B0 * phantom.gyro;
               water_bases{w} = beta * exp(2*pi*1i*f_water*t_vect) .* exp(-t_vect / T2_prime);
           end
           [Nx, Ny, Nz] = size(phantom.labels);
           x = linspace(-0.5, 0.5, Nx+1); kx = -voxrange(1)/2:1:voxrange(1)/2-1;
           x = x(1:end-1);
           y = linspace(-0.5, 0.5, Ny+1); ky = -voxrange(2)/2:1:voxrange(2)/2-1;
           y = y(1:end-1);
           bases_fid = zeros(Nx,Ny,seq_params.n);
           for sc = zslice
               for l = 1:numel(phantom.metabs_list)
                   metab_fid = reshape(metab_bases{l}.fids, 1, 1, []);
                   bases_fid = bases_fid + (phantom.metab_data(:,:,sc,l,1) .* metab_fid) .* exp(reshape(-t_vect,1,1,[]) ./ phantom.metab_data(:,:,sc,l,3));
               end
               for m = 1:numel(phantom.mm_list)
                   mm_fid = reshape(mm_bases{m}, 1, 1, []);
                   bases_fid = bases_fid + phantom.mm_data(:,:,sc,m,1) .* mm_fid;
               end
               Aw = rand(1,5);
               Aw = Aw / sum(Aw);
               for w = 1:5
                   water_fid = reshape(water_bases{w}, 1, 1, []);
                   bases_fid = bases_fid + Aw(w) * (phantom.water_data(:,:,sc,1) .* water_fid) .* exp(reshape(-t_vect,1,1,[]) ./ phantom.water_data(:,:,sc,2));
               end
           end
           for vx = 1:voxrange(1)
               for vy = 1:voxrange(2)
                   phase_x = squeeze(kx(vx) * x);
                   phase_y = squeeze(ky(vy) * y);
                   temp = bases_fid;
                   temp = temp .* exp(1i*2*pi*phase_x.');
                   temp = temp .* exp(1i*2*pi*phase_y);
                   temp = sum(sum(temp,1),2);
                   kspace(vx,vy,:) = temp;         
                end
           end
           noise = sigma*(randn(kspace_dims) + 1i*randn(kspace_dims));
           kspace = kspace + noise;
           t = metab_bases{1}.t;
           ppm = metab_bases{1}.ppm;
           imspace = fftshift(fftn(kspace) / numel(kspace));
       end
   end
end