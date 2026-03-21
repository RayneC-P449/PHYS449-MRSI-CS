classdef Phantom
   properties    
      labels
      metabs_list
      metab_data
      ppm_ref;
      mm_list
      mm_data
      mm_freqs
      mm_T2;
      water_data;
      gyro
      pd
   end
   methods
       function phantom = Phantom(cfg)
         if nargin > 0
            phantom.labels = cfg.labels;
            phantom.pd = cfg.pd;
            phantom.metabs_list = cfg.metabs_list;
            metabs_df = cfg.metabs;
            phantom.mm_list = cfg.mm_list;
            mm_json = cfg.mm;
            water = cfg.water;
            [I, J, K] = size(phantom.labels);
            L = numel(cfg.metabs_list);
            M = 3;
            metabs = zeros(2,L,3);
            mm = zeros(2,numel(phantom.mm_list),2);
            phantom.gyro = mm_json.B0_factor;
            for l = 1:L
                temp1 = metabs_df(metabs_df.Metabolite == phantom.metabs_list(l),:);
                for label = 1:2
                    temp2 = temp1(temp1.Label == label, :);
                    metabs(label,l,1) = fillmissing(temp2.Conc_mean,'constant',0);
                    metabs(label,l,2) = fillmissing(temp2.Conc_std,'constant',0);
                    metabs(label,l,3) = fillmissing(temp2.T2 / 1E3,'constant',Inf);
                end
            end
            phantom.mm_freqs = zeros(numel(phantom.mm_list), 1);
            for m = 1:numel(phantom.mm_list)
                for label = 1:2
                    mm(label,m,1) = mm_json.scale(m);
                end
                phantom.mm_freqs(m) = mm_json.ppm_position(m) - mm_json.ppm_offset;
                phantom.mm_T2(m) = mm_json.T2_wm(m);
            end
            phantom.metab_data = zeros(I,J,K,L,M);
            phantom.metab_data(:,:,:,:,3) = Inf;
            labelings = [2,3];
            for label = 1:2
                mask = (phantom.labels == labelings(label));
                temp = reshape(phantom.metab_data, [], L, M);
                temp(mask(:), :, :) = repmat(metabs(label,:,:), nnz(mask),1,1);
                phantom.metab_data = reshape(temp, size(phantom.metab_data));
            end
            phantom.mm_data = zeros(I,J,K,numel(phantom.mm_list), 2);
            for label = 1:2
                mask = (phantom.labels == label);
                temp = reshape(phantom.mm_data, [], numel(phantom.mm_list), 2);
                temp(mask(:), :,  :) = repmat(mm(label,:,:), nnz(mask), 1, 1);
                phantom.mm_data = reshape(temp, size(phantom.mm_data));
            end
            phantom.water_data = zeros(I,J,K,2);
            phantom.water_data(:,:,:,2) = Inf;
            for label = 1:2
                mask = (phantom.labels == label);
                temp = reshape(phantom.water_data, [], 2);
                temp(mask(:), :) = repmat(water(label, :), nnz(mask), 1);
                phantom.water_data = reshape(temp, size(phantom.water_data));
            end
            phantom.ppm_ref = mm_json.ppm_offset;
         end
      end
   end
end