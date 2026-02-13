classdef Phantom < handle
   properties    
      labels
      metabs_list
      data
   end
   methods
       % Function for constructing the digital phantom based on
       % concentration and T2 data from 
       % https://github.com/dennisvds/MRS-Digital-Phantom/tree/main/data/metabolites
       % and MRILab's tissue label map.
       function phantom = Phantom(cfg)
         if nargin > 0
            phantom.labels = cfg.labels;
            phantom.metabs_list = cfg.metabs_list;
            metabs_df = cfg.metabs;
            [I, J, K] = size(phantom.labels);
            L = numel(cfg.metabs_list);
            M = 3;
            metabs = zeros(2,L,3);
            for l = 1:L
                temp1 = metabs_df(metabs_df.Metabolite == phantom.metabs_list(l),:);
                for label = 1:2
                    temp2 = temp1(temp1.Label == label, :);
                    metabs(label,l,1) = fillmissing(temp2.Conc_mean,'constant',0);
                    metabs(label,l,2) = fillmissing(temp2.Conc_std,'constant',0);
                    metabs(label,l,3) = fillmissing(temp2.T2,'constant',0);
                end
            end
            phantom.data = zeros(I,J,K,L,M);
            labelings = [3,2]; % To make the dataframe data and the MRILab labelings match
            for label = 1:2
                mask = (phantom.labels == labelings(label));
                temp = reshape(phantom.data, [], L, M);
                temp(mask(:), :, :) = repmat(metabs(label,:,:), nnz(mask),1,1);
                phantom.data = reshape(temp, size(phantom.data));
            end
         end
      end
   end
end