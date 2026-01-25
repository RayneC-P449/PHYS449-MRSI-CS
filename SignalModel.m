classdef SignalModel
   properties    
      
   end
   methods
       function sm = SignalModel()

       end
       function signal = extract_signal(obj, phantom, seq_func, seq_params, voxrange)
           sig_dims = [voxrange, seq_params.n];
           signal = zeros(sig_dims);
           load spinSystems
           for l=1:numel(phantom.metabs_list)
               metab = phantom.metabs_list(l);
               seq_params.sys = eval(['sys' char(metab)]);
               args = struct2cell(seq_params);
               basis = seq_func(args{:});
           end
           for dim = 1:3
               N = vorange(dim);
               

           end
           



           
       end
   end
end