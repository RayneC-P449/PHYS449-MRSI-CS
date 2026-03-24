function mbasis = load_mbasis(metabs_list, seq_params, sequence)
    mbasis = cell(numel(metabs_list),1);
    for l=1:numel(metabs_list)
       metab = metabs_list(l);
       spin_sys = ['sys' char(metab)];
       load('spinSystems.mat', spin_sys);
       seq_params.sys = eval(spin_sys);
       switch sequence
           case 'press'
               basis = sim_press(seq_params.n, seq_params.sw, seq_params.B0, seq_params.lw, seq_params.sys, seq_params.tau1, seq_params.tau2);
           case 'se'
       end
       mbasis{l} = basis;
    end
end

