function bases = generate_bases(cfg)
    load spinSystems;
    bases = cell(1, length(cfg.metabs));
    for k = 1:length(cfg.metabs)
        metab = cfg.metabs{k};
        sys = eval(['sys' metab]);
        linewidth  = 2;           % Hz
        TE   = 30e-3;             % 30 ms in seconds
        tau1 = TE/4;              % PRESS definition
        tau2 = TE/4;
        basis = sim_press(cfg.n, cfg.sw, cfg.B0, linewidth, sys, tau1, tau2);
        bases{k} = basis;
    end 
end