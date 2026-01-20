function cfg = configure()
    addpath(genpath("../FID-A-master"));
    cfg = struct();
    cfg.n = 2048;
    cfg.sw = 2000;
    cfg.B0 = 3;
    cfg.metabs = { ...
        'NAA','NAAG','Cr','PCr','GPC','PCh', ...
        'Ins','Glu','Gln','GABA','GSH','Asp','Tau','PE'
    };
    data = readtable("data/metab_df.csv");
end



