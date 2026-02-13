classdef TestSuite
    properties

    end

    methods
        function res = t1(obj)
            cfg = struct();
            cfg.labels = niftiread('data/skeleton/labels.nii');
            yslice = 50:59;
            cfg.labels = cfg.labels(60, yslice, 131);

            cfg.metabs_list = ["NAA","NAAG","Cr","PCr"]; % short version for fast runs
            % cfg.metabs_list = ["NAA","NAAG","Cr","PCr","GPC","PCh", ...
            %         "Glu","Gln","GABA","GSH","Asp","Tau","PE", "Lac"];
            cfg.metabs = readtable('data/metabolites/metab_df.csv');
            phantom = Phantom(cfg);
            sm = SignalModel();
            seq_func = @sim_press;
            seq_params = struct();
            seq_params.n = 1024;
            seq_params.sw = 2000;
            seq_params.B0 = 3;
            seq_params.lw = 1;
            seq_params.sys = 0; 
            TE   = 30e-3;  
            seq_params.tau1 = TE/4;             
            seq_params.tau2 = TE/4;
            zslice = 1:1;
            voxrange = [1,10];
            [kspace, t, imspace, ppm] = sm.extract_kspace(phantom,@sim_press,seq_params, zslice, voxrange);
            test_plot = abs(squeeze(imspace(1,1,:))); 
            % plot(ppm, test_plot);
        end
        function res = t2(obj)
            cfg = struct();
            cfg.metabs = readtable('data/metabolites/metab_df.csv');
            cfg.metabs_list = ["NAA","NAAG","Cr","PCr"]; % short version for fast runs
            cfg.labels = zeros(2,2,1);
            cfg.labels(:,:,1) = [2 1; 0 1];
            phantom = Phantom(cfg);
            sm = SignalModel();
            seq_func = @sim_press;
            seq_params = struct();
            seq_params.n = 1024;
            seq_params.sw = 2000;
            seq_params.B0 = 3;
            seq_params.lw = 1;
            seq_params.sys = 0; 
            TE   = 30e-3;  
            seq_params.tau1 = TE/4;             
            seq_params.tau2 = TE/4;
            zslice = 1:1;
            voxrange = [2,2];
            [kspace, t, imspace, ppm] = sm.extract_kspace(phantom,@sim_press,seq_params, zslice, voxrange);
            test_plot = abs(squeeze(imspace(2,1,:))); 
            % plot(ppm, test_plot);
        end

        

    end

end