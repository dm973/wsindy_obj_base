    dr = '~/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/';
    pde_names = {'burgers.mat',...          %1 bias=0
        'burgers_vis.mat',...               %2 bias=0
        'visburgs_etdrk4.mat',...           %3 bias=0 -- 92
        'burgers_smallamp.mat',...          %4 bias=0 -- 52
        'burgers_vis0025_cubic2.mat',...    %5 bias=0 -- 128
        'KdV.mat',...                       %6 bias=0, variable -- 64
        'KS.mat',...                        %7 bias=0, -- 100
        'hyperKS.mat',...                   %8 bias=0, -- 120
        'lin_schrod.mat',...                %9 bias=0,H=0, -- 52
        'NLS.mat',...                       %10 %%% bias correct! cov helps a little, -- 48
        'NLS_long.mat',...                  %11 -- 120
        'transportDiff.mat',...             %12 bias=0,H=0, -- 48
        'advecdiff_exact.mat',...           %13 %%% bias=0,H=0, -- 36
        'AC_temp.mat',...                   %14 % - very bad, transitions to another equation?, -- 64
        'fkpp_tw.mat',...                   %15 % - very bad, transitions to another equation?, -- 46
        'sod.mat',...                       %16 %% -- 64
        'bw.mat',...                        %17 
        'bwEpstar.mat',...                  %18 
        'porous2.mat',...                   %19 bias=0
        'Sine_Gordon.mat',...               %20 %%% bias correct! cov + bias
        'wave2Du3.mat',...                  %21 
        'rxn_diff.mat',...                  %22 
        'full_or_old/rxn_diff_old.mat',...  %23 
        'Nav_Stokes.mat',...                %24 bias=0
        '2D_Blast_prat_90_r.mat',...        %25 
        'bwE.mat',...                       %26 %%2D
        'wave3D.mat',...                    %27 bias=0,H=0
        'wave3D_N128.mat',...               %28 bias=0,H=0
        '2D_Blast_prat_90_r_equi.mat',...   %29 
        };
