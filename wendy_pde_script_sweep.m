%% load data

pde_num = 10; % set to 0 to run on pre-loaded dataset
phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
tf_meth = 'FFT';
tf_param = 2;
subinds = -3;
cg = [-84 -84 -32];
eqs = ':';

nzs = [0.1 0.5 1 2];
cov_params = {[0 1],[1 0],[1 1],[2 0],[2 1]};
runs = 100;
wendy_params = {'maxits',100,'ittol',10^-4,'diag_reg',10^-inf,'trim_rows',1};

rng(1);
run_rng = randi(10^9,100,1);

stats = cell(length(nzs),length(cov_params),runs);

for ii=1:length(nzs)
    noise_ratio = nzs(ii);
    for jj=1:length(cov_params)
        for rr = 1:runs
            rng_seed = run_rng(rr);
            disp([ii jj rr])


            %%% choose PDE
            dr = '~/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/';
            pde_names = {'burgers.mat',...          %1
            'burgers_vis.mat',...               %2
            'visburgs_etdrk4.mat',...           %3
            'burgers_smallamp.mat',...          %4
            'burgers_vis0025_cubic2.mat',...    %5
            'KdV.mat',...                       %6
            'KS.mat',...                        %7
            'hyperKS.mat',...                   %8
            'lin_schrod.mat',...                %9 %%%
            'NLS.mat',...                       %10 %%% bias correct! cov helps a little
            'NLS_long.mat',...                  %11
            'transportDiff.mat',...             %12 
            'advecdiff_exact.mat',...           %13 %%%
            'AC_temp.mat',...                   %14 % - very bad, transitions to another equation?
            'fkpp_tw.mat',...                   %15 % - very bad, transitions to another equation?
            'sod.mat',...                       %16 %%
            'bw.mat',...                        %17
            'bwEpstar.mat',...                  %18
            'porous2.mat',...                   %19
            'Sine_Gordon.mat',...               %20 %%% bias correct! cov + bias
            'wave2Du3.mat',...                  %21
            'rxn_diff.mat',...                  %22
            'full_or_old/rxn_diff_old.mat',...  %23
            'Nav_Stokes.mat',...                %24
            '2D_Blast_prat_90_r.mat',...        %25
            'bwE.mat',...                       %26 %%2D
            'wave3D.mat',...                    %27
            'wave3D_N128.mat',...               %28
            '2D_Blast_prat_90_r_equi.mat',...   %29
            };
            
            if pde_num~=0
            pde_name = pde_names{pde_num};
            load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')
            end
            
            %%% get wsindy_data object
            
            %%% create data object
            Uobj = wsindy_data(U_exact,xs);
            nstates = Uobj.nstates;
            
            %%% coarsen data
            Uobj.coarsen(cg);
            
            %%% add noise
            Uobj.addnoise(noise_ratio,'seed',rng_seed);
            
            %%% get left-hand-side term
            lhsterms = lhs(eqs,:);
            numeq = size(lhsterms,1);
            
            %%% get library
            
            [lib,true_S] = true_lib(nstates,true_nz_weights(eqs));
            x_diffs = cell2mat(true_nz_weights');
            x_diffs = x_diffs(:,nstates+1:end-2);
            
            %%% get test function
            % phifun = 'pp';
            % tf_meth = 'FFT';
            % tau = 10^-10;
            % tauhat = 2;
            % tf_param = {[tau tauhat max(x_diffs(:))],[tau tauhat max(lhs(:,end))]};
            % 
            % % phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
            % % tf_meth = 'timefrac';
            % % tf_param = [0.3 0.2];
            % 
            % phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
            % tf_meth = 'FFT';
            % tf_param = 4;
            % 
            
            tf = arrayfun(@(i)testfcn(Uobj,'phifuns',phifun,'subinds',subinds,...
            'meth',tf_meth,'param',tf_param,'stateind',find(lhs(i,1:nstates),1)),(1:size(lhs,1))','uni',0);
            
            %%% build WSINDy linear system
            tic;
            WS = wendy_model(Uobj,lib,tf,cov_params{jj},'lhsterms',lhs);
            
            %%% WENDy solve
            [WS,w_its,~,~,CovW,RT] = WS_opt().wendy(WS,wendy_params{:});
            total_time_wendy = toc;
            
            %%% view governing equations, MSTLS loss, learned model accuracy
            Str_mod = WS.disp_mod;
            
            if exist('true_nz_weights','var')
                w_true = arrayfun(@(L)zeros(length(L.terms),1),WS.lib(:),'Un',0);
                if ~isempty(true_S{1})
                    for j=1:numeq;w_true{j}(true_S{j}(:,1)) = true_S{j}(:,2);end
                end
                w_true = cell2mat(w_true);
                E2 = norm(w_true-WS.weights)/norm(w_true);
                fprintf('\nCoeff err=%1.2e',E2)
            end
            disp(['runtime(s)=',num2str(total_time_wendy)])

            %%% display results
            errs = vecnorm(w_its-w_true,2,1)/norm(w_true);
            errs_OLS_WENDy = errs([1 end]);
            % disp(['runtime(s)=',num2str(total_time_wendy)])
            its_tot = size(w_its,2);
            err_all = abs(w_its(:,[1 end]) - w_true)./abs(w_true);
            res = (WS.G{1}*WS.weights-WS.b{1}+WS.bias)/norm(WS.b{1}-WS.bias);
            [~,p]=swtest(RT \ res);
            stats{ii,jj,rr} = {errs_OLS_WENDy, norm(res), p, its_tot, err_all, CovW};
        end
    end
end

%%

ind = [1 0];
cov_params = {[0 1],[1 0],[1 1],[2 0],[2 1]};

clf
clz = {'rx','bo','g.','kd','yx'};
for nz = 1:length(nzs)
for jj=1:length(cov_params)
    dat = squeeze(stats(nz,jj,:));
    dat = cellfun(@(d)d{ind(1)},dat,'uni',0);
    if ind(2)==0
        dat = cellfun(@(d) d(2)/d(1),dat);
    elseif ind(2)==1
        dat = cellfun(@(d) d(2),dat);
    end
    semilogy(nzs(nz),mean(dat),clz{jj},'markersize',10,'linewidth',3)
    xlabel('noise ratio')
    title(strrep(pde_name,'.mat',''))
    hold on
    if isequal(ind,[1 0])
        ylabel('E_2^{wendy}/E_2^{OLS}')
    elseif isequal(ind,[1 1])
        ylabel('E_2^{wendy}')
    end
    legend({'bias only','cov1 only','bias+cov1','cov2 only','bias+cov2'})
end
end
grid on
saveas(gcf,['~/Desktop/cov2_sweep_',strrep(pde_name,'.mat',''),num2str(ind),'.png'])