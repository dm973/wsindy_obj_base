%% load data

runs = 50;
nzs = [0.01 0.1 0.5 1];
tf_params = {{'meth','FFT','param',1},{'meth','FFT','param',2},{'meth','FFT','param',3},...
    {'meth','timefrac','param',0.1},{'meth','timefrac','param',0.2},{'meth','timefrac','param',0.3}};
pde_nums = [8 14 15];
coarse_param = [64 64 92 52 128 64 100 100 52 48 120 48 36 64 46 64];
subinds = -3;
eqs = ':';
cov_params = [1 1];
wendy_params = {'maxits',100,'ittol',10^-4,'diag_reg',10^-inf,'trim_rows',1};
phifun = @(v)exp(-9*[1./(1-v.^2)-1]);

for ll =1:length(pde_nums) 
    pde_num = pde_nums(ll);
    cg = -coarse_param(pde_nums(ll))*[1 1 1];
    rng(1);
    run_rng = randi(10^9,100,1);
    stats = cell(length(nzs),length(tf_params),runs);
    
    %%% choose PDE
    pde_names_script;
        
    pde_name = pde_names{pde_num};
    load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')

    for ii=1:length(nzs)
        noise_ratio = nzs(ii);
        for jj=1:length(tf_params)
            tf_param = tf_params{jj};
            parfor rr = 1:runs
                rng_seed = run_rng(rr);                
                %%% get wsindy_data object
                
                %%% create data object
                Uobj = wsindy_data(U_exact,xs);
                nstates = Uobj.nstates;
                
                %%% coarsen data
                Uobj.coarsen(cg);
                
                %%% add noise
                Uobj.addnoise(noise_ratio,'seed',rng_seed);
                
                tic;
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
                tf_param{:},'stateind',find(lhs(i,1:nstates),1)),(1:size(lhs,1))','uni',0);
                
                %%% build WSINDy linear system
                WS = wendy_model(Uobj,lib,tf,cov_params,'lhsterms',lhs);
                [WS,w_its,~,~,CovW,RT] = WS_opt().wendy(WS,wendy_params{:});
                total_time_wendy = toc;
                
                res = (WS.G{1}*WS.weights-WS.b{1}+WS.bias)/norm(WS.b{1}-WS.bias);
                [~,p]=swtest(RT \ res);
                stats{ii,jj,rr} = {w_its,CovW,norm(res),p,total_time_wendy};
                disp([ll ii jj rr])
                disp(total_time_wendy)

            end
        end
    end
    save(['~/Dropbox/Boulder/research/data/WENDy_PDE/',strrep(pde_name,'.mat',''),'_11sweep_',date,'.mat'],...
        'true_nz_weights','stats','nzs','tf_params','subinds','cg','eqs','cov_params','wendy_params','phifun')
end

%%

cap = [10^-4,10^2];
cov_params = [1 1];
pde_names_script;
for pde_num = 3:16
    try
        pde_name = pde_names{pde_num};
        dt = '07-Nov-2024';
        data_dr = '~/Dropbox/Boulder/research/data/WENDy_PDE/';
        ttl=[strrep(pde_name,'.mat',''),'_',strrep(num2str(cov_params),' ',''),'sweep_',dt];
        load([data_dr,ttl,'.mat'],...
            'true_nz_weights','stats','nzs','tf_params','subinds','cg','eqs','cov_params','wendy_params','phifun')
        true_vec = cell2mat(cellfun(@(x)x(:,end),true_nz_weights(:),'uni',0));
        errs_0 = cellfun(@(s)norm(s{1}(:,1)-true_vec)/norm(true_vec),stats);
        errs_w = cellfun(@(s)norm(s{1}(:,end)-true_vec)/norm(true_vec),stats);
        
        clf
        [x,y,r] = size(errs_0);
        for j=1:x
            subplot(ceil(x/2),2,j)
            for i=1:y 
                h1=scatter(i*ones(r,1),squeeze(errs_0(j,i,:)),'r','filled');
                hold on
                h2=scatter(i*ones(r,1),squeeze(errs_w(j,i,:)),'c');
                h3=scatter(i,mean(errs_0(j,i,:)),100,'kd','filled');
                h4=scatter(i,mean(errs_w(j,i,:)),100,'gd','filled');
            end
            set(gca,'Ylim',cap,'Yscale','log','Xtick',0:y,'Xticklabel',{'','fft1','fft2','fft3','tf.1','tf.2','tf.3'})
            grid on
            title(['noise level=',num2str(nzs(j))])
        end
        sgtitle([strrep(strrep(pde_name,'_',' '),'.mat',''),' grid=',num2str(-cg(1)),'x',num2str(-cg(2))] )
        saveas(gcf,[data_dr,ttl,'.png'])
        disp(['plot successful:',num2str(pde_num)])
    catch
        disp(['data not found:',num2str(pde_num)])
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