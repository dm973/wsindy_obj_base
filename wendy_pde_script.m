%% load data

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
    'NLS.mat',...                       %10
    'NLS_long.mat',...                  %11
    'transportDiff.mat',...             %12 
    'advecdiff_exact.mat',...           %13 %%%
    'AC_temp.mat',...                   %14
    'fkpp.mat',...                      %15
    'sod.mat',...                       %16 %%
    'bw.mat',...                        %17
    'bwEpstar.mat',...                  %18
    'porous2.mat',...                   %19
    'Sine_Gordon.mat',...               %20
    'wave2Du3.mat',...                  %21
    'rxn_diff.mat',...                  %22
    'full_or_old/rxn_diff_old.mat',...  %23
    'Nav_Stokes.mat',...                %24
    '2D_Blast_prat_90_r.mat',...        %25
    'bwE.mat',...                       %26 %%2D
    'wave3D.mat',...                    %27
    'wave3D_N128.mat',...               %28
    };

pde_num = 20; % set to 0 to run on pre-loaded dataset

if pde_num~=0
    pde_name = pde_names{pde_num};
    load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')
end

%% get wsindy_data object

%%% create data object
Uobj = wsindy_data(U_exact,xs);
nstates = Uobj.nstates;

%%% coarsen data
Uobj.coarsen([-64 -64 -100]);
fprintf('\ndata dims=');fprintf('%u ',Uobj.dims);fprintf('\n')

%%% add noise
noise_ratio = 0.01;
rng('shuffle') % comment out to reproduce results
rng_seed = rng().Seed; rng(rng_seed); 
Uobj.addnoise(noise_ratio,'seed',rng_seed);

tic;
%% get left-hand-side term

lhsterms = lhs;
numeq = size(lhsterms,1);

%% get library

[lib,true_S] = true_lib(nstates,true_nz_weights);
x_diffs = cell2mat(true_nz_weights');
x_diffs = x_diffs(:,nstates+1:end-2);

%% get test function

toggle_strong_form=0;
if toggle_strong_form==1
    phifun = 'delta';
    tf_meth = 'direct';
    tf_param = max(x_diffs);
else
    phifun = 'pp';
    tf_meth = 'FFT';
    tau = 10^-8;
    tauhat = 2;
    tf_param = {[tau tauhat max(x_diffs(:))],[tau tauhat max(lhs(:,end))]};

    phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    tf_meth = 'timefrac';
    tf_param = [0.15 0.15];

    subinds = -3;
end
tf = arrayfun(@(i)testfcn(Uobj,'phifuns',phifun,'subinds',subinds,...
    'meth',tf_meth,'param',tf_param,'stateind',find(lhs(i,1:nstates),1)),(1:size(lhs,1))','uni',0);
fprintf('\ntf rads=');fprintf('%u ',tf{1}.rads);fprintf('\n')

%% build WSINDy linear system

WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%% WENDy solve

[WS,w_its,res,res_0,CovW] = WS_opt().wendy(WS,'maxits',20,'ittol',10^-4,'diag_reg',10^-4,'trim_rows',1);
total_time_wendy = toc;

%% view governing equations, MSTLS loss, learned model accuracy
Str_mod = WS.disp_mod;
for j=1:numeq
    fprintf('\n----------Eq %u----------\n',j)
    cellfun(@(s)fprintf('%s\n',s),Str_mod{j})
end

if exist('true_nz_weights','var')
    w_true = arrayfun(@(L)zeros(length(L.terms),1),WS.lib(:),'Un',0);
    if ~isempty(true_S{1})
        for j=1:numeq;w_true{j}(true_S{j}(:,1)) = true_S{j}(:,2);end
    end
    w_true = cell2mat(w_true);
    fprintf('\nRel. resid. =')
    fprintf('%1.2e ',norm(WS.res))
    E2 = norm(w_true-WS.weights)/norm(w_true);
    fprintf('\nCoeff err=%1.2e',E2)
end

% figure(1);clf
m = 64;
colormap([copper(m);cool(m)])
for j=1:Uobj.nstates
    if Uobj.ndims==2
        subplot(Uobj.nstates,2,2*j-1)
        imagesc(Uobj.Uobs{j})
        subplot(Uobj.nstates,2,2*j)
        imagesc(U_exact{j})
    elseif Uobj.ndims==3
        for k=1:Uobj.dims(3)
            subplot(Uobj.nstates,2,2*j-1)
            imagesc(Uobj.Uobs{j}(:,:,k))
            subplot(Uobj.nstates,2,2*j)
            imagesc(U_exact{j}(:,:,k))
            drawnow
        end
    end
end

%%% display results
errs = vecnorm(w_its-w_true,2,1)/norm(w_true);
disp(['-----------------'])
disp(['------',strrep(pde_name,'.mat',''),'-------'])
disp([' '])
disp(['rel L2 errs (OLS, WENDy)=',num2str(errs([1 end]))])
disp(['runtime(s)=',num2str(total_time_wendy)])
disp(['num its=',num2str(size(w_its,2))])
disp(['params(true,OLS,WENDy):'])
disp([w_true w_its(:,[1 end])])
disp(['rel. param errors(OLS,WENDy):'])
disp(abs(w_its(:,[1 end]) - w_true)./abs(w_true))
disp(['-----------------'])