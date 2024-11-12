%% load data

%%% choose PDE
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

pde_num = 7; % set to 0 to run on pre-loaded dataset

if pde_num~=0
    pde_name = pde_names{pde_num};
    load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')
end

%% get wsindy_data object

%%% create data object
Uobj = wsindy_data(U_exact,xs);
nstates = Uobj.nstates;

%%% coarsen data
Uobj.coarsen([-128*[1 1] -48]);

%%% add noise
noise_ratio = 0.2;
rng('shuffle') % comment out to reproduce results
rng_seed = rng().Seed; rng(rng_seed); 
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%% get left-hand-side term

eqs = ':';
lhsterms = lhs(eqs,:);
numeq = size(lhsterms,1);

%% get library

[lib_true,true_S] = true_lib(nstates,true_nz_weights(eqs));

%% get library
use_true = 0; % use pre-loaded terms

%%% differential operators
x_diffs = [0:4];

%%% poly/trig functions
polys = [0:4];
trigs = [];

%%% custom terms
custom_terms = {};

if ~use_true==1
    [lib,true_S] = trig_poly_lib(polys,trigs,x_diffs,nstates,Uobj.ndims,numeq,true_nz_weights);
    lib = repelem(lib,1,numeq);
    if ~isempty(custom_terms); for j=1:numeq; lib(j).add_terms(custom_terms{j}); end;end
else
    [lib,true_S] = true_lib(nstates,true_nz_weights);
end

%% get test function

toggle_strong_form=0;
if toggle_strong_form==1
    phifun = 'delta';
    tf_meth = 'direct';
    tf_param = max(x_diffs);
else
    phifun = 'pp';
    tf_meth = 'FFT';
    tau = 10^-10;
    tauhat = 2;
    tf_param = {[tau tauhat max(x_diffs(:))],[tau tauhat max(lhs(:,end))]};

    % phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    % tf_meth = 'timefrac';
    % tf_param = [0.3 0.2];

    phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    tf_meth = 'FFT';
    tf_param = 3;

    % phifun = @(v)exp(-9*[1./(1-v.^2)aa-1]);
    % tf_meth = 'direct';
    % tf_param = [10 13];

    % phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    % tf_meth = 'timefrac';
    % tf_param = 0.1;

    subinds = -4;
end
tf = arrayfun(@(i)testfcn(Uobj,'phifuns',phifun,'subinds',subinds,...
    'meth',tf_meth,'param',tf_param,'stateind',find(lhs(i,1:nstates),1)),(1:size(lhs,1))','uni',0);

%% build WSINDy linear system
tic;
% WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs,'statcorrect',[0 1]);
WS = wendy_model(Uobj,lib,tf,[1 1],'lhsterms',lhs);

%% WENDy solve

toggle_wendy = 1;

MSTLS_args = {'lambda',10.^linspace(-4,0,100)};
WENDy_args = {'maxits_wendy',2,'diag_reg',10^-inf,'verbose',1};
MSTLS_WENDy_args = {'cov_thresh',0.5};

if toggle_wendy==0
    [WS,loss_wsindy,its,G,b] =  WS_opt().MSTLS(WS,MSTLS_args{:});
    w_its = WS.weights;
    res = WS.G{1}*WS.weights-WS.b{1};    
    res_0 = res;
    CovW = rms(res)^2*inv(WS.G{1}(:,WS.weights~=0)'*WS.G{1}(:,WS.weights~=0));
elseif toggle_wendy==1
    [WS,loss_wsindy,lambda,w_its,res,res_0,CovW,RT] =  WS_opt().MSTLS_WENDy(WS,WENDy_args{:},MSTLS_args{:},MSTLS_WENDy_args{:});
    disp(['wendy its at optimal lambda=',num2str(size(w_its,2))])
end

total_time_wendy = toc;

%% view governing equations, MSTLS loss, learned model accuracy
Str_mod = WS.disp_mod;
for j=1:numeq
    fprintf('\n----------Eq %u----------\n',j)
    cellfun(@(s)fprintf('%s\n',s),Str_mod{j})
end

fprintf('\ndata dims=');fprintf('%u ',Uobj.dims);fprintf('\n')
fprintf('\ntf rads=');fprintf('%u ',tf{1}.rads);fprintf('\n')

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

figure(1);
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
            imagesc(U_exact{j}(:,:,k*floor(length(xs{end})/Uobj.dims(end))))
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
% disp(['params(true,OLS,WENDy):'])
% disp([w_true w_its(:,[1 end])])
% disp(['rel. param errors(OLS,WENDy):'])
% disp(abs(w_its(:,[1 end]) - w_true)./abs(w_true))
% disp(['-----------------'])

%% plot wendy results

figure(2);
w_plot = w_its(:,end);
plot_wendy;

%% plot loss

figure(3);clf;
f = min(loss_wsindy(1,:));
g = min(loss_wsindy(2,loss_wsindy(1,:)==f));
for j=1:size(loss_wsindy,1)-1
    loglog(loss_wsindy(end,:),loss_wsindy(j,:),'o-',g,f,'rx')
    hold on
end
hold off