% clear all; close all; clc
rng('shuffle'); 
rng_seed = rng().Seed; rng(rng_seed);

%% load data

%%% choose PDE
pde_num = 16; % set to 0 to run on pre-loaded dataset
dr = '~/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/';
pde_names = {'burgers.mat',...          %1
    'burgers_vis.mat',...               %2
    'KS.mat',...                        %3
    'KdV.mat',...                       %4
    'transportDiff.mat',...             %5
    'hyperKS.mat',...                   %6
    'burgers_vis_cubic.mat',...         %7
    'advecdiff_exact.mat',...           %8
    'visburgs_etdrk4.mat',...           %9
    'fkpp.mat',...                      %10
    'sod.mat',...                       %11
    'lin_schrod.mat',...                %12
    'NLS.mat',...                       %13
    'porous2.mat',...                   %14
    'Sine_Gordon.mat',...               %15
    'Nav_Stokes.mat',...                %16
    };
if pde_num~=0
    pde_name = pde_names{pde_num};
    % load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')
    load([dr,pde_name],'U_exact','lhs','xs');true_nz_weights=[];
end

%% get wsindy_data object

%%% create data object
Uobj = wsindy_data(U_exact,xs);
nstates = Uobj.nstates;

%%% coarsen data
Uobj.coarsen([-200 -200 -200]);
fprintf('\ndata dims=');fprintf('%u ',Uobj.dims);

%%% add noise
noise_ratio = 0;
rng('shuffle') % comment out to reproduce results
rng_seed = rng().Seed; rng(rng_seed); 
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%%% scale data
Uobj.set_scales(1);
scales = Uobj.scales;

%% get lhs

%%% define left-hand-side term
lhsterms = lhs;
numeq = size(lhsterms,1);

%% get library
use_true = 0; % use pre-loaded terms

%%% differential operators
x_diffs = [0:2];

%%% poly/trig functions
polys = [0:2];
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
    tau = 10^-16;
    tauhat = 1;
    tf_param = {[tau tauhat max(x_diffs)],[tau tauhat 1]};
end
tf = arrayfun(@(i)testfcn(Uobj,'phifuns',phifun,...
    'meth',tf_meth,'param',tf_param,'stateind',find(lhsterms(i,1:nstates),1)),(1:size(lhsterms,1))','uni',0);
fprintf('\ntf rads=');fprintf('%u ',tf{1}.rads);fprintf('\n')

%% build WSINDy linear system
WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhsterms);

%% MSTLS solve

lambdas = 10.^linspace(-4,0,50);
[WS,loss_wsindy,its,G,b] = WS_opt().MSTLS(WS,'lambda',lambdas);

Mscale = arrayfun(@(L)L.get_scales(Uobj.scales),WS.lib(:),'un',0);%%% account for lhs!
lhs_scales = cellfun(@(t)t.get_scale(Uobj.scales),WS.lhsterms(:),'un',0);
Mscale_W = cell2mat(cellfun(@(M,L)M/L,Mscale,lhs_scales,'un',0));
WS.add_weights(WS.weights./Mscale_W);

%% view governing equations, MSTLS loss, learned model accuracy
Str_mod = WS.disp_mod;
for j=1:numeq
    fprintf('\n----------Eq %u----------\n',j)
    cellfun(@(s)fprintf('%s\n',s),Str_mod{j})
end

figure(2);clf;
f = min(loss_wsindy(1,:));
g = min(loss_wsindy(2,loss_wsindy(1,:)==f));
for j=1:size(loss_wsindy,1)-1
    loglog(loss_wsindy(end,:),loss_wsindy(j,:),'o-',g,f,'rx')
    hold on
end
hold off

if exist('true_nz_weights','var')
    w_true = arrayfun(@(L)zeros(length(L.terms),1),WS.lib(:),'Un',0);
    if ~isempty(true_S{1})
        for j=1:numeq;w_true{j}(true_S{j}(:,1)) = true_S{j}(:,2);end
    end
    w_true = cell2mat(w_true);
    Tps = tpscore(WS.weights,w_true);
    fprintf('\nTPR=%1.2f',Tps)
    res = cellfun(@(G,w,b,M,L) ((G.*M(:)')*w-(b*L))/norm(b*L),WS.G,WS.reshape_w,WS.b,Mscale,lhs_scales,'un',0);
    res_mags = cellfun(@(r)norm(r),res);
    fprintf('\nRel. resid. =')
    fprintf('%1.2e ',res_mags)
    E2 = norm(w_true-WS.weights)/norm(w_true);
    fprintf('\nCoeff err=%1.2e',E2)
end

figure(1);clf
m = 64;
colormap([copper(m);cool(m)])
for j=1:Uobj.nstates
    subplot(Uobj.nstates,2,2*j-1)
    imagesc(Uobj.Uobs{j})
    subplot(Uobj.nstates,2,2*j)
    imagesc(U_exact{j})
end


%% true model

[lib_true,true_S] = true_lib(nstates,true_nz_weights); 
WS_true = wsindy_model(Uobj,lib_true,tf,'lhsterms',lhsterms);
WS_true.weights = cell2mat(cellfun(@(t)t(:,end),true_nz_weights(:),'uni',0));
WS_true.cat_Gb;