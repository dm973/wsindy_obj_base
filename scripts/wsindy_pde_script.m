%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script applies WSINDy to PDE data. Default variables loaded for each PDE
% in pde_names are 
% - U_exact: cell array of N (D+1)dim solution fields over (D+1)dim
%           spacetime domain with convention of time along last axis
% - xs: cell array of D+1 1D grids defining spacetime domain
% - lhs: left-hand side operator given as a vector [p1 ... pn d1 ...
%           d(D+1)] denoting the term prod_i(d/dx_i)^di prod_j(u_j^p_j)
% - (optional) true_nz_weights: cell array of K matrices, each with rows corresponding
%           to terms in the given equation, same convention as lhs, but
%           with an extra column for the term coefficient
%           *** this variable is not strictly necessary

restart_run = true;

%% add wsindy_obj_base to path

fullPathToScript = mfilename('fullpath');
currentDir = fileparts(fullPathToScript);
parentDir = fileparts(currentDir);
addpath(genpath(parentDir))

if ~restart_run
    rng('shuffle')
    clear all;
    close all; 
end

set(0,'DefaultFigureWindowStyle','docked')

%% load data

pde_num = 3; % set to 0 to run on pre-loaded dataset

%%% choose PDE
dr = 'pde_data/';
pde_names = {'burgers.mat',...          
             'KS.mat',...                
             'NLS.mat',...               
             'porous2.mat',...     
             'sod_exact.mat',...
    };
if pde_num~=0
    pde_name = pde_names{pde_num};
    load([dr,pde_name],'U_exact','xs','lhs','true_nz_weights')
else
    pde_name = 'custom';
end

%% define wsindy_data object

Uobj = wsindy_data(U_exact,xs);

%%% Subsample data
Uobj.coarsen(4);

%%% add noise
noise_ratio = 0.25;
rng_seed = rng().Seed; rng(rng_seed);
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%%% plot data
figure(1)
Uobj.plotDyn;

%% define left-hand side

lhsterms = lhs;

%% define library

x_diffs = [0:4];            %%% differential operators
polys = [0:4]; trigs = [];  %%% poly/trig functions
custom_add =  {...          %%% custom terms using term algebra
        term('fHandle',@(u,v) exp(sin(u+u.^2))),...                                  % arbitrary term specified by function handle    
        % compterm(term('ftag',2), diffOp([1,0],'stateind',2)),...                   % term nonlinear in a derivative
        % prodterm(term('ftag',[-2i 2i]), diffOp([2,0],'stateind',1, 'nstates', 2)),...                 % product of two terms
        % addterm(diffOp([3,0],'stateind',1, 'nstates', 2), term('fHandle',@(u,v) tanh(u+v))),...              % sum of two terms
    };

custom_remove_f = {}; %{@(tag) all(tag(Uobj.nstates+1:Uobj.nstates+Uobj.ndims-1))};  % remove all cross derivatives
custom_remove_t = {}; %[1 0 0 1 0 0; 0 1 0 0 1 0];                                   % remove tags for divergence terms

lib = get_lib_pde(Uobj,polys,trigs,x_diffs,custom_add,custom_remove_f,custom_remove_t);

%% define testfcn 

phifun = 'pp';
tau = 10^-10; tauhat = 1;
tf_param = {[tau tauhat max(x_diffs)]};
tf_args = {'phifuns',phifun,'meth','FFT','param',tf_param,'subinds',-3};
tf = testfcn(Uobj,tf_args{:});

%% scale data, redefine testfunction

Uobj.set_scales([],'lib',lib,'tf',tf);
tf = testfcn(Uobj,tf_args{:});

%% define WSINDy model

WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%% optimize coefficients

%%% get coefficient scale vector
Mscale = arrayfun(@(L)L.get_scales(Uobj.scales),WS.lib(:),'un',0);
lhs_scales = cellfun(@(t)t.get_scale(Uobj.scales),WS.lhsterms(:),'un',0);
Mscale = cellfun(@(M,L)M/L,Mscale,lhs_scales,'un',0);
Mscale_W = cell2mat(Mscale);

%%% optimization parameters
lambdas = 10.^linspace(-4,0,25);
threshold_scheme = 1;

[WS,loss_wsindy,its,G,b] = WS_opt().MSTLS_0(WS,'lambdas',lambdas,'M_diag',Mscale,'toggle_jointthresh',threshold_scheme,'alpha',[]);

%%% non-dimensionalized coefficients
W_nd = cellfun(@(w,m)w./m,WS.reshape_w,Mscale,'un',0); 

%% view results
clc

fprintf('\ndata dims=');fprintf('%u ',Uobj.dims);
fprintf('\ntf rads=');fprintf('%u ',WS.tf{1}{1}.rads);
fprintf('\nsize G=');fprintf('%u ',size(WS.Gs{1}{1}));fprintf('\n')

%%% display model
Str_mod = WS.disp_mod;
for j=1:WS.numeq
    fprintf('\n----------Eq %u----------\n',j)
    fprintf('%s=',WS.lhsterms{j}.get_str)
    cellfun(@(s)fprintf('%s \n',s),Str_mod{j})
end
cellfun(@(G)fprintf('cond(g)=%1.2e \n',cond(G)),WS.G)

if exist('true_nz_weights','var')
    w_true = inject_true_weights(WS,true_nz_weights);
    Tps = tpscore(WS.weights,w_true);
    fprintf('\nTPR=%1.2f',Tps)
    E2 = norm(w_true-WS.weights)/norm(w_true);
    fprintf('\nCoeff err=%1.2e',E2)
    fprintf('\nsupp rec=%i\n',isequal(find(w_true),find(WS.weights)))
else
    fprintf('no model to compare to')
end

%%% display data
figure(1);clf;
n=1;
subplot(2,1,1)
imagesc(Uobj.Uobs{n}(:,:,1))
title('observed')
subplot(2,1,2)
imagesc(U_exact{n}(:,:,1))
title('ground truth')

%%% plot MSTLS loss
figure(2);clf;
f = min(loss_wsindy(1,:));
g = min(loss_wsindy(2,loss_wsindy(1,:)==f));
for j=1:size(loss_wsindy,1)-1
    loglog(loss_wsindy(end,:),loss_wsindy(j,:),'o-',g,f,'rx')
    hold on
end
hold off
legend({'MSTLS loss','optimal lambda'});

figure(3);clf
%%% plot residual
for j=1:WS.numeq
    subplot(WS.numeq,1,j)
    plot([WS.bs{1}{j} WS.Gs{1}{j}*W_nd{j}])
    legend('b','G*w')
    title(['||G*w-b||/||b||=',num2str(norm(WS.Gs{1}{j}*W_nd{j}-WS.bs{1}{j})/norm(WS.bs{1}{j}))])
end