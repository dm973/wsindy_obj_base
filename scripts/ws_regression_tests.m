%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script tests various regression frameworks:
% OLS
% OLS with constraints
% OLS with constraints and bias
% WENDy with constraints and bias

%% boiler plate

%%% add wsindy_obj_base to path
scriptsdir = fileparts(matlab.desktop.editor.getActiveFilename);
repodir = fileparts(scriptsdir);
addpath(genpath(repodir));

%% load data

pde_num = 4; % set to 0 to run on workspace U_exact, xs, lhs variables

dr = 'pde_data/';
pde_names = {'burgers.mat',...
             'KS.mat',...           
             'NLS.mat',...
             'porous2.mat',...
             'sod_exact.mat'
    };
if pde_num~=0
    pde_name = pde_names{pde_num};
    load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')
else
    pde_name = 'custom';
end

%% create data object
Uobj = wsindy_data(U_exact,xs);

%%% coarsen spacetime grid
subsample = -54;
Uobj.coarsen(subsample);

%%% add noise
noise_ratio = 1.0;
Uobj.addnoise(noise_ratio);

%%% set testfcn 
phifun = 'pp'; tau = 1e-8; tauhat = 3; maxdiffs = 2;
tf_param = {[tau tauhat maxdiffs]};
tf_args = {'phifuns',phifun,'meth','FFT','param',tf_param,'subinds',-2};
tf = testfcn(Uobj,tf_args{:});

%%% scale data
Uobj.set_scales(1);
tf = testfcn(Uobj,tf_args{:}); % must recompute test function weights

fprintf('\n-----------%s---------------',pde_name)
fprintf('\ndata dims='); fprintf('%i ',Uobj.dims)
fprintf('\ndata scales='); fprintf('%i ',Uobj.scales)
fprintf('\ntf rads='); fprintf('%i ',tf.rads)
fprintf('\ntf powers='); fprintf('%i ',cellfun(@(p)p(end),tf.param))
fprintf('\n')

%% OLS
disp('-------------------OLS-------------------')

%%% define libary
lib = true_lib(Uobj.nstates,true_nz_weights);

%%% define wsindy_model
WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%%% get coefficients
WS = WS_opt().ols(WS);

%%% display model
print_model(WS,true_nz_weights)
Uobj.plotDyn

%% OLS with constraints
disp('-------------------OLS with constraints-------------------')

%%% define libary
lib = true_lib(Uobj.nstates,true_nz_weights);

%%% define wsindy_model
WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%%% define linear constraints
minimum_viscosity = 0.0001;
[Ai_cell,bi_cell] = stable_diffusion_constraints(WS,minimum_viscosity);
linregargs = cellfun(@(A,b) {'Aineq', A, 'bineq', b, 'verbose', 'None'}, Ai_cell, bi_cell, 'un', 0 );
linregargs_bd = WS_opt().lra_to_blkdiag(linregargs);

%%% get coefficients
WS = WS_opt().ols(WS, 'linregargs', linregargs);

%%% display model
print_model(WS,true_nz_weights)

%% OLS with constraints and bias
disp('-------------------OLS with constraints and bias correction-------------------')

%%% define libary
lib = true_lib(Uobj.nstates,true_nz_weights);

%%% define wsindy_model
WS = wendy_model(Uobj,lib,tf,[0,1],'lhsterms',lhs);

%%% define linear constraints
minimum_viscosity = -1;
[Ai_cell,bi_cell] = stable_diffusion_constraints(WS,minimum_viscosity);
linregargs = cellfun(@(A,b) {'Aineq', A, 'bineq', b, 'verbose', 'None'}, Ai_cell, bi_cell, 'un', 0 );

%%% get coefficients
WS = WS_opt().ols(WS, 'linregargs', linregargs);

%%% display model
print_model(WS,true_nz_weights)

%% WENDy with constraints and bias correct
disp('-------------------WENDy with constraints and bias correct-------------------')

%%% define libary
lib = true_lib(Uobj.nstates,true_nz_weights);

%%% define wsindy_model
WS = wendy_model(Uobj,lib,tf,[1,1],'lhsterms',lhs);

%%% define linear constraints
minimum_viscosity = -1;
[Ai_cell,bi_cell] = stable_diffusion_constraints(WS,minimum_viscosity);
linregargs = cellfun(@(A,b) {'Aineq', A, 'bineq', b, 'verbose', 'None'}, Ai_cell, bi_cell, 'un', 0 );

%%% get coefficients
[WS,w_its,res,res_0,CovW] = WS_opt().wendy(WS,'linregargs', linregargs, 'verbose',1, 'ittol', 1e-4);

%%% display model
print_model(WS,true_nz_weights)
w_true = cell2mat(cellfun(@(t)t(:,end),true_nz_weights(:),'un',0));
w_plot = WS.weights;
plot_wendy;

%% MSTLS
disp('-------------------MSTLS-------------------')

%%% define libary
x_diffs = [1];%%% differential operators
polys = [1:4]; trigs = [];%%% poly/trig functions
custom_add = {}; custom_remove_f = {@(t)max(t([1,3]))>1, @(t)max(t)>3, @(t)and(t(1),t(3)) }; custom_remove_t = [];
lib = get_lib(Uobj,polys,trigs,x_diffs, custom_add, custom_remove_f, custom_remove_t);

%%% define wsindy_model
WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%%% get coefficients
[WS,loss_wsindy] = WS_opt().MSTLS_0(WS,'lambdas', 10.^linspace(-4,0,100), 'toggle_jointthresh', 4);

% WS = WS_opt().ols(WS);

%%% display model
print_model(WS,true_nz_weights)

cond(WS.Gs{1}{1})

%%% plot MSTLS loss
if ~isempty(loss_wsindy)
    figure(2);clf;
    f = min(loss_wsindy(1,:));
    g = min(loss_wsindy(2,loss_wsindy(1,:)==f));
    for j=1:size(loss_wsindy,1)-1
        loglog(loss_wsindy(end,:),loss_wsindy(j,:),'o-',g,f,'rx')
        hold on
    end
    hold off
    legend;
end

%% MSTLS with constraints
disp('-------------------MSTLS with constraints -------------------')

%%% define libary
x_diffs = [0:2];%%% differential operators
polys = [0:4]; trigs = [];%%% poly/trig functions
custom_add = {}; custom_remove_f = {}; custom_remove_t = [];
lib = get_lib(Uobj,polys,trigs,x_diffs, custom_add, custom_remove_f, custom_remove_t);

%%% define wsindy_model
WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%%% define linear constraints
minimum_viscosity = 0.0001;
[Ai_cell,bi_cell] = stable_diffusion_constraints(WS,minimum_viscosity);
linregargs =  cellfun(@(A,b) {'Aineq', A, 'bineq', b, 'verbose', 'None'}, Ai_cell, bi_cell, 'un', 0 );

%%% get coefficients
[WS,loss_wsindy] = WS_opt().MSTLS_0(WS, 'linregargs', linregargs);

%%% display model
print_model(WS,true_nz_weights)

%%% plot MSTLS loss
if ~isempty(loss_wsindy)
    figure(2);clf;
    f = min(loss_wsindy(1,:));
    g = min(loss_wsindy(2,loss_wsindy(1,:)==f));
    for j=1:size(loss_wsindy,1)-1
        loglog(loss_wsindy(end,:),loss_wsindy(j,:),'o-',g,f,'rx')
        hold on
    end
    hold off
    legend;
end

%% MSTLS with bias correct
%%%%%%%%%%%%% currently needs to work through MSTLS_WENDy
disp('-------------------MSTLS with bias correct-------------------')

%%% define libary
x_diffs = [0:2];%%% differential operators
polys = [0:4]; trigs = [];%%% poly/trig functions
custom_add = {}; custom_remove_f = {}; custom_remove_t = [];
lib = get_lib(Uobj,polys,trigs,x_diffs, custom_add, custom_remove_f, custom_remove_t);

%%% define wsindy_model
WS = wendy_model(Uobj,lib,tf,[0,1],'lhsterms',lhs);

%%% get coefficients
WS = WS_opt().MSTLS_WENDy(WS);

%%% display model
print_model(WS,true_nz_weights)

%% MSTLS with constraints
disp('-------------------WENDy with constraints and bias correct-------------------')

%%% define libary
x_diffs = [0:2];%%% differential operators
polys = [0:4]; trigs = [];%%% poly/trig functions
custom_add = {}; custom_remove_f = {}; custom_remove_t = [];
lib = get_lib(Uobj,polys,trigs,x_diffs, custom_add, custom_remove_f, custom_remove_t);

%%% define wsindy_model
WS = wendy_model(Uobj,lib,tf,[1,1],'lhsterms',lhs);

%%% define linear constraints
minimum_viscosity = 0.0001;
[Ai_cell,bi_cell] = stable_diffusion_constraints(WS,minimum_viscosity);
linregargs ={};  cellfun(@(A,b) {'Aineq', A, 'bineq', b, 'verbose', 'None'}, Ai_cell, bi_cell, 'un', 0 );

%%% get coefficients
WS = WS_opt().MSTLS_0(WS, 'linregargs', linregargs);

%%% display model
print_model(WS,true_nz_weights)

function lib = get_lib(Uobj,polys,trigs,x_diffs,custom_add,custom_remove_f,custom_remove_t)    
    nstates = Uobj.nstates;
    ndims = Uobj.ndims;
    
    tags = get_tags(polys,trigs,nstates);
    lib = library('nstates',nstates);
    
    diff_tags = get_tags(x_diffs,[],ndims);
    diff_tags = diff_tags(diff_tags(:,end)==0,:);
    for j=1:size(tags,1)
        for i=1:size(diff_tags,1)
            if all([~and(sum(diff_tags(i,:))>0,...
                    isequal(tags(j,:),zeros(1,nstates))),...
                    ~cellfun(@(b)b([tags(j,:) diff_tags(i,:)]),custom_remove_f),...
                    ~ismember_rows([tags(j,:) diff_tags(i,:)],custom_remove_t)])
                lib.add_terms(term('ftag',tags(j,:),'linOp',diff_tags(i,:)));
            end
        end
    end
    lib.add_terms(custom_add);
end


function print_model(WS,true_nz_weights)

    Str_mod = WS.disp_mod;
    for j=1:WS.numeq
        fprintf('----------Eq %i----------\n',j)
        fprintf('%s=\n',WS.lhsterms{j}.get_str)
        cellfun(@(s)fprintf('%s\n',s),Str_mod{j})
    end
    
    if exist('true_nz_weights','var')
        w_true = inject_true_weights(WS,true_nz_weights);
        Tps = tpscore(WS.weights,w_true);
        fprintf('\nTPR=%1.2f',Tps)
        E2 = norm(w_true-WS.weights)/norm(w_true);
        fprintf('\nCoeff err=%1.2e',E2)
        fprintf('\n')
    end
end
