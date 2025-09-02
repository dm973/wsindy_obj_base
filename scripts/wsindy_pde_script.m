%% add wsindy_obj_base to path
addpath(genpath('../'))

%% load data

%%% choose PDE
dr = 'pde_data/';
pde_names = {'burgers.mat',...          
             'KS.mat',...                
             'NLS.mat',...               
             'porous2.mat',...     
             'sod_exact.mat',...
    };

pde_num = 3; % set to 0 to run on pre-loaded dataset

if pde_num~=0
    pde_name = pde_names{pde_num};
    load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')
end

%% create data object
Uobj = wsindy_data(U_exact,xs);

%%% coarsen spacetime grid
Uobj.coarsen(3);

%%% add noise
Uobj.addnoise(0.4);

%%% set library
x_diffs = [0:4];%%% differential operators
polys = [0:4]; trigs = [];%%% poly/trig functions
custom_add =  {...  %%% custom terms using term algebra
        term('fHandle',@(u,v) exp(sin(u+u.^2))),...                               % arbitrary term specified by function handle    
        compterm(term('ftag',2), diffOp([1,0],'stateind',2)),...                   % term nonlinear in a derivative
        prodterm(term('ftag',[-2i 2i]), diffOp([2,0],'stateind',1, 'nstates', 2)),...                 % product of two terms
        addterm(diffOp([3,0],'stateind',1, 'nstates', 2), term('fHandle',@(u,v) tanh(u+v))),...              % sum of two terms
    };

custom_remove_f = {}; %{@(tag) all(tag(Uobj.nstates+1:Uobj.nstates+Uobj.ndims-1))};  % remove all cross derivatives
custom_remove_t = {}; %[1 0 0 1 0 0; 0 1 0 0 1 0];                              % remove tags for divergence terms

lib = get_lib(Uobj,polys,trigs,x_diffs,custom_add,custom_remove_f,custom_remove_t);

%%% set testfcn 
phifun = 'pp';
tau = 10^-10; tauhat = 1;
tf_param = {[tau tauhat max(x_diffs)]};
tf_args = {'phifuns',phifun,'meth','FFT','param',tf_param,'subinds',-3};
tf = testfcn(Uobj,tf_args{:});

%%% scale data
Uobj.set_scales([],'lib',lib,'tf',tf);
tf = testfcn(Uobj,tf_args{:});

WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhs);

%% solve for coefficients

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

w_true = WS.reshape_w; w_true = cellfun(@(w)w*0,w_true,'un',0);
for i=1:WS.numeq
    tags = WS.lib(i).tags(:);
    ii = ~cellfun(@(tt) isnumeric(tt),tags);
    tags(ii) = repmat({zeros(1,Uobj.nstates+Uobj.ndims)},length(find(ii)),1);
    w_true{i}(ismember(cell2mat(tags),true_nz_weights{i}(:,1:end-1),'rows')) = ...
        true_nz_weights{i}(:,end);
end
cellfun(@(w,v)fprintf('coeff err=%1.2e\n',norm(w-v)/norm(v)),w_true,WS.reshape_w)
cellfun(@(w,v)fprintf('supp rec=%i\n',isequal(find(w),find(v))),w_true,WS.reshape_w)


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

%% functions

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
