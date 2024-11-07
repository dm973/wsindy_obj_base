%% load data
addpath(genpath('wsindy_obj_base'))
dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/'; 
% % filename = 'wendydata_Van_der_Pol.mat';
filename = 'wendydata_Lotka_Volterra.mat';
% % filename = 'wendydata_Hindmarsh-Rose.mat';
load([dr,filename],'t','x','x0','true_prod_tags','true_vec','params','rhs_p','features');
% tol_ode = 10^-10;
% options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
% t=linspace(0,10,10000);
% [t,x]=ode45(@(t,x)rhs_p(x,params),t,x0,options_ode_sim);
% plot(t,x)
xcell = {x};

%% get wsindy_data object
nx = ones(1,size(xcell{1},2));
xred = cellfun(@(x)x.*nx,xcell,'uni',0);
ntraj = length(xred);
tred = t;
Uobj = arrayfun(@(i)wsindy_data(xred{i},tred(:)),(1:ntraj)');
Uobj.coarsen(ceil(Uobj.dims/500));

nstates = Uobj.nstates;
M = Uobj.dims;

noise_ratio = 0;
rng('shuffle')
rng_seed = rng().Seed; rng(rng_seed);
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%% get lib tags

polys = 0:2;
trigs = [];
tags_1 = get_tags(polys,[],nstates);
% tags_1 = unique([tags_1;tags_1.*[-1 1];tags_1.*[1 -1]],'rows');
% tags_1 = tags_1(min(tags_1,[],2)>-3,:);
tags_2 = get_tags([],trigs,nstates);
tags = prodtags(tags_1,tags_2);
lib = library('tags',tags);
tags_comp = tags_1(and(max(tags_1,[],2)<3,max(tags_1,[],2)>0),:);
tags_comp = mat2cell(tags_comp,ones(size(tags_comp,1),1),size(tags_comp,2));
% lib.complib({@(x)x./(0.0001+x.^2)},tags_comp);
% lib.complib({@(x)1./(0.0001+x.^2)},{[0 1];[1 0];[1 1]});

%% get test function

phifuns = {optTFcos(2,0)};
param_tf = {'meth','param','phifuns'};
param_vals = {{'FFT'},num2cell([1]),phifuns};

s = cell(1,length(param_tf));
e = cellfun(@(p)1:length(p),param_vals,'uni',0);
[s{:}] = ndgrid(e{:});
s = cell2mat(cellfun(@(ss)ss(:),s,'uni',0));
tf = [];
for j=1:size(s,1)
    for ll=1:ntraj
        opt = {Uobj(ll)};%,'subinds',ceil(size(s,1)*1.5)};
        for k=1:size(s,2)
            opt = [opt,{param_tf{k},param_vals{k}{s(j,k)}}];
        end
        tf = [tf;{repelem({testfcn(opt{:})},1,nstates)}];
    end
end
% tf = testfcn(Uobj,'phifuns','delta','meth','direct','param',5,'mtmin',1,'subinds',2);

%% build WSINDy linear system

WS = wsindy_model(repmat(Uobj,size(s,1),1),lib,tf,'coarsen_L',1,'multitest',1);

%% solve

optm = WS_opt();
toggle_wendy = 0;

if toggle_wendy==0
    [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS);
elseif toggle_wendy==1
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',100,'regmeth','MSTLS');
elseif toggle_wendy==2
    [WS,loss_wsindy,its,G,b] = optm.MSTLS(WS,'alpha',0.01);
    [WS,w_its,res,res_0,CovW] = optm.wendy(WS,'maxits',20);
elseif toggle_wendy==3
    [WS,loss_wsindy,lambda,w_its,res,res_0,CovW] = optm.MSTLS_WENDy(WS,'maxits_wendy',50,'lambda',10.^linspace(-4,-1,50),'verbose',1);
    disp(['wendy its at optimal lambda=',num2str(size(w_its,2))])
end

%% simulate learned and true reduced systems

toggle_compare = 1:ntraj;
if ~isempty(toggle_compare)
    w_plot = WS.weights;
    rhs_learned = WS.get_rhs;
    tol_dd = 10^-12;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));
    try
        params_sparse = arrayfun(@(tt) zeros(length(tt.tags),1), WS.lib, 'uni', 0);
        ind=1;
        for i=1:length(params)
            for j=1:length(params{i})
                true_inds = gettrueinds(vertcat(WS.lib(i).tags),true_prod_tags{ind});
                params_sparse{i}(true_inds) = true_vec(ind);
                ind=ind+1;
            end
        end
        params_sparse = cell2mat(params_sparse);
        w_true = params_sparse(:);
        disp(['TPR=',num2str(tpscore(WS.weights,w_true))])
    catch
        w_true = w_plot*0;
    end

    for i=toggle_compare
        t_train = Uobj.grid{1};% linspace(tred(1),tred(end),2000);
        x0_reduced = Uobj(i).get_x0([]);
        [t_learned,xH0_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0_reduced,options_ode_sim);
        figure(7+i-1);clf
        for j=1:nstates
            subplot(nstates,1,j)
            plot(Uobj(i).grid{1},Uobj(i).Uobs{j},'b-o',t_learned,xH0_learned(:,j),'r-.','linewidth',2)
            try
                title(['rel err=',num2str(norm(xH0_learned(:,j)-Uobj(i).Uobs{j})/norm(Uobj(i).Uobs{j}))])
            catch
            end
            legend({'data','learned'})
        end
    end
    
    if exist('w_its','var')
        figure(8)
        plot_wendy;
    end
end