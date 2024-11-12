%% load data
% 
% dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WENDy_data/ode_data/'; 
% odes={'wendydata_Logistic_Growth.mat',...
%     'wendydata_Lotka_Volterra.mat',...
%     'wendydata_FitzHugh-Nagumo.mat',...
%     'wendydata_Hindmarsh-Rose.mat',...
%     'wendydata_biochemM1.mat',...
%     'wendydata_alphapinene.mat',};
% filename = odes{4};
% load([dr,filename],'t','x','x0','true_prod_tags','true_vec','params','rhs_p','features');
% w_true = true_vec;

dr = '/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_ODE/'; 
odes={'cubic_longtime.mat', ...
    'lorenz_longtime.mat',...
    'lotka_longtime.mat',...
    'vanderpol_longtime.mat'};
filename = odes{1};
load([dr,filename],'t','x','weights');
w_true = cell2mat(cellfun(@(w)w(:,end),weights(:),'uni',0));
true_prod_tags = cellfun(@(w)w(:,1:end-1),weights,'uni',0);

%% get wsindy_data object

Uobj = wsindy_data(x,t);
Uobj.coarsen(-1024);

noise_ratio = 0.3;
rng('shuffle')
rng_seed = rng().Seed; rng(rng_seed);
Uobj.addnoise(noise_ratio,'seed',rng_seed,'uniform',0);

%% get lhs

eq_order = 1;
eqs = ':';
E = eye(Uobj.nstates);
lhs = arrayfun(@(i)term('ftag',E(i,:),'linOp',eq_order),(1:Uobj.nstates)','uni',0);
lhsterms = lhs(eqs);
numeq = length(lhsterms);

%% get lib tags

try
    lib = cellfun(@(tm)library('tags',tm),true_prod_tags(eqs));
catch
    lib = cellfun(@(tm)library('terms',tm),features(eqs));
end

%% get tf SVD

nstates = Uobj.nstates;
M = Uobj.dims;
tobs = Uobj.grid{1};
num_eq = length(lhs);
subinds = 1;

eta = 9;
phifuns = {@(x) exp(-eta*(1-x.^2).^(-1)), @(x) (1-x.^2).^eta};
center_scheme = 'uni';
toggle_VVp_svd = 0.999;%NaN; % default NaN. 0, no SVD reduction; in (0,1), truncates Frobenious norm; NaN, truncates SVD according to cornerpoint of cumulative sum of singular values

phifun = phifuns{1};                % defined in wendy_snf_params.m
tf_meth = 'FFT';                     % 'mtmin','FFT','direct','timefrac': default 'mtmin'
mt_params = 2.^(0:3);               % see get_rad.m:  default 2.^(0:3)

tf = cellfun(@(ls)arrayfun(@(tf_param)testfcn(Uobj,'phifuns',phifun,'subinds',subinds,...
    'meth',tf_meth,'param',tf_param,'stateind',find(ls.ftag)),mt_params),lhs,'uni',0);
mt = cell2mat(cellfun(@(t)[t.rads],tf(:),'uni',0));

K_max = 5000;
K_min = length(w_true);
mt_max = max(floor((M-1)/2)-K_min,1);

[cm,cn] = size(mt);
K = min(floor(K_max/nstates/cm), M);

V_cell = cell(num_eq,1);
Vp_cell = cell(num_eq,1);
for nn=1:num_eq
    V_cell{nn} = [];
    Vp_cell{nn} = [];
    for j=1:length(mt_params)
        V_cell{nn} = [V_cell{nn};tf{nn}(j).get_testmat(0)];
    end
    [V_cell{nn},Vp_cell{nn}] = VVp_svd(full(V_cell{nn}),K_min,tobs,toggle_VVp_svd);
end

tf = cellfun(@(V,Vp)testfcn_Vcell(Uobj,{{V,-Vp}},{[0 1]}),V_cell,Vp_cell,'uni',0);
% 
%% get test function

toggle_strong_form=0;
if toggle_strong_form==1
    phifun = 'delta';
    tf_meth = 'direct';
    tf_param = 1;
else
    phifun = 'pp';
    tf_meth = 'FFT';
    tau = 10^-10;
    tauhat = 2;
    tf_param = {[tau tauhat 1]};

    % phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    % tf_meth = 'timefrac';
    % tf_param = [0.15];

    phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    tf_meth = 'FFT';
    tf_param = 2;

    % phifun = @(v)exp(-9*[1./(1-v.^2)-1]);
    % tf_meth = 'direct';
    % tf_param = 15;

    subinds = -3;
end
tf = cellfun(@(ls)testfcn(Uobj,'phifuns',phifun,'subinds',subinds,...
    'meth',tf_meth,'param',tf_param,'stateind',find(ls.ftag)),lhs,'uni',0);

fprintf('\ntf rads=');fprintf('%u ',tf{1}.rads);fprintf('\n')

%% build WSINDy linear system
tic;
WS = wendy_model(Uobj,lib,tf,[2 1],'lhsterms',lhsterms);

%% solve
[WS,w_its,res,res_0,CovW,RT] = WS_opt().wendy(WS,'maxits',100,'ittol',10^-4,'diag_reg',10^-inf,'trim_rows',1);
total_time_wendy = toc;

%% results

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
%%% display results
errs = vecnorm(w_its-w_true,2,1)/norm(w_true);
disp(['-----------------'])
disp([' '])
disp(['rel L2 errs (OLS, WENDy)=',num2str(errs([1 end]))])
disp(['runtime(s)=',num2str(total_time_wendy)])
disp(['num its=',num2str(size(w_its,2))])
% disp(['params(true,OLS,WENDy):'])
% disp([w_true w_its(:,[1 end])])
% disp(['rel. param errors(OLS,WENDy):'])
% disp(abs(w_its(:,[1 end]) - w_true)./abs(w_true))
% disp(['-----------------'])

%% simulate learned and true reduced systems

toggle_compare = 1;
if toggle_compare==1
    w_plot = w_its(:,end);
    rhs_learned = WS.get_rhs('w',w_plot);
    tol_dd = 10^-12;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,Uobj.nstates));

    t_train = Uobj.grid{1};% linspace(tred(1),tred(end),2000);
    x0_reduced = Uobj.get_x0([]);
    [t_learned,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0_reduced,options_ode_sim);
    figure(7);clf
    for j=1:Uobj.nstates
        subplot(Uobj.nstates,1,j)
        plot(Uobj.grid{1},Uobj.Uobs{j},'b-o',t_learned,x_learned(:,j),'r-.','linewidth',2)
        try
            title(['rel err=',num2str(norm(x_learned(:,j)-Uobj.Uobs{j})/norm(Uobj.Uobs{j}))])
        catch
        end
        legend({'data','learned'})
    end
    
    if exist('w_its','var')
        try
            figure(8)
            plot_wendy;
        catch
            w_true = w_plot*0;
        end
    end
end

%%
%%%%%%%%%%%% WENDy: covariance-corrected ODE parameter estimation
%%%%%%%%%%%% Copyright 2023, All Rights Reserved
%%%%%%%%%%%% Code by Daniel Ames Messenger

function [V,Vp] = VVp_svd(V,K_min,t,toggle_VVp_svd)
    m = length(t);
    dt = mean(diff(t));
    [U,S,~] = svd(V',0);
    sings = diag(S);
    if toggle_VVp_svd>0
        s = find(cumsum(sings.^2)/sum(sings.^2)>toggle_VVp_svd^2,1);
        if isempty(s)
            s = min(K,size(V,1));
        end
    else
        corner_data = cumsum(sings)/sum(sings);
        s = getcorner_svd(corner_data,(1:length(corner_data))');%
        s = min(max(K_min,s),size(V,1));
    end
    inds = 1:s;
    V = U(:,inds)'*dt;

    Vp = V';
    Vp_hat = fft(Vp);
    if mod(m,2)==0
        k = [0:m/2 -m/2+1:-1]';
    else
        k = [0:floor(m/2) -floor(m/2):-1]';
    end
    Vp_hat = ((2*pi/m/dt)*1i*k).*Vp_hat;
    if mod(m,2)==0
        Vp_hat(m/2) = 0;        
    end
    Vp = real(ifft(Vp_hat))';

end

function tstarind = getcorner_svd(Ufft,xx)
    NN = length(Ufft);
    Ufft = Ufft/max(abs(Ufft))*NN;
    errs = zeros(NN,1);
    for k=1:NN
        [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k);
        errs(k) = sum(abs((L1-Ufft_av1)./Ufft_av1)) + sum(abs((L2-Ufft_av2)./Ufft_av2)); % relative l1           
    end
    [~,tstarind] = min(errs);
end

function [L1,L2,m1,m2,b1,b2,Ufft_av1,Ufft_av2]=build_lines(Ufft,xx,k)
   NN = length(Ufft);
   subinds1 = 1:k;
   subinds2 = k:NN;
   Ufft_av1 = Ufft(subinds1);
   Ufft_av2 = Ufft(subinds2);
   
   [m1,b1,L1]=lin_regress(Ufft_av1,xx(subinds1));
   [m2,b2,L2]=lin_regress(Ufft_av2,xx(subinds2));
end

function [m,b,L]=lin_regress(U,x)
   m = (U(end)-U(1))/(x(end)-x(1));
   b = U(1)-m*x(1);
   L = U(1)+m*(x-x(1));
end