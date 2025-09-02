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

%% get wsindy_data object

%%% create data object
Uobj = wsindy_data(U_exact,xs);
nstates = Uobj.nstates;

%%% coarsen data
Uobj.coarsen(3);
fprintf('\ndata dims=');fprintf('%u ',Uobj.dims);fprintf('\n')

%%% add noise
noise_ratio = 0.4;
rng('shuffle') % comment out to reproduce results
rng_seed = rng().Seed; rng(rng_seed); 
Uobj.addnoise(noise_ratio,'seed',rng_seed);

%% get library

[lib,true_S] = true_lib(nstates,true_nz_weights);
x_diffs = cell2mat(true_nz_weights');
x_diffs = x_diffs(:,nstates+1:end-2);

%% get test function

toggle_strong_form = 0;
if toggle_strong_form==1
    phifun = 'delta';
    tf_meth = 'direct';
    tf_param = max(x_diffs);
else
    phifun = 'pp';
    tf_meth = 'FFT';
    tau = 10^-10;
    tauhat = 1;
    tf_param = {[tau tauhat max(x_diffs(:))],[tau tauhat max(lhs(:,end))]};

    subinds = -3;
end
tf = arrayfun(@(i)testfcn(Uobj,'phifuns',phifun,'subinds',subinds,...
    'meth',tf_meth,'param',tf_param,'stateind',find(lhs(i,1:nstates),1)),(1:size(lhs,1))','uni',0);
fprintf('\ntf rads=');fprintf('%u ',tf{1}.rads);fprintf('\n')

%% build WSINDy linear system
tic;
WS = wendy_model(Uobj,lib,tf,[2 1],'lhsterms',lhs);

%% WENDy solve
[WS,w_its,res,res_0,CovW,RT] = WS_opt().wendy(WS,'maxits',100,'ittol',10^-4,'diag_reg',10^-inf,'trim_rows',1);
total_time_wendy = toc;

%% view governing equations, MSTLS loss, learned model accuracy
Str_mod = WS.disp_mod;
for j=1:WS.numeq
    fprintf('\n----------Eq %u----------\n',j)
    cellfun(@(s)fprintf('%s\n',s),Str_mod{j})
end

if exist('true_nz_weights','var')
    w_true = arrayfun(@(L)zeros(length(L.terms),1),WS.lib(:),'Un',0);
    if ~isempty(true_S{1})
        for j=1:WS.numeq;w_true{j}(true_S{j}(:,1)) = true_S{j}(:,2);end
    end
    w_true = cell2mat(w_true);
    fprintf('\nRel. resid. =')
    fprintf('%1.2e ',norm(WS.res))
    E2 = norm(w_true-WS.weights)/norm(w_true);
    fprintf('\nCoeff err=%1.2e',E2)
end

toggle_plot = 1;
if toggle_plot
figure(1);
m = 64;

colormap(bone(m))
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
            imagesc(U_exact{j}(:,:,1+(k-1)*ceil(length(xs{end})/Uobj.dims(end))))
            drawnow
        end
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

%% plot wendy statistics
figure(2);
w_plot = w_its(:,end);
plot_wendy;