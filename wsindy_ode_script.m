%% load data
nstates = 4;
tol_ode = 10^-10;
t=linspace(0,10,100);
p = [-1,1,3,2,nstates];
rhs_true = @(x) rhs_p(t,x,p);

ntraj = 50;
xcell = cell(ntraj,1);
tcell = cell(ntraj,1);
for i=1:ntraj
    x0 = randn(nstates,1)*1;
    % x0 = (rand(nstates,1)-0.5)*6;
    options_ode_sim = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
    [t,x]=ode45(@(t,x)rhs_true(x),t,x0,options_ode_sim);
    xcell{i} = x;
    tcell{i} = t;
end

%% get wsindy_data object
Uobj = cellfun(@(x,t)wsindy_data(x,t),xcell,tcell);

maxtp = 50;
arrayfun(@(U)U.coarsen(ceil(U.dims/maxtp)),Uobj);

ntraj = length(Uobj);
nstates = Uobj.nstates;
M = Uobj.dims;

noise_ratio = 0.01;
rng('shuffle')
rng_seed = rng().Seed; rng(rng_seed);
arrayfun(@(U)U.addnoise(noise_ratio,'seed',rng_seed),Uobj);

scales = mean(cell2mat(arrayfun(@(U)cellfun(@(u)norm(u,inf)/2,U.Uobs),Uobj(:),'un',0)));
scales = [scales,1];
arrayfun(@(U)U.set_scales(scales),Uobj);

figure(1)
for i=1:ntraj
    Uobj(i).plotDyn;
    hold on
end

figure(2)
datall = cell2mat(arrayfun(@(U)[U.Uobs{:}],Uobj(:),'un',0));
for i=1:nstates
    subplot(ceil(nstates/2),2,i)
    histogram(datall(:,i),ceil(size(datall,1)^(2/5)))
end

%% get lib tags

polys = 0:4;
tags = get_tags(polys,[],nstates);
lib = library('tags',tags);

%% get test function

tf = arrayfun(@(U)testfcn(U,'phifuns',optTFcos(3,0),'meth','FFT','param',1,'subinds',-3),Uobj,'un',0);

%% build WSINDy linear system
WS = wsindy_model(Uobj,lib,tf);

%%% set up rescaling transform
Mscale = arrayfun(@(L)L.get_scales(scales),WS.lib(:),'un',0);
lhs_scales = cellfun(@(t)t.get_scale(scales),WS.lhsterms(:),'un',0);
Mscale = cellfun(@(M,L)M/L,Mscale,lhs_scales,'un',0);
Mscale_W = cell2mat(Mscale);

%%% solve
lambdas = 10.^linspace(-4,0,200);
toggle_jointthresh = 4;
[WS,loss_wsindy,its,G,b] = WS_opt().MSTLS_0(WS,'lambda',lambdas,'toggle_jointthresh',toggle_jointthresh, 'M_diag',Mscale);
W_nd = cellfun(@(w,m)w./m,WS.reshape_w,Mscale,'un',0);

%%% Diagnose
fprintf('\ndata dims =');
arrayfun(@(U)fprintf('%u ',U.dims),Uobj);
fprintf('\n')

Str_mod = WS.disp_mod;
for j=1:WS.numeq
    fprintf('----------Eq %i----------\n',j)
    fprintf('%s=\n',WS.lhsterms{j}.get_str)
    cellfun(@(s)fprintf('%s\n',s),Str_mod{j})
end

fprintf('\n')
resids = cellfun(@(G,w,b)norm(G*w-b)/norm(b),WS.Gs{1},W_nd,WS.bs{1});
arrayfun(@(r)disp(['rel resid=',num2str(r)]),resids);
cellfun(@(s)disp(['sparsity=',num2str(length(s))]),WS.get_supp)
fprintf('\ntf rads=');
cellfun(@(tf)fprintf('%u ',tf.rads),WS.tf{1});
fprintf('\nsize G =')
cellfun(@(G) fprintf('%d ',size(G)), WS.G)
if ~exist('coltrim_fac','var')
    coltrim_fac = 0;
end
if coltrim_fac
    fprintf('\nsize G after coltrim =')
    r_inds = arrayfun(@(j)col_trim_inds{j,loss_wsindy(j,:)==min(loss_wsindy(j,:))},1:WS.numeq,'un',0);
    cellfun(@(l)fprintf('%d ',length(l)),r_inds)
end
fprintf('\ncond G =')
fprintf('%1.1e ', cellfun(@(G)cond(G),WS.G));

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

%% simulate learned and true reduced systems

toggle_compare = 1:ntraj;
toggle_display = 0;
if ~isempty(toggle_compare)
    rhs_learned = WS.get_rhs('w',cell2mat(W_nd));
    tol_dd = 10^-12;
    options_ode_sim = odeset('RelTol',tol_dd,'AbsTol',tol_dd*ones(1,nstates));
    for i=toggle_compare
        fprintf('\n------------dataset %i-------------',i)
        t_train = Uobj(i).grid{1};
        x0_data = Uobj(i).get_x0([]);
        [t_learned,x_learned]=ode15s(@(t,x)rhs_learned(x),t_train,x0_data,options_ode_sim);
        try
            rel_forward_err = arrayfun(@(j)norm(x_learned(:,j)-Uobj(i).Uobs{j})/norm(Uobj(i).Uobs{j}),1:nstates);
            fprintf('\nrel forward err=')
            fprintf('%1.2e ',rel_forward_err)
        catch
            rel_forward_err = NaN;
        end

        if toggle_display
            figure(7+i);clf
            for j=1:nstates
                subplot(nstates,1,j)
                plot(Uobj(i).grid{1},Uobj(i).Uobs{j},'b-o',t_learned,x_learned(:,j),'r-.','linewidth',2)
                legend({'data','learned'})
                title(['rel err=',num2str(rel_forward_err(j))])
            end
        end
    end
end

function dxdt = rhs_p(t,x,p)
    n = p(end);
    dxdt = p(1)*x(:).^p(3);
    for i=1:n-1
        dxdt(i) = dxdt(i) + p(2)*x(i).*x(i+1).^p(4);
    end
    dxdt(n) = dxdt(n) + p(2)*x(n).*x(1).^p(4);
end