%% load data

dr = 'C:\Users\385570\Desktop\data\wsindy_pde_data\';
pde_name = 'Nav_Stokes.mat';
load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')

%% create data object

Uobj = wsindy_data(U_exact,xs);

%%% coarsen spacetime grid
Uobj.coarsen([2 2 1]);

%%% add noise
Uobj.addnoise(0.0);

%%% scale data
Uobj.set_scales([]);

%% set lhs 

lhsterms = term('ftag',[0,0,1],'linOp',[0,0,1]); %%% term objects are of the form L*f(u,x) with f specified by 'ftag' of 'fHandle', and L by 'linOp'

%% set library

%%% differential operators
x_diffs = [0:2];

%%% poly/trig functions
polys = [0:2]; trigs = [];

%%% customize
custom_add = {...
    term('fHandle',@(u,v,w) exp(sin(u+v))),...                              % an exotic term specified by function handle    
    compterm(@(u)u.^2,diffOp([1,1,0])),...                                  % a term nonlinear in a derivative
    prodterm(term('ftag',[0,0,-1i]),diffOp([0,0,2])),...                    % a product of two terms
    addterm(diffOp([4,0,0],'stateind',3),diffOp([0,4,0],'stateind',3)),...  % add two terms together
    };

custom_remove_f = {@(tag) all(tag(4:5))};                   % remove all cross derivatives
custom_remove_t = [1 0 0 1 0 0; 0 1 0 0 1 0];             % remove tags for divergence terms

lib = get_lib(Uobj,polys,trigs,x_diffs,custom_add,custom_remove_f,custom_remove_t);

%% set testfcn 

tf = testfcn(Uobj,'meth','FFT','param',2,'stateind',3,'subinds',-2);

%% get wsindy model

WS = wsindy_model(Uobj,lib,tf,'lhsterms',lhsterms);

%% solve for coefficients

[WS,loss_wsindy,its,G,b] = WS_opt().MSTLS(WS,'lambdas',10.^linspace(-4,0,100));

%% view results

%%% display model
Str_mod = WS.disp_mod;
for j=1:WS.numeq
    fprintf('\n----------Eq %u----------\n',j)
    fprintf('%s=',WS.lhsterms{j}.get_str)
    cellfun(@(s)fprintf('%s \n',s),Str_mod{j})
end

%%% plot MSTLS loss
figure(1);clf;
f = min(loss_wsindy(1,:));
g = min(loss_wsindy(2,loss_wsindy(1,:)==f));
for j=1:size(loss_wsindy,1)-1
    loglog(loss_wsindy(end,:),loss_wsindy(j,:),'o-',g,f,'rx')
    hold on
end
hold off
legend;

figure(2);clf
%%% plot residual
for j=1:WS.numeq
    subplot(WS.numeq,1,j)
    plot([WS.bs{1}{j} WS.Gs{1}{j}*WS.reshape_w{j}])
    legend('b','G*w')
    title(['rel. resid=',num2str(norm(WS.Gs{1}{j}*WS.reshape_w{j}-WS.bs{1}{j})/norm(WS.bs{1}{j}))])
end

function lib = get_lib(Uobj,polys,trigs,x_diffs,custom_add,custom_remove_f,custom_remove_t)    
    nstates = Uobj.nstates;
    ndims = Uobj.ndims;
    
    tags_poly = get_tags(polys,[],nstates);
    tags_trig = get_tags([],trigs,nstates);
    tags = prodtags(tags_poly,tags_trig);
    lib = library('nstates',nstates);
    
    diff_tags = get_tags(x_diffs,[],ndims);
    diff_tags = diff_tags(diff_tags(:,end)==0,:);
    for i=1:size(diff_tags,1)
        for j=1:length(tags)
            if all([~and(sum(diff_tags(i,:))>0,...
                    isequal(tags{j},zeros(1,nstates))),...
                    ~cellfun(@(b)b([tags{j} diff_tags(i,:)]),custom_remove_f),...
                    ~ismember_rows([tags{j} diff_tags(i,:)],custom_remove_t)])
                lib.add_terms(term('ftag',tags{j},'linOp',diff_tags(i,:)));
            end
        end
    end
    lib.add_terms(custom_add);
end
