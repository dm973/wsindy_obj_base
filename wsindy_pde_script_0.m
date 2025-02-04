%% load data

dr = 'C:\Users\385570\Desktop\data\wsindy_pde_data\';
pde_name = 'Nav_Stokes.mat';
load([dr,pde_name],'U_exact','lhs','true_nz_weights','xs')

%% create data object

Uobj = wsindy_data(U_exact,xs);

Uobj.coarsen([2 2 2]);
Uobj.addnoise(0.1);

%% lhs 

lhsterms = term('ftag',[0,0,1],'linOp',[0,0,1]); %%% D^\alpha f(u,x)
numeq = size(lhsterms,1);

%% get library

%%% differential operators
x_diffs = [0:2];

%%% poly/trig functions
polys = [0:3]; trigs = [];
custom_terms = {};

lib = trig_poly_lib(polys,trigs,x_diffs,Uobj.nstates,Uobj.ndims,numeq,[]);

%% get testfcn 

tf = testfcn(Uobj);

%% get wsindy model

WS = wsindy_model(Uobj,lib,tf);

%% solve for coefficients

[WS,loss_wsindy,its,G,b] = WS_opt().MSTLS(WS);

WS.disp_mod{:}

function [lib,true_S] = trig_poly_lib(polys,trigs,x_diffs,nstates,ndims,numeq,true_nz_weights)
    tags_poly = get_tags(polys,[],nstates);
    % tags_1 = tags_1(tags_1(:,3)>0,:);
    tags_trig = get_tags([],trigs,nstates);
    tags = prodtags(tags_poly,tags_trig);
    lib = library('nstates',nstates);
    
    diff_tags = get_tags(x_diffs,[],ndims);
    diff_tags = diff_tags(diff_tags(:,end)==0,:);
    % diff_tags = diff_tags(sum(diff_tags>0,2)<=1,:);
    true_S = cell(numeq,1);
    for i=1:size(diff_tags,1)
        for j=1:length(tags)
            if ~and(sum(diff_tags(i,:))>0,isequal(tags{j},zeros(1,nstates)))
                lib.add_terms(term('ftag',tags{j},'linOp',diff_tags(i,:)));

                %%% get true weights if specified
                if exist('true_nz_weights','var')
                    if ~isempty(true_nz_weights)
                        S = cellfun(@(tnz) ismember(tnz(:,1:end-1),[tags{j} diff_tags(i,:)],'rows'), true_nz_weights(:), 'Un',0);
                        for jj=1:length(S)
                            if any(S{jj})
                                true_S{jj} = [true_S{jj};[length(lib.terms) true_nz_weights{jj}(S{jj},end)]];
                            end
                        end
                    end
                end
            end
        end
    end
end

