classdef WS_opt < handle

    properties
        meth
        diag_reg
    end

    methods
    
        function obj = WS_opt(varargin)
            default_meth = @obj.MSTLS;
    
            p = inputParser;
            addParameter(p,'meth',default_meth);
            parse(p,varargin{:})
    
            obj.meth = p.Results.meth;
        end

        function WS = run_default_meth(obj,WS)
            WS = obj.meth(WS);
        end

        function WS = ols(obj,WS,varargin)
            WS.cat_Gb;
            p = inputParser;
            addParameter(p,'S',[]);
            addParameter(p,'linregargs',{});
            parse(p,varargin{:})
            S = p.Results.S;
            linregargs = p.Results.linregargs;
            ii = cellfun(@(a)isequal(a,'S'),linregargs);
            if any(ii)
                S = linregargs{find(ii)+1};
            end
            if isempty(S)
                S = cellfun(@(g)true(size(g,2),1),WS.G,'uni',0);
            elseif and(isequal(S,0),~isempty(WS.weights))
                if isequal(WS.catm,'blkdiag')
                    S = {WS.weights~=0};    
                else
                    S = cellfun(@(w)w~=0,WS.reshape_w,'uni',0);
                end
            end
            w = cell2mat(cellfun(@(g,b,s)obj.linreg(g(:,s),b,linregargs{:},'S',s), WS.G, WS.b, S, 'uni',0));
            WS.add_weights(w,'toggle_cov',0);
        end

        function y = inject_sparse(obj,w,S)
            y = S*0;
            y(S) = w;
        end

        function WS = ols_tf(obj,WS,varargin)
            WS = obj.ols(WS);

            default_res_tol = 1.5; % if negative, set directly, otherwise use WS.get_theoryres
            default_m_inc = 0.1;
            default_subinds = [];
            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'res_tol',default_res_tol);
            addParameter(p,'subinds',default_subinds);
            addParameter(p,'m_inc',default_m_inc);
            parse(p,WS,varargin{:});
            res_tol = p.Results.res_tol;
            m_inc = p.Results.m_inc;
            subinds = p.Results.subinds;

            for i=1:WS.ntraj
                for j=1:WS.nstates
                    if res_tol>0
                        tol = res_tol*WS.get_theoryres(i,j);
                    else
                        tol = -res_tol;
                    end
                    w = WS.Gs{i}{j} \ WS.bs{i}{j};                    
                    check = norm(WS.Gs{i}{j}*w-WS.bs{i}{j})/norm(WS.bs{i}{j});
                    iter = 0;
%                     disp([iter check tol])
                    foo = 0;
                    while and(check>tol,~any(foo))
                        WS.tf{i}{j}.subinds = subinds;
                        WS.tf{i}{j}.meth = 'direct';
                        WS.tf{i}{j}.param = ceil((1+m_inc)*WS.tf{i}{j}.rads);
                        WS.tf{i}{j}.get_rads(WS.dat(i));
                        WS.tf{i}{j}.get_subinds(WS.dat(i));
                        WS.tf{i}{j}.Cfs = cell(WS.ndims,1);
                        WS.get_Gb_ij(i,j);
                        disp(['rad (',num2str(i),',',num2str(j),')=',num2str(WS.tf{i}{j}.rads)]);
                        w = WS.Gs{i}{j} \ WS.bs{i}{j};
                        check = norm(WS.Gs{i}{j}*w-WS.bs{i}{j})/norm(WS.bs{i}{j});
                        foo = WS.tf{i}{j}.rads==WS.tf{i}{j}.mtmax;
                        if res_tol>0
                            tol = res_tol*WS.get_theoryres(i,j);
                        else
                            tol = -res_tol;
                        end
                        iter = iter +1;
                        disp([iter check tol])
                    end
                end
            end
            WS = obj.ols(WS);
        end

        function [WS,loss_wsindy,its,G,b,col_trim_inds] = MSTLS_0(obj,WS,varargin)
            if isempty(WS.G)
                WS.cat_Gb;
            end
            G = WS.G; b = WS.b;
            if isempty(WS.weights)
                WS.weights = zeros(sum(cellfun(@(G)size(G,2), WS.Gs{1})), 1 );
            end
            default_M_diag = cellfun(@(G) ones(size(G,2),1),G,'uni',0);

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'lambdas',[]);
            addParameter(p,'maxits',inf);
            addParameter(p,'alpha',0.01);
            addParameter(p,'gamma',0);
            addParameter(p,'M_diag',default_M_diag);
            addParameter(p,'toggle_jointthresh',1);
            addParameter(p,'linregargs',repmat({{'verbose','none'}},WS.numeq,1));
            addParameter(p,'incl_inds',cell(WS.numeq,1));
            addParameter(p,'coltrim',0);
            addParameter(p,'subset_eq',1:length(G));
            addParameter(p,'toggle_discrep',0);

            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            alpha = arrayfun(@(L) (p.Results.alpha*length(L.terms)+1)^-1,WS.lib);
            gamma = p.Results.gamma;
            M_diag = p.Results.M_diag;
            if isempty(M_diag)
                M_diag = default_M_diag;
            end
            toggle_jointthresh = p.Results.toggle_jointthresh;
            linregargs = p.Results.linregargs;
            incl_inds = p.Results.incl_inds;
            toggle_coltrim = p.Results.coltrim;
            toggle_discrep = p.Results.toggle_discrep;
            subset_eq = p.Results.subset_eq;

            lambdas = p.Results.lambdas;
            if isempty(lambdas)
                if toggle_jointthresh==1
                    lambdas = 10.^(linspace(max(log10(abs(b{1}'*G{1}./vecnorm(G{1}).^2)'))-4,max(log10(abs(b{1}'*G{1}./vecnorm(G{1}).^2)')),100));
                elseif toggle_jointthresh==4
                    lambdas = abs(b{1}'*G{1}./vecnorm(G{1}).^2)';
                else
                    lambdas = 10.^linspace(-4,0,100);
                end
            end

            if and(toggle_discrep==1,~isempty(WS.weights))                
                b = cellfun(@(g,b,w) b-g*w,G,b,WS.reshape_w,'uni',0);
                G = cellfun(@(g,w) g(:,~w),G,WS.reshape_w,'uni',0);
                % fix!!!
                % if ~isempty(Aeq)
                %     deq = cellfun(@(d,A,w,M) d-(A.*M')*w,deq,Aeq,WS.reshape_w,M_diag,'uni',0);
                %     Aeq = cellfun(@(A,w) A(:,~w),Aeq,WS.reshape_w,'uni',0);
                % end
                M_diag = cellfun(@(M,w)M(~w),M_diag,WS.reshape_w,'uni',0);
            elseif and(toggle_discrep==2,~isempty(WS.weights))
                b = cellfun(@(g,b,w) b-g*w,G,b,WS.reshape_w,'uni',0);
            end

            W_ls = cellfun(@(g,b,LRA) obj.linreg(g,b,LRA{:}), G, b, linregargs, 'uni',0);
            GW_ls = cellfun(@(g,w) norm(g*w), G,W_ls, 'uni',0);            

            W = cellfun(@(G) zeros(size(G,2),1), G,'un',0);
            its = zeros(length(G),1);
            loss_wsindy = zeros(length(G)+1,length(lambdas));
            col_trim_inds = {};
            for k=subset_eq
                W_all = zeros(size(G{k},2),length(lambdas));
                for l=1:length(lambdas)
                    if toggle_coltrim
                        if l>1
                            inds = coltrim(G{k}(:,col_trim_inds{k,l-1}),1-toggle_coltrim*lambdas(l),b{k},incl_inds_temp);
                            inds = col_trim_inds{k,l-1}(inds);
                        else
                            inds = coltrim(G{k},1-toggle_coltrim*lambdas(l),b{k},incl_inds{k});
                        end
                        col_trim_inds{k,l} = inds;
                        G_temp = G{k}(:,inds);
                        M_temp = M_diag{k}(inds);
                        if ~isequal(incl_inds{k},'all')
                            incl_inds_temp = find(ismember(inds,incl_inds{k}));
                        else
                            incl_inds_temp = incl_inds{k};
                        end
                    else
                        G_temp = G{k};
                        inds = 1:size(G_temp,2);
                        M_temp = M_diag{k};
                        incl_inds_temp = incl_inds{k};
                    end
                    [w_temp,~] = obj.sparsifyDynamics(G_temp, b{k}, lambdas(l), gamma, M_temp, maxits, toggle_jointthresh, linregargs{k},incl_inds_temp);
                    W_all(inds,l) = w_temp;
                end
                proj_cost = alpha(k)*vecnorm(G{k}*(W_all./M_diag{k}-W_ls{k}))/GW_ls{k};
                overfit_cost = (1-alpha(k))*arrayfunvec(W_all,@(w)length(find(w)),1)/length(find(W_ls{k}));
                lossvals = proj_cost + overfit_cost;
                W{k} = W_all(:,find(lossvals == min(lossvals),1));
                loss_wsindy(k,:) = lossvals;
                weight_inds = sum(cellfun(@(G)size(G,2),WS.Gs{1}(1:k-1)))+1:sum(cellfun(@(G)size(G,2),WS.Gs{1}(1:k)));

                if and(toggle_discrep==1,any(WS.weights(weight_inds)))
                    wtemp = WS.weights(weight_inds);
                    wtemp(~WS.weights(weight_inds)) = W{k};
                    WS.weights(weight_inds) = wtemp;
                elseif and(toggle_discrep==2,any(WS.weights(weight_inds)))
                    WS.weights(weight_inds) = WS.weights(weight_inds)+W{k};
                else
                    WS.weights(weight_inds) = W{k};
                end

            end
            loss_wsindy(end,:) = lambdas;
        end

        function [WS,loss_wsindy,its,lambda_mins] = MSTLS_group(obj,WS,varargin)
            Gs = {}; bs = {}; default_M_diags = {};
            P = length(WS);
            for W = WS
                if isempty(W.G)
                    W.cat_Gb;
                end
                if isempty(W.weights)
                    W.weights = zeros(sum(cellfun(@(G)size(G,2), W.Gs{1})), 1 );
                end
                Gs = [Gs,{W.G}]; bs = [bs,{W.b}];
                default_M_diags = [default_M_diags ,{cellfun(@(G) ones(size(G,2),1),Gs,'uni',0)}];
            end

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'lambdas',0);
            addParameter(p,'maxits',inf);
            addParameter(p,'alpha',0.01);
            addParameter(p,'gamma',0);
            addParameter(p,'toggle_sign',false);
            addParameter(p,'M_diag',default_M_diags);
            addParameter(p,'toggle_jointthresh',1);
            addParameter(p,'linregargs',[]);
            addParameter(p,'incl_inds',cell(WS(1).numeq,1));
            addParameter(p,'subset_eq',1:WS(1).numeq);
            addParameter(p,'toggle_discrep',0);

            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            toggle_sign = p.Results.toggle_sign;
            alpha = arrayfun(@(L) (p.Results.alpha*length(L.terms)+1)^-1,WS(1).lib);
            gamma = p.Results.gamma;
            M_diags = p.Results.M_diag;
            if isempty(M_diags)
                M_diags = default_M_diags;
            end
            toggle_jointthresh = p.Results.toggle_jointthresh;
            linregargs = p.Results.linregargs;
            incl_inds = p.Results.incl_inds;
            toggle_discrep = p.Results.toggle_discrep;
            subset_eq = p.Results.subset_eq;

            lambdas = p.Results.lambdas;

            if and(toggle_discrep==1,~isempty(WS(1).weights))
                for p=1:P
                    bs{p} = cellfun(@(g,b,w) b-g*w,Gs{p},bs{p},WS(p).reshape_w,'uni',0);
                    Gs{p} = cellfun(@(g,w) g(:,~w),Gs{p},WS(p).reshape_w,'uni',0);
                % fix!!!
                % if ~isempty(Aeq)
                %     deq = cellfun(@(d,A,w,M) d-(A.*M')*w,deq,Aeq,WS.reshape_w,M_diag,'uni',0);
                %     Aeq = cellfun(@(A,w) A(:,~w),Aeq,WS.reshape_w,'uni',0);
                % end
                    M_diags{p} = cellfun(@(M,w)M(~w),M_diags{p},WS(p).reshape_w,'uni',0);
                end
            elseif and(toggle_discrep==2,~isempty(WS(1).weights))
                for p=1:P
                    bs{p} = cellfun(@(g,b,w) b-g*w,Gs{p},bs{p},WS(p).reshape_w,'uni',0);
                end
            end

            W_ls = cellfun(@(G,b,LRA) cellfun(@(g,b,L) obj.linreg(g,b,L{:}), G,b,LRA,'un',0), Gs, bs, linregargs, 'un',0);
            GW_ls = cellfun(@(G,W) cellfun(@(g,w) norm(g*w), G, W, 'uni',0), Gs, W_ls,'un',0);            

            W = cellfun(@(G) zeros(size(G,2),P), Gs{1},'un',0);
            its = zeros(length(Gs{1}),1);
            loss_wsindy = zeros(WS(1).numeq+1,length(lambdas));
            lambda_mins = zeros(WS(1).numeq,1);
            col_trim_inds = {};
            for k=subset_eq
                W_all = zeros(size(Gs{1}{k},2),P,length(lambdas));
                G_k = cellfun(@(G)G{k},Gs,'Un',0);
                b_k = cellfun(@(b)b{k},bs,'Un',0);
                M_k = cellfun(@(M)M{k},M_diags,'Un',0);
                GW_ls_k = cellfun(@(GW)GW{k},GW_ls,'Un',0);
                W_ls_k = cellfun(@(W)W{k},W_ls,'Un',0);
                linregargs_k = cellfun(@(L)L{k},linregargs,'Un',0);
                incl_inds_k = incl_inds{k};
                
                for l=1:length(lambdas)
                    [w_temp,~] = obj.sparsifyGroupDynamics(G_k, b_k, lambdas(l), gamma, M_k, maxits, toggle_jointthresh, linregargs_k,incl_inds_k,toggle_sign);
                    W_all(:,:,l) = w_temp;
                end

                proj_cost = 0; overfit_cost = 0;
                for p=1:P
                    proj_cost = proj_cost + vecnorm(G_k{p}*(squeeze(W_all(:,p,:))./M_k{p}-W_ls_k{p}))/GW_ls_k{p};
                    overfit_cost = overfit_cost + arrayfunvec(squeeze(W_all(:,p,:)),@(w)length(find(w)),1)/length(find(W_ls_k{p}));
                end

                proj_cost = alpha(k)*proj_cost;
                overfit_cost = (1-alpha(k))*overfit_cost;
                lossvals = proj_cost + overfit_cost;
                W{k} = W_all(:,:,find(lossvals == min(lossvals),1));
                lambda_mins(k) = lambdas(find(lossvals == min(lossvals),1));
                loss_wsindy(k,:) = lossvals;

                for p=1:P
                    weight_inds = sum(cellfun(@(G)size(G,2),WS(p).Gs{1}(1:k-1)))+1:sum(cellfun(@(G)size(G,2),WS(p).Gs{1}(1:k)));
                    if and(toggle_discrep==1,any(WS(p).weights(weight_inds)))
                        wtemp = WS(p).weights(weight_inds);
                        wtemp(~WS.weights(weight_inds)) = W{k}(:,p);
                        WS(p).weights(weight_inds) = wtemp;
                    elseif and(toggle_discrep==2,any(WS(p).weights(weight_inds)))
                        WS(p).weights(weight_inds) = WS(p).weights(weight_inds)+W{k}(:,p);
                    else
                        WS(p).weights(weight_inds) = W{k}(:,p);
                    end
                end

            end
            loss_wsindy(end,:) = lambdas;
        end

        function [WS,loss_wsindy,its,G,b,r_inds] = MSTLS(obj,WS,varargin)
            if isempty(WS.G)
                WS.cat_Gb;
            end
            G = WS.G; b = WS.b;

            if WS.toggleH
                default_M_diag = {ones(length(WS.lib.terms),1)};
            else
                default_M_diag = cellfun(@(G) ones(size(G,2),1),G,'uni',0);
            end

            default_lambdas = [];
            default_maxits = inf;
            default_alpha = 0.01;
            default_gamma = 0;
            default_Hreg = 0;
            default_applycov = 0; 
            default_toggle_discrep = 0;
            default_toggle_jointthresh = 1;
            default_reg0 = 0;

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'lambdas',default_lambdas);
            addParameter(p,'maxits',default_maxits);
            addParameter(p,'alpha',default_alpha);
            addParameter(p,'gamma',default_gamma);
            addParameter(p,'M_diag',default_M_diag);
            addParameter(p,'Hreg',default_Hreg);
            addParameter(p,'applycov',default_applycov);
            addParameter(p,'toggle_discrep',default_toggle_discrep);
            addParameter(p,'toggle_jointthresh',default_toggle_jointthresh);
            addParameter(p,'reg0',default_reg0);
            addParameter(p,'Aeq',[]);
            addParameter(p,'deq',[]);
            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            alpha = arrayfun(@(L) (p.Results.alpha*length(L.terms)+1)^-1,WS.lib);
            gamma = p.Results.gamma;
            M_diag = p.Results.M_diag;
            if isempty(M_diag)
                M_diag = default_M_diag;
            end
            Hreg = p.Results.Hreg;
            applycov = p.Results.applycov;
            toggle_discrep = p.Results.toggle_discrep;
            toggle_jointthresh = p.Results.toggle_jointthresh;
            reg0 = p.Results.reg0;
            Aeq = p.Results.Aeq;
            deq = p.Results.deq;

            if applycov==1
                for k=1:length(G)
                    [G{k},b{k}] = WS.apply_cov(G{k},b{k},obj.diag_reg);
                end
            end

            lambdas = p.Results.lambdas;

            if isempty(lambdas)
                if toggle_jointthresh==1
                    lambdas = 10.^(linspace(max(log10(abs(b{1}'*G{1}./vecnorm(G{1}).^2)'))-4,max(log10(abs(b{1}'*G{1}./vecnorm(G{1}).^2)')),100));
                else
                    lambdas = 10.^linspace(-4,0,100);
                end
            end

            if and(Hreg>0,WS.toggleH)
                Hlib = cellfun(@(tm) tm.fHandle, WS.lib.terms, 'uni',0);
                for j=1:WS.ntraj
                    regmat = conv2(WS.tf{j}{1}.Cfs{1}(2,:),1,cell2mat(cellfun(@(H)H(WS.dat(j).Uobs{:}),Hlib(:)','uni',0)),'valid');
                    regmat = Hreg*regmat(1:mean(diff(WS.tf{j}{1}.subinds{1})):end,:);
                    regb = zeros(size(regmat,1),1);
                    G = {[G{1};regmat]};b = {[b{1};regb]};
                end
            end

            if and(toggle_discrep==1,~isempty(WS.weights))                
                b = cellfun(@(g,b,w) b-g*w,G,b,WS.reshape_w,'uni',0);
                G = cellfun(@(g,w) g(:,~w),G,WS.reshape_w,'uni',0);
                if ~isempty(Aeq)
                    deq = cellfun(@(d,A,w,M) d-(A.*M')*w,deq,Aeq,WS.reshape_w,M_diag,'uni',0);
                    Aeq = cellfun(@(A,w) A(:,~w),Aeq,WS.reshape_w,'uni',0);
                end
                M_diag = cellfun(@(M,w)M(~w),M_diag,WS.reshape_w,'uni',0);
            elseif and(toggle_discrep==2,~isempty(WS.weights))
                b = cellfun(@(g,b,w) b-g*w,G,b,WS.reshape_w,'uni',0);
            end 

            if ~isempty(Aeq)
                W_ls = cellfun(@(g,b,A,d) lineqcorrect(g,b,A,d), G,b,Aeq,deq, 'uni',0);
            else
                W_ls = cellfun(@(g,b) g \ b, G,b, 'uni',0);
            end

            GW_ls = cellfun(@(g,w) norm(g*w), G,W_ls, 'uni',0);            

            if isnumeric(reg0) & reg0>0 % initial guess: truncated Least-squares
                reg_inds = cellfun(@(g,b) abs(b'*g)./vecnorm(g)/norm(b),G,b,'uni',0);
                for i=1:length(G)
                    [~,reg_inds{i}] = sort(reg_inds{i},'descend');
                    reg_inds{i} = reg_inds{i}(1:min(reg0,end));
%                     disp(['log10(cond(G(reg0)))=',num2str(log10(cond(G{i}(:,reg_inds{i}))))])
%                     disp(['initial coeff. magnitude range (log10 scale)=',num2str(range(log10(abs(G{i}(:,reg_inds{i}) \ b{i}))))])
                end
            elseif isequal(reg0,'coltrim')
                reg_inds = repmat({'coltrim'},1,length(G));
            else
                reg_inds = repmat({[]},1,length(G));
            end

            W = cell(length(G),1);
            r_inds = cell(length(G),length(lambdas));
            its = zeros(length(G),1);
            loss_wsindy = zeros(length(G)+1,length(lambdas));
            for k=1:length(G)
                W_all = zeros(size(G{k},2),length(lambdas));
                for l=1:length(lambdas)
                    if isempty(Aeq)
                        [W_all(:,l),~,r_inds{k,l}] = sparsifyDynamics(G{k}, b{k}, lambdas(l), 1, gamma, M_diag{k}, maxits, toggle_jointthresh, reg_inds{k});
                    else
                        [W_all(:,l),~,r_inds{k,l}] = sparsifyDynamics_lineq(G{k}, b{k}, lambdas(l), 1, gamma, M_diag{k}, maxits, toggle_jointthresh, reg_inds{k},Aeq{k},deq{k});
                    end
                end
                if ~isempty(M_diag)
                    proj_cost = alpha(k)*vecnorm(G{k}*(W_all./M_diag{k}-W_ls{k}))/GW_ls{k};
                else
                    proj_cost = alpha(k)*vecnorm(G{k}*(W_all-W_ls{k}))/GW_ls{k};
                end
                overfit_cost = (1-alpha(k))*arrayfunvec(W_all,@(w)length(find(w)),1)/length(find(W_ls{k}));
                lossvals = proj_cost + overfit_cost;            
                W{k} = W_all(:,find(lossvals == min(lossvals),1));
                loss_wsindy(k,:) = lossvals;
            end
            loss_wsindy(end,:) = lambdas;
            if and(toggle_discrep==1,~isempty(WS.weights))
                wtemp = WS.weights;
                wtemp(~WS.weights) = cell2mat(W);
                WS.weights = wtemp;
            elseif and(toggle_discrep==2,~isempty(WS.weights))
                WS.weights = WS.weights+cell2mat(W);
            else
                WS.weights = cell2mat(W);
            end
        end

        function [W,G,b] = MSTLS_param(obj,WS,P,varargin)
            %%% P is a ntraj x numparambasis matrix
            G = cell(1,WS.nstates);
            for j=1:WS.nstates
                for i=1:WS.ntraj
                    G{j} = [G{j};kron(WS.Gs{i}{j},P(i,:))];
                end
            end
            b = arrayfun(@(i)cell2mat(cellfun(@(b)b{i},WS.bs,'uni',0)),1:WS.numeq,'uni',0);

            default_M_diag = cellfun(@(G) ones(size(G,2),1),G,'uni',0);

            default_lambdas = [];
            default_maxits = inf;
            default_alpha = 0.01;
            default_gamma = 0;
            default_toggle_jointthresh = 1;
            default_reg0 = 0;

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'lambdas',default_lambdas);
            addParameter(p,'maxits',default_maxits);
            addParameter(p,'alpha',default_alpha);
            addParameter(p,'gamma',default_gamma);
            addParameter(p,'M_diag',default_M_diag);
            addParameter(p,'toggle_jointthresh',default_toggle_jointthresh);
            addParameter(p,'reg0',default_reg0);
            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            alpha = (p.Results.alpha*mean(arrayfun(@(L)length(L.terms)*size(P,2),WS.lib))+1)^-1;
            gamma = p.Results.gamma;
            M_diag = p.Results.M_diag;
            toggle_jointthresh = p.Results.toggle_jointthresh;
            reg0 = p.Results.reg0;

            lambdas = p.Results.lambdas;

            if isempty(lambdas)
                if toggle_jointthresh==1
                    lambdas = 10.^(linspace(max(log10(abs(b{1}'*G{1}./vecnorm(G{1}).^2)'))-4,max(log10(abs(b{1}'*G{1}./vecnorm(G{1}).^2)')),100));
                else
                    lambdas = 10.^linspace(-4,0,100);
                end
            end

            W_ls = cellfun(@(g,b) g \ b, G,b, 'uni',0);
            GW_ls = cellfun(@(g,w) norm(g*w), G,W_ls, 'uni',0);            

            if reg0>0 % initial guess: truncated Least-squares
                reg_inds = cellfun(@(g,b) abs(b'*g)./vecnorm(g)/norm(b),G,b,'uni',0);
                for i=1:length(G)
                    [~,reg_inds{i}] = sort(reg_inds{i},'descend');
                    reg_inds{i} = reg_inds{i}(1:min(reg0,end));
%                     disp(['log10(cond(G(reg0)))=',num2str(log10(cond(G{i}(:,reg_inds{i}))))])
%                     disp(['initial coeff. magnitude range (log10 scale)=',num2str(range(log10(abs(G{i}(:,reg_inds{i}) \ b{i}))))])
                end
            else
                reg_inds = repmat({[]},1,length(G));
            end

            W = cell(length(G),1);
            its = zeros(length(G),1);
            loss_wsindy = zeros(length(G)+1,length(lambdas));
            for k=1:length(G)
                W_all = zeros(size(G{k},2),length(lambdas));
                for l=1:length(lambdas)
                    [W_all(:,l),~] = sparsifyDynamics(G{k}, b{k}, lambdas(l), 1, gamma, M_diag{k}, maxits, toggle_jointthresh, reg_inds{k});
                end
                if ~isempty(M_diag)
                    proj_cost = alpha*vecnorm(G{k}*(W_all./M_diag{k}-W_ls{k}))/GW_ls{k};
                else
                    proj_cost = alpha*vecnorm(G{k}*(W_all-W_ls{k}))/GW_ls{k};
                end
                overfit_cost = (1-alpha)*arrayfunvec(W_all,@(w)length(find(w)),1)/length(find(W_ls{k}));
                lossvals = proj_cost + overfit_cost;            
                W{k} = W_all(:,find(lossvals == min(lossvals),1));
                loss_wsindy(k,:) = lossvals;
            end
            loss_wsindy(end,:) = lambdas;
            W = cellfun(@(w)reshape(w,size(P,2),[])',W,'uni',0);
        end

        function x = linreg(obj,A,b,varargin)
            p = inputParser;
            addRequired(p,'A');
            addRequired(p,'b');
            addParameter(p,'S',true(size(A,2),1));
            addParameter(p,'x0',[]);
            addParameter(p,'Aineq',[]);
            addParameter(p,'bineq',[]);
            addParameter(p,'Aeq',[]);
            addParameter(p,'beq',[]);
            addParameter(p,'LB',[]);
            addParameter(p,'UB',[]);
            addParameter(p,'consttol',10^-10);
            addParameter(p,'opttol',10^-10);
            addParameter(p,'maxits',1000);
            addParameter(p,'verbose','none');
            parse(p,A,b,varargin{:})

            x0 = p.Results.x0;
            Aineq = p.Results.Aineq;
            bineq = p.Results.bineq;
            if isempty(Aineq)
                Aineq = [];
                bineq = [];
            end
            Aeq = p.Results.Aeq;
            beq = p.Results.beq;
            if isempty(Aeq)
                Aeq = [];
                beq = [];
            end
            LB = p.Results.LB;
            UB = p.Results.UB;
            consttol = p.Results.consttol;
            opttol = p.Results.opttol;
            maxits = p.Results.maxits;
            S = p.Results.S;
            verbosity = p.Results.verbose;

            if any(S)
                if size(A,2)~=length(find(S))
                    A = A(:,S);
                end

                if isempty(x0)
                    if diff(size(A))>=0
                        % 
                        % reg0 = rank(A,norm(A)*10^-4);
                        % reg_inds = abs(b'*A)./vecnorm(A)/norm(b);
                        % [~,reg_inds] = sort(reg_inds,'descend');
                        % reg_inds = reg_inds(1:min(reg0,end));
                        % x0 = zeros(size(A,2),1);
                        % x0(reg_inds) = A(:,reg_inds) \ b;
                        % x0 = obj.inject_sparse(x0,S);

                        x0 = obj.inject_sparse(lsqminnorm(A,b),S);
                    else
                        x0 = obj.inject_sparse(A\b,S);
                    end
                end
    
                if any([~isempty(Aineq) ~isempty(bineq) ~isempty(Aeq) ~isempty(beq) ~isempty(LB) ~isempty(UB)])
                    if ~isempty(Aineq)
                        Aineq = Aineq(:,S);
                    end
                    if ~isempty(Aeq)
                        Aeq = Aeq(:,S);
                    end
                    if ~isempty(LB)
                        LB = LB(S);
                    end
                    if ~isempty(UB)
                        UB = UB(S);
                    end
                    % if isequal(verbosity,'None')
                    %     N1 = null(Aeq);
                    %     try
                    %         N = null((Aineq*N1)');
                    %         e = max(abs(bineq'*N));
                    %         if e > 0
                    %             disp(['NO FEASIBLE BOUNDARY POINT: e=',num2str(e)])
                    %         end
                    %     end
                    % end
                    options = optimoptions('quadprog','Display',verbosity,'ConstraintTolerance',consttol,'OptimalityTolerance',opttol,'MaxIterations',maxits);
                    x = quadprog((A'*A),-(A'*b),Aineq,bineq,Aeq,beq,LB,UB,x0,options);
                    if isempty(x)
                        x = zeros(size(A,2),1);
                    end
                    x = obj.inject_sparse(x,S);
                else
                    x = x0;
                end
            else
                x = S*0;
            end

        end

        function [WS,w_its,res,res_0,CovW,RT] = wendy(obj,WS,varargin)
            % options: maxits,ittol,diag_reg,w,regmeth

            default_maxits = 20;
            default_ittol = 10^-4;
            default_diag_reg = 10^-6;
            default_w = WS.weights;
            default_regmeth = 'ols';

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'maxits',default_maxits);
            addParameter(p,'trim_rows',0);
            addParameter(p,'ittol',default_ittol);
            addParameter(p,'diag_reg',default_diag_reg);
            addParameter(p,'w',default_w);
            addParameter(p,'regmeth',default_regmeth);
            addParameter(p,'verbose',0);
            addParameter(p,'linregargs',{});
            addParameter(p,'CovW_type','OLS');
            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            ittol = p.Results.ittol;
            obj.diag_reg = p.Results.diag_reg;
            w = p.Results.w;
            regmeth = p.Results.regmeth;
            verbosity = p.Results.verbose;
            linregargs = p.Results.linregargs;
            tr = p.Results.trim_rows;
            CovW_type = p.Results.CovW_type;

            if verbosity
                tic,
            end

            if ~isempty(w)
                if isequal(w,0)
                    WS = obj.ols(WS,'S',0,'linregargs',linregargs);                    
                else
                    WS.add_weights(w,'toggle_cov',1);
                end
            else
                WS = obj.ols(WS,'linregargs',linregargs);
            end

            if tr>1
                WS.get_Lfac;
                WS.trim_rows('trim_factor',tr);
            end
            
            check = 1;
            sparse_inds = WS.weights~=0;
            w_its = WS.weights;
            res_0 = [];
            res = [];
            its = 0;
            WS.cat_Gb('cat','blkdiag');
            G_0 = WS.G{1};
            b_0 = WS.b{1};


            while and(check,its<maxits)
                % disp(['iter=',num2str(its)])
                if isequal(regmeth,'ols')
                    [G,b,RT] = WS.apply_cov(G_0(:,sparse_inds),b_0,obj.diag_reg);
                    w = obj.linreg(G,b,linregargs{:},'S',sparse_inds);
                    WS.add_weights(w,'toggle_cov',1);
                elseif isequal(regmeth,'MSTLS')
                    [WS,~,~,G,b] = obj.MSTLS(WS,'applycov',1);
                    G = blkdiag(G{:}); b = cell2mat(b);
                    sparse_inds = WS.weights~=0;
                    G = G(:,sparse_inds);
                    RT = chol(WS.cov)';
                end

                w_its = [w_its WS.weights];
                check = norm(diff(w_its(:,end-1:end),[],2))/norm(w_its(:,end-1))>ittol;

                res_0 = [res_0 G_0(:,sparse_inds)*w_its(sparse_inds,end)-b_0];

                res = [res G*w_its(sparse_inds,end)-b];
                its = size(w_its,2)-1;
            end

            Ginv = pinv(G_0(:,WS.weights~=0));
            if ~exist('RT','var')
                try
                    RT = speye(size(Ginv,2))*norm(res_0(:,end))/sqrt(size(res_0,1)-1);
                catch
                    RT = speye(size(Ginv,2));
                end
            end
            
            if isequal(CovW_type, 'C-R')
            %%% Cramer-Rao inspired covariance (more consistent with WENDy parameter distribution)
                CovW = inv(G'*G);
            elseif isequal(CovW_type, 'OLS')
            %%% OLS left-inverse to compute parameter covariance
                Ginv = Ginv*RT;
                CovW = Ginv*Ginv';
            end
            
            if verbosity
                disp(['wendy iter time=',num2str(toc),'; sparsity=',num2str(length(find(WS.weights))),'; its=',num2str(its)])
            end
        end

        function [WS,w_its,res,res_0,CovW,RT] = wendy2(obj,WS,varargin)
            % options: maxits,ittol,diag_reg,w,regmeth

            default_maxits = 20;
            default_ittol = 10^-4;
            default_diag_reg = 10^-6;
            default_w = WS.weights;
            default_regmeth = 'ols';

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'maxits',default_maxits);
            addParameter(p,'trim_rows',0);
            addParameter(p,'ittol',default_ittol);
            addParameter(p,'diag_reg',default_diag_reg);
            addParameter(p,'w',default_w);
            addParameter(p,'regmeth',default_regmeth);
            addParameter(p,'verbose',0);
            addParameter(p,'linregargs',{});
            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            ittol = p.Results.ittol;
            obj.diag_reg = p.Results.diag_reg;
            w = p.Results.w;
            regmeth = p.Results.regmeth;
            verbosity = p.Results.verbose;
            linregargs = p.Results.linregargs;
            tr = p.Results.trim_rows;

            if verbosity
                tic,
            end

            if ~isempty(w)
                if isequal(w,0)
                    WS = obj.ols(WS,'S',0,'linregargs',linregargs);                    
                else
                    WS.add_weights(w,'toggle_cov',1);
                end
            else
                WS = obj.ols(WS,'linregargs',linregargs);
            end

            if tr>1
                WS.get_Lfac;
                WS.trim_rows('trim_factor',tr);
            end
            
            check = 1;
            sparse_inds = WS.weights~=0;
            w_its = WS.weights;
            res_0 = [];
            res = [];
            its = 0;
            WS.cat_Gb('cat','blkdiag');
            G_0 = WS.G{1};
            b_0 = WS.b{1};

            while and(check,its<maxits)
                if isequal(regmeth,'ols')
                    [G,b,RT] = WS.apply_cov(G_0(:,sparse_inds),b_0,obj.diag_reg);

                    v = RT \ (G*WS.weights-b);
                    A = G - (2*RT) \ WS.apply_gradC(v);

                    w = obj.linreg(A'*G,A'*b,linregargs{:},'S',sparse_inds);
                    WS.add_weights(w,'toggle_cov',1);
                end

                w_its = [w_its WS.weights];
                check = norm(diff(w_its(:,end-1:end),[],2))/norm(w_its(:,end-1))>ittol;

                res_0 = [res_0 G_0(:,sparse_inds)*w_its(sparse_inds,end)-b_0];

                res = [res G*w_its(sparse_inds,end)-b];
                its = size(w_its,2)-1;
            end
            
            % CovW = inv(G'*G);
            Ginv = pinv(G_0(:,WS.weights~=0));
            if ~exist('RT','var')
                RT = speye(size(Ginv,2))*norm(res_0(:,end))/sqrt(size(res_0,1)-1);
            end
            Ginv = Ginv*RT;
            CovW = Ginv*Ginv';
            if verbosity
                disp(['wendy iter time=',num2str(toc),'; sparsity=',num2str(length(find(WS.weights))),'; its=',num2str(its)])
            end
        end

        function [WS,loss_wsindy,lambda,w_its,res,res_0,CovW,RT] = MSTLS_WENDy(obj,WS,varargin)

            default_lambdas = 10.^linspace(-4,0,100);
            default_maxits = inf;
            default_alpha = 0.01;

            default_maxits_wendy = 20;
            default_ittol = 10^-6;
            default_diag_reg = 10^-6;

            p = inputParser;
            %%%% MSTLS params
            addParameter(p,'lambdas',default_lambdas);
            addParameter(p,'maxits',default_maxits);
            addParameter(p,'alpha',default_alpha);
            addParameter(p,'M_diag',ones(sum(arrayfun(@(L)length(L.terms),WS.lib)),1));

            %%%% MSTLS params
            addParameter(p,'cov_thresh',0);

            %%%% wendy params
            addRequired(p,'WS');
            addParameter(p,'maxits_wendy',default_maxits_wendy);
            addParameter(p,'ittol',default_ittol);
            addParameter(p,'diag_reg',default_diag_reg);
            addParameter(p,'verbose',0);
            addParameter(p,'linregargs',{});

            parse(p,WS,varargin{:})

            maxits = p.Results.maxits;
            alpha = (p.Results.alpha*mean(arrayfun(@(L)length(L.terms),WS.lib))+1)^-1;
            lambdas = p.Results.lambdas;
            M_diag =  p.Results.M_diag;
            linregargs = p.Results.linregargs;
            
            cov_thresh = p.Results.cov_thresh;

            maxits_wendy = p.Results.maxits_wendy;
            ittol = p.Results.ittol;
            diag_reg = p.Results.diag_reg;
            verbosity = p.Results.verbose;

            if maxits_wendy>0
                vw = {'maxits',maxits_wendy,'ittol',ittol,'diag_reg',diag_reg,'verbose',verbosity,'linregargs',linregargs};
            else
                vw = {'maxits',maxits_wendy,'ittol',ittol,'diag_reg',diag_reg,'verbose',verbosity,'w',0,'linregargs',linregargs};
            end
            if isempty(WS.G)
                WS.cat_Gb('cat','blkdiag');
            end
            G_0 = WS.G{1};
            b_0 = WS.b{1};
            W_ls = obj.linreg(G_0,b_0,linregargs{:});
            bnds = norm(b_0)./vecnorm(G_0)';
            GW_ls = norm(G_0*W_ls);
            WS.add_weights(W_ls,'toggle_cov',1);
            
            proj_cost = []; overfit_cost = []; lossvals = [];
            
            Wmat = zeros(length(W_ls),length(lambdas));
            for l=1:length(lambdas)
                WS.weights = W_ls;
                [WS,its] = obj.sparsifyDynamics_wendy(WS,lambdas(l),M_diag,bnds,maxits,cov_thresh,vw);
                proj_cost = [proj_cost alpha*norm(G_0*(WS.weights-W_ls))/GW_ls];
                overfit_cost = [overfit_cost (1-alpha)*length(find(WS.weights))/length(W_ls)];
                lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
                Wmat(:,l) = WS.weights;
            end
            l = find(lossvals == min(lossvals),1);  
            lambda = lambdas(l);
            WS.weights = Wmat(:,l);
            % if any(Wmat(:,l)~=0)
            %     WS.weights = obj.linreg(G_0(:,Wmat(:,l)~=0),b_0,linregargs{:},'S',Wmat(:,l)~=0);
            % else
            %     WS.weights = Wmat(:,l);
            % end

            [WS,w_its,res,res_0,CovW,RT] = obj.wendy(WS,vw{:});

            WS.add_weights(WS.weights.*M_diag,'toggle_cov',1);
            loss_wsindy = zeros(2,length(lambdas));
            loss_wsindy(1,:) = lossvals;
            loss_wsindy(end,:) = lambdas;
        end

        function [WS,loss_wsindy,its] = MSTLSQP(obj,WS,varargin)
            if isempty(WS.G)
                WS.cat_Gb;
            end
            G = WS.G; 
            b = WS.b;

            if WS.toggleH
                default_M_diag = {ones(length(WS.lib.terms),1)};
            else
                default_M_diag = cellfun(@(G) ones(size(G,2),1),G,'uni',0);
            end

            default_lambdas = 10.^linspace(-4,0,100);
            default_maxits = inf;
            default_alpha = 0.01;
            default_gamma = 0;

            defaultexcl_inds=repmat({[]},1,WS.numeq);
            defaultAineq=repmat({[]},1,WS.numeq);
            defaultbineq=repmat({[]},1,WS.numeq);
            defaultopt_tol = min(cellfun(@(G)1/cond(G'*G),WS.G));
            defaultconst_tol = defaultopt_tol;
            defaultmaxQPits = 1000;
            defaultdispQP='off';
            defaultauto = [];
            
            inp = inputParser;
            addRequired(inp,'WS');
            addParameter(inp,'lambdas',default_lambdas);
            addParameter(inp,'maxits',default_maxits);
            addParameter(inp,'alpha',default_alpha);
            addParameter(inp,'gamma',default_gamma);
            addParameter(inp,'M_diag',default_M_diag);

            addParameter(inp,'excl_inds',defaultexcl_inds);
            addParameter(inp,'Aineq',defaultAineq);
            addParameter(inp,'bineq',defaultbineq);
            addParameter(inp,'opt_tol',defaultopt_tol);
            addParameter(inp,'const_tol',defaultconst_tol);
            addParameter(inp,'maxQPits',defaultmaxQPits);
            addParameter(inp,'dispQP',defaultdispQP);
            addParameter(inp,'auto',defaultauto);
            
            parse(inp,WS,varargin{:});  
            
            lambdas = inp.Results.lambdas;
            maxits = inp.Results.maxits;
            alpha = (inp.Results.alpha*mean(arrayfun(@(L)length(L.terms),WS.lib))+1)^-1;
            gamma = inp.Results.gamma;
            M_diag = inp.Results.M_diag;

            excl_inds = inp.Results.excl_inds;
            Aineq = inp.Results.Aineq;
            bineq = inp.Results.bineq;
            opt_tol = inp.Results.opt_tol;
            const_tol = inp.Results.const_tol;
            maxQPits = inp.Results.maxQPits;
            dispQP = inp.Results.dispQP;
            auto = inp.Results.auto;

            if isequal(auto,'weakLyap')
                Aineq = cell(3,1);
                bineq = cell(3,1);
                E = eye(3);
                for j=1:WS.numeq
                    tt = term('ftag',E(j,:),'linOp',1);
                    v = WS.tf{1}{j}.test(WS.dat,tt);
                    Aineq{j} = v(:).*WS.G{j};
                    bineq{j} = zeros(size(Aineq{j},1),1);
                end
            end

            wtemp = cell(WS.numeq,1);
            its = zeros(WS.numeq,1);
            loss_wsindy = zeros(WS.numeq+1,length(lambdas));
            for i=1:WS.numeq
                [wtemp{i},resid,its(i),lossvals,thrs_EL] = obj.wsindy_pde_RGLS_seq_qp(lambdas,gamma,G{i},b{i},M_diag{i},maxits,alpha,Aineq(i),bineq(i),excl_inds(i),opt_tol,const_tol,maxQPits,dispQP);
                loss_wsindy(i,:) = lossvals(1,:);
            end
            loss_wsindy(end,:) = lambdas;
            WS.weights = cell2mat(wtemp);
        end

        function [w,its,thrs_EL] = sparsifyDynamics(obj,G,b,lambda,gamma,M,maxits,toggle_jointthresh,linregargs,incl_inds)

            [~,nn] =size(G);
            n = size(b,2);
            if isempty(M)
                M = ones(nn,1);
            end
            if isequal(incl_inds,'all')
                incl_inds = 1:nn;
            end
            if  gamma ~= 0
                G = [G;gamma*eye(nn)];
                b = [b;zeros(nn,n)];
            end
            
            w = M.*obj.linreg(G,b,linregargs{:});
            if toggle_jointthresh == 1
                % threshold based on JCP paper
                bnds = norm(b)./vecnorm(G)'.*M;
                LBs = lambda*max(1,bnds);
                UBs = 1/lambda*min(1,bnds);
            elseif toggle_jointthresh == 2
                % threshold only on term magnitude
                bnds = norm(b)./vecnorm(G)'.*M;
                LBs = lambda*bnds;
                UBs = 1/lambda*bnds;
            elseif toggle_jointthresh == 3
                % threshold based on JCP but with term projection
                bnds = norm(b)^2./abs(b'*G)'.*M;
                bnds2 = norm(b)./vecnorm(G)'.*M;
                UBs = 1/lambda*bnds2; % upper bound by term magnitude
                LBs = lambda*bnds; % lower bound by projection
            elseif toggle_jointthresh == 4
                % threshold only on term projection
                bnds = norm(b)^2./abs(b'*G)'.*M;
                LBs = lambda*bnds;
                UBs = 1/lambda*bnds;
            else
                % threshold only Hamiltonian coarse-graining - should be
                % robust to small coefficients
                bnds = norm(b)./vecnorm(G)'.*M;
                w0 = abs(b'*G);
                nrms = vecnorm(G);
                alpha = max(w0./nrms.^2.*M');
                beta = max(w0./nrms/norm(b));
                LBs = lambda*max(alpha,bnds*beta); 
                UBs = 1/lambda*min(alpha,bnds*beta);
            end
            thrs_EL = [LBs bnds UBs];
            
            smallinds = 0*w;
            for j=1:min(nn,maxits)
                smallinds_new = or(abs(w)<LBs,abs(w)>UBs);
                smallinds_new(incl_inds) = 0;
                if all(smallinds_new(:)==smallinds(:))
                    its = j;
                    return
                else
                    smallinds = smallinds_new;
                    w(smallinds)=0;    
                    for ind=1:n
                        w(:,ind) = M.*obj.linreg(G(:,~smallinds),b(:,ind),linregargs{:},'S',~smallinds);
                    end
                end
            end
            its = j;
        end

        function [ws,its,thrs_ELs] = sparsifyGroupDynamics(obj,Gs,bs,lambda,gamma,Ms,maxits,toggle_jointthresh,linregargss,incl_inds,toggle_sign)
            %%% designed for single RHS vector only

            gs_norm = 1;

            ws = [];
            UBss = [];
            LBss = []; 
            for p = 1:length(Gs)
                G = Gs{p};
                b = bs{p};
                M = Ms{p};
                linregargs = linregargss{p};
        
                [~,nn] =size(G);
                if isempty(M)
                    M = ones(nn,1);
                end
                if isequal(incl_inds,'all')
                    incl_inds = 1:nn;
                end
                if  gamma ~= 0
                    G = [G;gamma*eye(nn)];
                    b = [b;zeros(nn,1)];
                end
                
                w = M.*obj.linreg(G,b,linregargs{:});
                if toggle_jointthresh == 1
                    % threshold based on JCP paper
                    bnds = norm(b)./vecnorm(G)'.*M;
                    LBs = lambda*max(1,bnds);
                    UBs = 1/lambda*min(1,bnds);
                elseif toggle_jointthresh == 2
                    % threshold only on term magnitude
                    bnds = norm(b)./vecnorm(G)'.*M;
                    LBs = lambda*bnds;
                    UBs = 1/lambda*bnds;
                elseif toggle_jointthresh == 3
                    % threshold based on JCP but with term projection
                    bnds = norm(b)^2./abs(b'*G)'.*M;
                    bnds2 = norm(b)./vecnorm(G)'.*M;
                    UBs = 1/lambda*bnds2; % upper bound by term magnitude
                    LBs = lambda*bnds; % lower bound by projection
                elseif toggle_jointthresh == 4
                    % threshold only on term projection
                    bnds = norm(b)^2./abs(b'*G)'.*M;
                    LBs = lambda*bnds;
                    UBs = 1/lambda*bnds;
                end
                ws = [ws w];
                Gs{p} = G;
                bs{p} = b;
                UBss = [UBss UBs];
                LBss = [LBss LBs];
            end

            UBs = vecnorm(UBss,gs_norm,2);
            LBs = vecnorm(LBss,gs_norm,2);
            thrs_ELs = [UBs LBs];

            w_comb = zeros(size(ws,1),1);
            smallinds = 0*w_comb ;
            for j=1:min(nn,maxits)

                %%% combine coeffs
                w_comb = vecnorm(ws,gs_norm,2);
                w_sign = abs(std(sign(ws),[],2));

                %%% threshold based on combined coeffs
                if ~toggle_sign
                    smallinds_new = or(w_comb<LBs,w_comb>UBs);
                elseif isequal(toggle_sign,true)
                    smallinds_new = any([w_comb<LBs,w_comb>UBs,w_sign],2);
                elseif isnumeric(toggle_sign)
                    if sum(~smallinds)<toggle_sign
                        smallinds_new = any([w_comb<LBs,w_comb>UBs,w_sign],2);
                    else
                        smallinds_new = or(w_comb<LBs,w_comb>UBs);
                    end
                end
                smallinds_new(incl_inds) = 0;
                if all(smallinds_new(:)==smallinds(:))
                    its = j;
                    return
                else
                    smallinds = smallinds_new;
                    ws(smallinds,:)=0;
                    for p = 1:length(Gs)
                        G = Gs{p};
                        b = bs{p};
                        M = Ms{p};
                        linregargs = linregargss{p};
                        ws(:,p) = M.*obj.linreg(G(:,~smallinds),b,linregargs{:},'S',~smallinds);
                    end
                end
            end
            its = j;

        end

        function [WS,its] = sparsifyDynamics_wendy(obj,WS,lambda,M,bnds,maxits,cov_thresh,vw)
            LBs = lambda*max(1./M,bnds);
            UBs = 1/lambda*min(1./M,bnds);
            smallinds = WS.weights*0;
            n = length(smallinds);
            for j=1:min(n,maxits)
                smallinds_new = or(abs(WS.weights)<LBs,abs(WS.weights)>UBs);
                if all(smallinds_new(:)==smallinds(:))
                    its = j;
                    return
                else
                    smallinds = smallinds_new;
                    w = WS.weights;
                    w(smallinds) = 0;
                    WS.add_weights(w,'toggle_cov',1);
                    if any(w)
                        [WS,~,~,~,C] = obj.wendy(WS,vw{:});
                        w = WS.weights;
                        inds = find(w);
                        if ~isempty(inds)
                            I = abs(w(inds)) < sqrt(diag(C))*cov_thresh;
                            w(inds(I)) = 0;
                            WS.add_weights(w,'toggle_cov',1);
                        end
                    end
                end
            end
            its = j;
        end

        function [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq_qp(obj,lambdas,gamma,G,b,M,maxits,alpha,A,c,excl_inds,opt_tol,const_tol,max_its,disp_opt)

            maxits=min(maxits,size(G,2));
            
            [~,m] = size(G);
            [~,num_eq] = size(b);
            
            W_ls = [G;gamma*eye(m)] \ [b;zeros(m,num_eq)];
            GW_ls = norm(G*W_ls);
            
            proj_cost = [];
            overfit_cost = [];
            lossvals = [];
            
            if isempty(lambdas)
                lam_max = max(max(abs(G'*b),[],2)./vecnorm(G).^2');
                lam_min = min(vecnorm(G*W_ls))/size(G,2)/max(vecnorm(G));
                lambdas = 10.^linspace(log10(lam_min), log10(lam_max),50);
            end
            
            if and(length(lambdas)==1,all(lambdas<0))
                lam_max = max(max(abs(G'*b),[],2)./vecnorm(G).^2');
                lam_min = min(vecnorm(G*W_ls))/size(G,2)/max(vecnorm(G));
                lambdas = 10.^linspace(log10(lam_min), log10(lam_max),-lambdas);
            end
           
            W = zeros(m,num_eq);
            for l=1:length(lambdas)
                lambda = lambdas(l);
                for k=1:num_eq
                    if isempty(M)
                        [W(:,k),its,~] = obj.sparsifyDynamics_qp(G,b(:,k),lambda,gamma,[],A{k},c{k},find(excl_inds{k}),opt_tol,const_tol,max_its,disp_opt,maxits);
                    else
                        [W(:,k),its,~] = obj.sparsifyDynamics_qp(G,b(:,k),lambda,gamma,M(:,k),A{k},c{k},find(excl_inds{k}),opt_tol,const_tol,max_its,disp_opt,maxits);
                        W(:,k) = W(:,k)./M(:,k);
                    end
                end
                proj_cost = [proj_cost alpha*norm(G*(W-W_ls))/GW_ls];
                overfit_cost = [overfit_cost (1-alpha)*length(find(W))/length(find(W_ls))];
                lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
            end
            
            l = find(lossvals == min(lossvals),1);
            lambda = lambdas(l);
            its_all = zeros(num_eq,1);
            
            resid = b*0;
            for k=1:num_eq
                if isempty(M)
                    [W(:,k),its,thrs_EL] = obj.sparsifyDynamics_qp(G,b(:,k),lambda,gamma,[],A{k},c{k},find(excl_inds{k}),opt_tol,const_tol,max_its,disp_opt,maxits);
                    resid(:,k) = (b(:,k) - G*W(:,k))/norm(b(:,k)); 
                else
                    [W(:,k),its,thrs_EL] = obj.sparsifyDynamics_qp(G,b(:,k),lambda,gamma,M(:,k),A{k},c{k},find(excl_inds{k}),opt_tol,const_tol,max_its,disp_opt,maxits);
                    resid(:,k) = (b(:,k) - G*(W(:,k)./M(:,k)))/norm(b(:,k)); 
                end
                its_all(k) = its;
            end
            lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
        end
        
        function [Xi,its,thrs_EL] = sparsifyDynamics_qp(obj,Theta,dXdt,lambda,gamma,M,A,b,excl_inds,opt_tol,const_tol,max_its,disp_opt,max_its_stls)
        % Copyright 2015, All Rights Reserved
        % Code by Steven L. Brunton
        % For Paper, "Discovering Governing Equations from Data: 
        %        Sparse Identification of Nonlinear Dynamical Systems"
        % by S. L. Brunton, J. L. Proctor, and J. N. Kutz
        %
        % modified by Daniel A. Messenger, 2020 to prevent return of zero vector
        % and include regularization
        %
        % compute Sparse regression: sequential least squares
            if isempty(disp_opt)
                disp_opt='none';
            end
            if isempty(const_tol)
                const_tol=10^-16;
            end
            options = optimoptions('quadprog','Display',disp_opt,'ConstraintTolerance',const_tol,'OptimalityTolerance',opt_tol,'MaxIterations',max_its);
            n = min(size(dXdt));
            nn = size(Theta,2);
        
            dXdt_conj = Theta'*dXdt;
            Theta_conj = Theta'*Theta;
        
            if  gamma ~= 0
                Theta_conj = Theta_conj+gamma^2*eye(nn);
            end
        
            Xi = quadprog(Theta_conj,-dXdt_conj,A,b,[],[],[],[],[],options);  % initial guess: Least-squares
            if isempty(M)
                thrs_EL = [];
            else
                Xi = M.*Xi;
                bnds = norm(dXdt)./vecnorm(Theta)'.*M; 
                LBs = lambda*max(1,bnds);
                UBs = 1/lambda*min(1,bnds);
                thrs_EL = [LBs bnds UBs];
            end
        
            smallinds = 0*Xi;
            its = 0;
            while its < max_its_stls
                if ~isempty(M)
                    smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
                    smallinds_new(excl_inds) = 0;
                    if or(length(find(smallinds_new))==nn,all(smallinds_new(:)==smallinds(:)))
                        return
                    else
                        smallinds = smallinds_new;
                        Xi(smallinds)=0;
                        for ind=1:n
                            biginds = ~smallinds(:,ind);
                            Xi(biginds,ind) = M(biginds).*quadprog(Theta_conj(biginds,biginds),-dXdt_conj(biginds,ind),A(:,biginds(1:min(size(A,2),end))),b,[],[],[],[],[],options);
                        end
                    end
                else
                    smallinds_new = (abs(Xi)<lambda);
                    smallinds_new(excl_inds) = 0;
                    if or(all(smallinds_new(:)==smallinds(:)),length(find(smallinds_new))==length(Xi))
                        its = j;
                        return
                    else
                        smallinds = smallinds_new;
                        Xi(smallinds)=0;
                        for ind = 1:n        
                            biginds = ~smallinds(:,ind);
                            Xi(biginds,ind) = quadprog(Theta_conj(biginds,biginds),-dXdt_conj(biginds,ind),A(:,biginds(1:min(size(A,2),end))),b,[],[],[],[],[],options);
                        end
                    end
                end
                its=its+1;
            end
        end
    
        function [WS,CV_all,supports_all] = subspacePursuitCV(obj,WS,varargin)

            inp = inputParser;
            addParameter(inp,'s',15);
            addParameter(inp,'ncv',25);
            addParameter(inp,'toggle_discrep',0);
            addParameter(inp,'M_diag',[]);
            addParameter(inp,'linregargs',{});

            parse(inp,varargin{:});  
            
            s = inp.Results.s;
            ncv = inp.Results.ncv;
            toggle_discrep = inp.Results.toggle_discrep;
            M_diag = inp.Results.M_diag;
            linregargs = inp.Results.linregargs;

            if isequal(class(WS),'wsindy_model')
                WS.cat_Gb;            
                G = WS.G; b = WS.b;
                if and(toggle_discrep==1,~isempty(WS.weights))
                    b = cellfun(@(g,b) b-g*WS.weights,G,b,'uni',0);
                    G = cellfun(@(g) g(:,~WS.weights),G,'uni',0);
                elseif and(toggle_discrep==2,~isempty(WS.weights))
                    b = cellfun(@(g,b) b-g*WS.weights,G,b,'uni',0);
                end 
            elseif isequal(class(WS),'cell')
                G = WS{1};
                b = WS{2};
                toggle_discrep = 0;
            end

            if isempty(linregargs)
                linregargs = cellfun(@(g)[],G,'un',0);
            end

            wtemp = cell(length(b),1);
            supports_all = cell(length(b),1);
            CV_all = cell(length(b),1);
            for j=1:length(b)
                supportList            = cell(min(s,size(G{j},2)),1);
                crossValidationErrList = zeros(min(s,size(G{j},2)),1);
                WColumnNorm            = vecnorm(G{j});
            
                for i=1:min(s,size(G{j},2))
                    support       = obj.SPV2(G{j} * diag(1./WColumnNorm), b{j} ./norm(b{j},2), i, linregargs{j});
                    supportList{i}   = support;
                    crossValidationErrList(i) = obj.computeCrossValidationErrV2_dam(support, G{j}, b{j}, ncv);
                end
                [~, CrossIdx]=min(crossValidationErrList);
                supportPred = supportList{CrossIdx}';
    
                S = logical(zeros(size(G{j},2),1));
                S(supportPred) = 1;
                wtemp{j} = obj.linreg(G{j}(:,supportPred),b{j},linregargs{j}{:},'S',S);
                % if isequal(class(WS),'wsindy_model') 
                %     if and(toggle_discrep==1,~isempty(WS.weights))
                %         wtemp{j} = WS.reshape_w{j};
                %         wtemp2 = zeros(size(G{j},2),1);
                %         wtemp2(supportPred) = G{j}(:,supportPred) \ b{j};
                %         wtemp(~wtemp) = wtemp2;
                %     elseif and(toggle_discrep==2,~isempty(WS.weights))
                %         wtemp{j} = WS.reshape_w{j};
                %         wtemp{j}(supportPred) = wtemp{j}(supportPred) + G{j}(:,supportPred) \ b{j};
                %     else
                %         wtemp{j} = zeros(size(G{j},2),1);
                %         wtemp{j}(supportPred) = G{j}(:,supportPred) \ b{j};
                %     end
                % elseif isequal(class(WS),'cell')
                %     wtemp{j} = zeros(size(G{j},2),1);
                %     wtemp{j}(supportPred) = G{j}(:,supportPred) \ b{j};
                % end
                supports_all{j} = supportList;
                CV_all{j} = crossValidationErrList;
            end
            if ~isempty(M_diag)
                wtemp = cellfun(@(w,M) M.*w, wtemp,M_diag,'un',0);
            end
            if isequal(class(WS),'wsindy_model')
                WS.weights = cell2mat(wtemp);
            elseif isequal(class(WS),'cell')
                WS = wtemp(:)';
            end
        end

        function support = SPV2(obj,W,b,sparsity,linregargs)
            % "Subspace Pursuit for Compressive Sensing: Closing the
            %  Gap Between Performance and Complexity"%
            
            % INPUT
            % W        :  feature matrix
            % sparsity :  sparsity level
            
            % OUTPUT
            % support  :  a list of index
            
            itermax = 15;
            [~,N]=size(W);
            
            cv = abs( b'*W );
            [~, cv_index] = sort(cv,'descend');
            
            Lda = cv_index(1:sparsity);
            Phi_Lda = W(:,Lda);
            
            % x = (Phi_Lda'*Phi_Lda)\(Phi_Lda' * b);

            S = logical(zeros(N,1));
            S(Lda) = 1;
            x = obj.linreg(Phi_Lda,b,linregargs{:},'S',S);

            r = b - Phi_Lda*x(Lda);
            res = norm(r);
            iter = 0;
            
            if (res < 1e-12)
                X = x;
                support = find(X);
                return
            end
            
            usedlda = zeros(1,N);
            usedlda(Lda)=1;
            
            for iter = 1:itermax
                res_old = res;
                %%Step1 find T^{\prime} and add it to \hat{T}
                cv = abs( r'*W );
                [~, cv_index] = sort(cv,'descend');
                Sga = union(Lda, cv_index(1:sparsity));
                Phi_Sga = W(:,Sga);
                
                %%% find the most significant K indices
                S = logical(zeros(N,1));
                S(Sga) = 1;
                x_temp = obj.linreg(Phi_Sga,b,linregargs{:},'S',S);
                x_temp = x_temp(S);
                % x_temp = (Phi_Sga'*Phi_Sga)\(Phi_Sga' * b);        
                [~, x_temp_index] = sort( abs(x_temp) , 'descend' );
                Lda = Sga(x_temp_index(1:sparsity));
                Phi_Lda = W(:,Lda);
                usedlda(Lda)=1;
                
                %%% calculate the residue
                S = logical(zeros(N,1));
                S(Lda) = 1;
                x = obj.linreg(Phi_Lda,b,linregargs{:},'S',S);
                % x = (Phi_Lda'*Phi_Lda)\(Phi_Lda' * b);
                r = b - Phi_Lda*x(Lda);
                res = norm(r);
            
                X=x;
                if ( res/res_old >= 1 || res < 1e-12)
                    support = find(X);
                    return
                end
            end
            support = find(X);
        
        end

        function [err] = computeCrossValidationErrV2_dam(obj,ind, A, b, ntrials)
            % script for computing cross validation error
            errCrossAccu = zeros(1,ntrials);
            for jj = 1:ntrials
                e=obj.computeCVErrV4(ind,A, b);
                errCrossAccu(jj) =  e;
            end
            err = mean(errCrossAccu)  + std(errCrossAccu);
            end
            
        function [err]=computeCVErrV4(obj,support,W,b,ratio)
            if ~exist('ratio','var')
                ratio = 1/100;
            end
            
            n     = size(b,1);
            inds  = randperm(n);
            W     = W(inds,:);
            b     = b(inds);
            
            % split the data into two parts
            endOfPart1 = floor(n*ratio);
            
            IdxPart1  = 1:endOfPart1;
            IdxPart2  = endOfPart1+1:n;
            
            
            % compute e1
            coeff            = zeros(size(W,2),1);
            coeff(support)   = W(IdxPart1,support)\b(IdxPart1);
            e1               = norm(W(IdxPart2,:)*coeff - b(IdxPart2) ,2)/norm(b(IdxPart2),2);
            
            % compute e2
            coeff            = zeros(size(W,2),1);
            coeff(support)   = W(IdxPart2,support)\b(IdxPart2);
            e2               = norm(W(IdxPart1,:)*coeff - b(IdxPart1) ,2)/norm(b(IdxPart1),2);
            
            % compute weighted error
            err              = e1 * (1-ratio) + e2 * ratio;
            
        end

    end

end