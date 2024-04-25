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

        function [WS,loss_wsindy,its,G,b] = MSTLS(obj,WS,varargin)
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
            alpha = (p.Results.alpha*mean(arrayfun(@(L)length(L.terms),WS.lib))+1)^-1;
            gamma = p.Results.gamma;
            M_diag = p.Results.M_diag;
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
                    if isempty(Aeq)
                        [W_all(:,l),~] = sparsifyDynamics(G{k}, b{k}, lambdas(l), 1, gamma, M_diag{k}, maxits, toggle_jointthresh, reg_inds{k});
                    else
                        [W_all(:,l),~] = sparsifyDynamics_lineq(G{k}, b{k}, lambdas(l), 1, gamma, M_diag{k}, maxits, toggle_jointthresh, reg_inds{k},Aeq{k},deq{k});
                    end
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
            Aeq = p.Results.Aeq;
            beq = p.Results.beq;
            LB = p.Results.LB;
            UB = p.Results.UB;
            consttol = p.Results.consttol;
            opttol = p.Results.opttol;
            maxits = p.Results.maxits;
            S = p.Results.S;
            verbosity = p.Results.verbose;

            if any(S)
                if isempty(x0)
                    if diff(size(A))>=0
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
                    options = optimoptions('quadprog','Display',verbosity,'ConstraintTolerance',consttol,'OptimalityTolerance',opttol,'MaxIterations',maxits);
                    x = quadprog((A'*A),-(A'*b),Aineq,bineq,Aeq,beq,LB,UB,x0,options);
                    x = obj.inject_sparse(x,S);
                else
                    x = x0;
                end
            else
                x = S*0;
            end

        end

        function [WS,w_its,res,res_0,CovW] = wendy(obj,WS,varargin)
            % options: maxits,ittol,diag_reg,w,regmeth

            default_maxits = 20;
            default_ittol = 10^-6;
            default_diag_reg = 10^-6;
            default_w = WS.weights;
            default_regmeth = 'ols';

            p = inputParser;
            addRequired(p,'WS');
            addParameter(p,'maxits',default_maxits);
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

            if verbosity
                tic,
            end

            if ~isempty(w)
                if isequal(w,0)
                    WS = obj.ols(WS,'S',0,'linregargs',linregargs);                    
                else
                    WS.add_weights(w);
                end
            else
                WS = obj.ols(WS,'linregargs',linregargs);
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
                    WS.add_weights(w);
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
            if exist('RT','var')
                Ginv = Ginv*RT;
            end
            CovW = Ginv*Ginv';
            if verbosity
                disp(['wendy iter time=',num2str(toc),'; sparsity=',num2str(length(find(WS.weights))),'; its=',num2str(its)])
            end
        end

        function [WS,loss_wsindy,lambda,w_its,res,res_0,CovW] = MSTLS_WENDy(obj,WS,varargin)

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
            WS.add_weights(W_ls);
            
            proj_cost = []; overfit_cost = []; lossvals = [];
            
            Wmat = zeros(length(W_ls),length(lambdas));
            for l=1:length(lambdas)
                WS.weights = W_ls;
                [WS,its] = obj.sparsifyDynamics_wendy(WS,lambdas(l),M_diag,bnds,maxits,vw);
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

            [WS,w_its,res,res_0,CovW] = obj.wendy(WS,vw{:});

            WS.add_weights(WS.weights.*M_diag);
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

        function WS = subspacePursuitCV(obj,WS,s,ncv,toggle_discrep)
            if ~exist('toggle_discrep','var')
                toggle_discrep = 0;
            end

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

            wtemp = cell(length(b),1);
            for j=1:length(b)
                supportList            = cell(min(s,size(G{j},2)),1);
                crossValidationErrList = zeros(min(s,size(G{j},2)),1);
                WColumnNorm            = vecnorm(G{j});
            
                for i=1:min(s,size(G{j},2))
                    support       = SPV2(G{j} * diag(1./WColumnNorm), b{j} ./norm(b{j},2),i);
                    supportList{i}   = support;
                    crossValidationErrList(i) = computeCrossValidationErrV2_dam(support, G{j}, b{j}, ncv);
                end
                [~, CrossIdx]=min(crossValidationErrList);
                supportPred = supportList{CrossIdx}';
    
                if isequal(class(WS),'wsindy_model') 
                    if and(toggle_discrep==1,~isempty(WS.weights))
                        wtemp{j} = WS.reshape_w{j};
                        wtemp2 = zeros(size(G{j},2),1);
                        wtemp2(supportPred) = G{j}(:,supportPred) \ b{j};
                        wtemp(~wtemp) = wtemp2;
                    elseif and(toggle_discrep==2,~isempty(WS.weights))
                        wtemp{j} = WS.reshape_w{j};
                        wtemp{j}(supportPred) = wtemp{j}(supportPred) + G{j}(:,supportPred) \ b{j};
                    else
                        wtemp{j} = zeros(size(G{j},2),1);
                        wtemp{j}(supportPred) = G{j}(:,supportPred) \ b{j};
                    end
                elseif isequal(class(WS),'cell')
                    wtemp{j} = zeros(size(G{j},2),1);
                    wtemp{j}(supportPred) = G{j}(:,supportPred) \ b{j};
                end
            end
            if isequal(class(WS),'wsindy_model')
                WS.weights = cell2mat(wtemp);
            elseif isequal(class(WS),'cell')
                WS = wtemp(:)';
            end

        end

        function [WS,its] = sparsifyDynamics_wendy(obj,WS,lambda,M,bnds,maxits,vw)
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
                    WS.add_weights(w);
                    if any(w)
                        [WS,~,~,~,C] = obj.wendy(WS,vw{:});
                        % w = WS.weights;
                        % inds = find(w);
                        % if ~isempty(inds)
                        %     I = abs(w(inds)) < sqrt(diag(C))/4;
                        %     w(inds(I)) = 0;
                        %     WS.add_weights(w);
                        % end
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

    end

end

