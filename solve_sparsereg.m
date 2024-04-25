function [W,resid,its_all,lossvals,thrs_EL,lambda_hat,G_fin,b_fin] = solve_sparsereg(G,b,varargin)
    defaultmeth = 'MSTLS';
    defaultkeeprows = 1:size(G,1);
    defaultcov = ones(size(G,1),1);
    defaultlambda = 10.^(linspace(-4,0,50));
    defaultgamma = 0;
    defaultM = [];
    defaultmaxits = size(G,2);
    defaultsparsity_scale=0;
    defaultalpha=0.5;
    defaultexcl_inds=repmat({[]},1,size(b,2));
    defaultAineq=repmat({[]},1,size(b,2));
    defaultbineq=repmat({[]},1,size(b,2));
    defaulttol=[];
    defaultmaxQPits = 1000;
    defaultdispQP='off';
    defaultmuDR=1;
    
    inp = inputParser;
    addParameter(inp,'meth',defaultmeth);
    addParameter(inp,'keeprows',defaultkeeprows);
    addParameter(inp,'cov',defaultcov);
    addParameter(inp,'lambda',defaultlambda);
    addParameter(inp,'gamma',defaultgamma);
    addParameter(inp,'M',defaultM);
    addParameter(inp,'maxits',defaultmaxits);
    addParameter(inp,'sparsity_scale',defaultsparsity_scale);
    addParameter(inp,'alpha',defaultalpha);
    addParameter(inp,'excl_inds',defaultexcl_inds);
    addParameter(inp,'Aineq',defaultAineq);
    addParameter(inp,'bineq',defaultbineq);
    addParameter(inp,'tol',defaulttol);
    addParameter(inp,'maxQPits',defaultmaxQPits);
    addParameter(inp,'dispQP',defaultdispQP);
    addParameter(inp,'muDR',defaultmuDR);
    
    parse(inp,varargin{:});  
    
    meth = inp.Results.meth;
    inds_keep = inp.Results.keeprows;
    cov = inp.Results.cov;
    lambda = inp.Results.lambda;
    gamma = inp.Results.gamma;
    M = inp.Results.M;
    maxits = inp.Results.maxits;
    sparsity_scale = inp.Results.sparsity_scale;
    alpha = inp.Results.alpha;
    excl_inds = inp.Results.excl_inds;
    Aineq = inp.Results.Aineq;
    bineq = inp.Results.bineq;
    tol = inp.Results.tol;
    maxQPits = inp.Results.maxQPits;
    muDR = inp.Results.muDR;
    dispQP = inp.Results.dispQP;
    
    G_fin = (1./cov(inds_keep)).*G(inds_keep,:);
    b_fin = (1./cov(inds_keep)).*b(inds_keep,:);
    
    if isequal(meth,'MSTLS')
        if and(length(lambda)==1,all(lambda>=0))
            [W,resid,its_all,thrs_EL] = wsindy_pde_RGLS(lambda,gamma,G_fin,b_fin,M,maxits);
            lambda_hat = lambda;
            lossvals = [];
        else
            if isequal(sparsity_scale,0)
                [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq(lambda,gamma,G_fin,b_fin,M,maxits,alpha);
                lambda_hat_ind=find(lossvals(min(end,4),:)>0,1,'last');
                lambda_hat = lossvals(min(end,4),lambda_hat_ind);
            elseif isequal(sparsity_scale,1)
                [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq2(lambda,gamma,G_fin,b_fin,M,maxits,alpha);
                lambda_hat_ind=find(lossvals(min(end,4),:)>0,1,'last');
                lambda_hat = lossvals(min(end,4),lambda_hat_ind);
            else
                [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq(lambda,gamma,G_fin,b_fin,M,maxits,alpha);
                lambda_hat_ind=find(lossvals(min(end,4),:)>0,1,'last');
                loss_min = lossvals(1,lambda_hat_ind);
                [W1,resid1,its_all1,lossvals1,thrs_EL1] = wsindy_pde_RGLS_seq2(lambda,gamma,G_fin,b_fin,M,maxits,alpha);
                lambda_hat_ind1=find(lossvals1(min(end,4),:)>0,1,'last');
                loss_min1 = lossvals1(1,lambda_hat_ind1);
                if loss_min1<loss_min
                    W=W1;
                    resid=resid1;
                    its_all=its_all1;
                    lossvals=lossvals1;
                    thrs_EL=thrs_EL1;
                    lambda_hat_ind=lambda_hat_ind1;
                end
                lambda_hat = lossvals(min(end,4),lambda_hat_ind);
            end
        end
        if gamma>0
            for i=1:size(b_fin,2)
                if ~isempty(M)
                    W(W(:,i)~=0,i) = M(W(:,i)~=0,i).*(G_fin(:,W(:,i)~=0) \ b_fin(:,i));
                else
                    W(W(:,i)~=0,i) = G_fin(:,W(:,i)~=0) \ b_fin(:,i);
                end
            end
        end
    elseif isequal(meth,'MSTLSQP')
        [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq_qp(lambda,gamma,G_fin,b_fin,M,maxits,alpha,Aineq,bineq,excl_inds,tol,maxQPits,dispQP);
        lambda_hat = lossvals(min(end,4),lossvals(min(end,4),:)>0);
        if ~isempty(lambda_hat)
            lambda_hat = lambda_hat(end);
        else
            lambda_hat=lambda(1);
        end
    elseif isequal(meth,'DR')
        [W,its_all] = dougrach(G_fin,b_fin,gamma,lambda,muDR,maxits,tol,alpha,M);
        if ~isempty(M)
            W = W.*M;
        end
        resid=[];its_all=1;
        lossvals=[];
        thrs_EL=[];lambda_hat=[];
    elseif isequal(meth,'LOOCV')% leave one out cross val
        runs=gamma;trainfrac=lambda;
        [W,routs,inds_keep] = sparsesingletermthresh(G_fin,b_fin,trainfrac,runs,tol);
        if ~isempty(M)
            W = W.*M;
        end
        resid=[];its_all=1;
        lossvals=[routs';10.^(1:length(routs))];
        thrs_EL=[];lambda_hat=[];
    elseif isequal(meth,'QR')
        [W,Wcell,res,~,~] = qr_sparse_reg(G_fin,b_fin,M,lambda,gamma,tol);
        resid=[];its_all=1;
        lossvals=[res;10.^(1:length(res))];
        thrs_EL=[];lambda_hat=[];
    %     figure(10)
    %     subplot(1,2,1)
    %     plot(res,'o-')
    %     subplot(1,2,2)
    %     spy(Wcell{1})
    elseif isequal(meth,'MSTLSBF')
        W=stlsbf(G,b,M,alpha,muDR,tol,maxQPits,maxits,lambda,gamma);
        resid=[];its_all=1;
        lossvals=[];
        thrs_EL=[];lambda_hat=[];
    elseif isequal(meth,'MIO')
        W=MIO(G,b,M,lambda,tol,alpha,gamma,muDR);
        %Ws=MIO(G,b,M,k_start,logtol,restol,eps1,eps2)
        its_all=1;
        lossvals=[];
        thrs_EL=[];
        lambda_hat=[];
    elseif isequal(meth,'SP')
        W = subspacePursuitCV(G,b,M,lambda,maxits);
        its_all=1;
        lossvals=[];
        thrs_EL=[];
        lambda_hat=[];
    end
    
    if ~isempty(M)
        resid = (G_fin*(W./M)-b_fin)./vecnorm(b_fin); 
    else
        resid = (G_fin*W-b_fin)./vecnorm(b_fin); 
    end

end

function Wfin=stlsbf(G,b,M,epinc,max_rr_fac,bflim,maxits,max_its_stls,lambda_init,gamma)
    if epinc==0
        epinc=0.99;
    end
    if max_rr_fac<1
        max_rr_fac=1;
    end
    lambdainc=10^(4/100);
    num_eq = size(b,2);
    Wfin=zeros(size(G,2),num_eq);
    epfin=zeros(num_eq,1);
    for n=1:num_eq
        i=1;
        S=size(G,2);
        choices=inf;
        while S>bflim
            W=sparsifyDynamics(G,b(:,n),lambda_init,gamma,M(:,n),max_its_stls);
            supp=find(W~=0);
            S=length(supp);
            i=i+1;
            lambda_init=lambda_init*lambdainc;
        end

        Gred=G(:,supp)./M(supp,n)';
        [Q,R]=qr(Gred,0);
        bred=Q'*b(:,n);
        inds=logical(de2bi(1:2^S-1));
        Ws= zeros(S,2^S-1);
        for j=1:2^S-1
            Ws(inds(j,:),j)=R(:,inds(j,:))\bred;
        end
        rr=vecnorm(Gred*Ws-b(:,n))/norm(b(:,n));
        max_ep=rr(end)*max_rr_fac;
        ep=max_ep;
        ii=1;
        while and(choices>1,ii<maxits)
            inds_keep=find(rr<ep);
            inds_fin=inds_keep(sum(inds(inds_keep,:),2)==min(sum(inds(inds_keep,:),2)));
            choices=length(inds_fin);
            if choices==0
                ep=ep/epinc;
                epinc=epinc^2;
                choices=2;
            else
                ep=ep*epinc;
            end
            ii=ii+1;
        end
        if ii==maxits
            disp('reached maximum iteration')
            inds_keep=find(rr<=max_ep);
            inds_fin=inds_keep(sum(inds(inds_keep,:),2)==min(sum(inds(inds_keep,:),2)));
            epfin(n)=max_ep;
            Wfin(supp,n)=Ws(:,inds_fin(:,1));
        else            
            epfin(n)=ep/epinc;
            Wfin(supp,n)=Ws(:,inds_fin);
        end
    end
end

function [W,resid,its_all,thrs_EL] = wsindy_pde_RGLS(lambda,gamma,G,b,M,maxits)
    [~,J] = size(G);
    [~,num_eq] = size(b);
    W = zeros(J,num_eq);
    
    its_all = zeros(num_eq,1);
    resid = b*0;
    for k=1:num_eq
        if isempty(M)
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M,maxits);
            resid(:,k) = (b(:,k) - G*W(:,k))/norm(b(:,k)); 
        else
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,k), maxits);
            resid(:,k) = (b(:,k) - G*(W(:,k)./M(:,k)))/norm(b(:,k)); 
        end
        its_all(k) = its;
    end
end

function [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq(lambdas,gamma,G,b,M,maxits,alpha)

    [~,m] = size(G);
    [~,num_eq] = size(b);
    
    if or(isequal(gamma,0),isempty(gamma))
        W_ls = G \ b;
    else    
        W_ls = [G;gamma*eye(m)] \ [b;zeros(m,num_eq)];
    end
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
                [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, [], maxits);
            else
                [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,k), maxits);
                W(:,k) = W(:,k)./M(:,k);
            end
        end
        proj_cost = [proj_cost 2*alpha*norm(G*(W-W_ls))/GW_ls];
        overfit_cost = [overfit_cost 2*(1-alpha)*length(find(W~=0))/length(find(W_ls~=0))];
        lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
    end
    
    l = find(lossvals == min(lossvals),1);
    
    lambda = lambdas(l);
    its_all = zeros(num_eq,1);
    
    resid = b*0;
    for k=1:num_eq
        if isempty(M)
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, []);
            resid(:,k) = (b(:,k) - G*W(:,k))/norm(b(:,k)); 
        else
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,k), maxits);
            resid(:,k) = (b(:,k) - G*(W(:,k)./M(:,k)))/norm(b(:,k)); 
        end
        its_all(k) = its;
    end
    lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
end

function [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq2(lambdas,gamma,G,b,M,maxits,alpha)

    [~,m] = size(G);
    [~,num_eq] = size(b);
        
    if or(isequal(gamma,0),isempty(gamma))
        W_ls = G \ b;
    else    
        W_ls = [G;gamma*eye(m)] \ [b;zeros(m,num_eq)];
    end
    GW_ls = norm(G*W_ls);
    
    if length(lambdas)<=1
    %     if isempty(lambdas)
    %         lambdas = 100;
    %     end
    %     Wtemp = sort(abs(G'*b./vecnorm(G).^2'));
    %     ind = 3;%findchangepts(log10(Wtemp));
    %     lambdas = 10.^linspace(log10(mean(Wtemp))-3, log10(mean(Wtemp)),lambdas);
        if isempty(lambdas)
            num_lam = 100;
        else
            num_lam = -lambdas;
        end
        lam_max = min(max(max(abs(G'*b),[],2)./vecnorm(G).^2'),1);
        lam_min = min(vecnorm(G*W_ls))/size(G,2)/max(vecnorm(G));
        lambdas = 10.^linspace(log10(lam_min), log10(lam_max),num_lam);
    end
    
    proj_cost = [];
    overfit_cost = [];
    lossvals = [];
    
    W = zeros(m,num_eq);
    
    for l=1:length(lambdas)
        lambda = lambdas(l);
        for k=1:num_eq
            [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, ones(m,1), maxits);
        end    
        proj_cost = [proj_cost 2*alpha*norm(G*(W-W_ls))/GW_ls];
        overfit_cost = [overfit_cost 2*(1-alpha)*length(find(W~=0))/length(W(:))];
        lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
    end
    
    l = find(lossvals == min(lossvals),1);
    
    lambda = lambdas(l);
    its_all = zeros(num_eq,1);
    
    resid = b*0;
    for k=1:num_eq
        [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, ones(m,1));
        resid(:,k) = (b(:,k)-G*W(:,k))/norm(b(:,k));
        if ~isempty(M)
            W(:,k) = W(:,k).*M(:,k);
        end
        its_all(k) = its;
    end
    lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
end

function [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq_qp(lambdas,gamma,G,b,M,maxits,alpha,A,c,excl_inds,const_tol,max_its,disp_opt)

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
                [W(:,k),~,~] = sparsifyDynamics_qp(G,b(:,k),lambda,gamma,[],A{k},c{k},find(excl_inds{k}),const_tol,max_its,disp_opt,maxits);
            else
                [W(:,k),~,~] = sparsifyDynamics_qp(G,b(:,k),lambda,gamma,M(:,k),A{k},c{k},find(excl_inds{k}),const_tol,max_its,disp_opt,maxits);
                W(:,k) = W(:,k)./M(:,k);
            end
        end
        proj_cost = [proj_cost 2*alpha*norm(G*(W-W_ls))/GW_ls];
        overfit_cost = [overfit_cost 2*(1-alpha)*length(find(W~=0))/length(find(W_ls~=0))];
        lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
    end
    
    l = find(lossvals == min(lossvals),1);
    lambda = lambdas(l);
    its_all = zeros(num_eq,1);
    
    resid = b*0;
    for k=1:num_eq
        if isempty(M)
            [W(:,k),its,thrs_EL] = sparsifyDynamics_qp(G,b(:,k),lambda,gamma,[],A{k},c{k},find(excl_inds{k}),const_tol,max_its,disp_opt,maxits);
            resid(:,k) = (b(:,k) - G*W(:,k))/norm(b(:,k)); 
        else
            [W(:,k),its,thrs_EL] = sparsifyDynamics_qp(G,b(:,k),lambda,gamma,M(:,k),A{k},c{k},find(excl_inds{k}),const_tol,max_its,disp_opt,maxits);
            resid(:,k) = (b(:,k) - G*(W(:,k)./M(:,k)))/norm(b(:,k)); 
        end
        its_all(k) = its;
    end
    lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
end

function [w,num_its] = dougrach(A,b,gamma,lambda,mu,maxits,tol,alph,M)

    normz = vecnorm(A,inf);
    A = A./normz;
    [m,n] = size(A);
    w = (A \ b);%+randn(size(A,2),1);
    
    num_its=1;
    check=2*tol;
    
    if isempty(M)
        M = w*0+1;
    end
    
    while and(num_its < maxits, check>tol)
	    wprox = max(abs(w)-gamma*lambda,0).*sign(w);
	    wprox = 2*wprox-w;
	    wstar = ([sqrt(gamma)*A;eye(n)]) \ [sqrt(gamma)*b; wprox];
	    wstar = (1-mu/2)*w+mu/2*(2*wstar-wprox);
    % 	wstar(M.*abs(wstar)./normz'<alph)=0;
	    check = norm(wstar-w)/norm(w);
	    w=wstar;
	    num_its = num_its+1;
    end 
    
    for i=1:size(b,2)
        for k=1:2
            smallinds = abs(w(:,i)./normz'.*M(:,i))<alph;
            w(smallinds,i) = 0;
            inds = w(:,i)~=0;
            w(inds,i) = (A(:,inds) \ b(:,i))./normz(inds)';
        end
    end

end

function [wsparse,Wcell,res,Q,R] = qr_sparse_reg(G,b,M,lambda,gamma,tol)
    [Q,R] = qr(G,0);
    wsparse = zeros(size(G,2),size(b,2));
    Wcell = repmat({R*0},size(b,2),1);
    for nn=1:size(b,2)
        w = Q \ b(:,nn);
        [~,inds] = sort(abs(w),'descend');
        W = [];
        Wbool = w*0;
        res = [];
        restemp = 1;
        i = 1;
        check = tol+1;
        while and(check>tol,i<=size(G,2))
            W = [W w*0];
            res = [res restemp];
            proj = Q(:,inds(i))'*G(:,1:inds(i))./vecnorm(G(:,1:inds(i))).^2;
%             proj = R(inds(i),1:inds(i))./vecnorm(G(:,1:inds(i))).^2;
            [~,a] = max(abs(proj));
            Wbool(:,i) = Wbool(:,end);
            Wbool(a,i) = Wbool(a,i) + 1;
            W(Wbool(:,i)~=0,i) = G(:,Wbool(:,i)~=0) \ b(:,nn);
%             W(Wbool(:,i)~=0,i) = [G(:,Wbool(:,i)~=0);gamma*norm(G(:,Wbool(:,i)~=0))*eye(length(find(Wbool(:,i)~=0)))] \ [b(:,nn);zeros(length(find(Wbool(:,i)~=0)),1)];
            W(Wbool(:,i)~=0,i) = R(:,Wbool(:,i)~=0) \ w;
            restemp = norm(G*W(:,i)-b(:,nn))/norm(b(:,nn));
            if i>=2
                if ~isequal(find(Wbool(:,end)~=0),find(Wbool(:,end-1)~=0))
                    check = res(end-1)-res(end);
                else
                    check = tol+1;
                end
            end
            i = i+1;
        end
        istar = findchangepts(res,'Statistic','linear');
        if isempty(istar)
            istar = size(W,2);
        end
        wsparse_nn = W(:,istar);
        [wsparse_nn(wsparse_nn~=0),~] = sparsifyDynamics(G(:,wsparse_nn~=0),b(:,nn),lambda,gamma,M(wsparse_nn~=0,nn),inf);
        wsparse(:,nn) = wsparse_nn;
        Wcell{nn} = W;
    end
end

function W = subspacePursuitCV(G,b,M,lambda,maxits)

    supportList            = cell(lambda,1);
    crossValidationErrList = zeros(lambda,1);
    WColumnNorm            = vecnorm(G);
    
    for i=1:lambda
        support       = SPV2(G * diag(1./WColumnNorm), b ./norm(b,2),i);
        supportList{i}   = support;
        crossValidationErrList(i) = computeCrossValidationErrV2_dam(support, G, b, maxits);
    end
    [~, CrossIdx]=min(crossValidationErrList);
    supportPred = supportList{CrossIdx}';

    W = zeros(size(G,2),1);
    W(supportPred) = G(:,supportPred) \ b;
    if ~or(isempty(W),sum(M(:))==0)
        W = W.*M;
    end

end

function Ws=MIO(G,b,M,k_start,logtol,restol,eps1,eps2)

    J = size(G,2);
    n = size(b,2);
    bnd = max(abs(b'*G)./vecnorm(G).^2,2)';

    k = k_start;
    Wcell = cell(k,1);
    res = zeros(k,n);
    check = restol-1;
    while and(check < restol,k>=1)
        Ws = zeros(J,n);
        for nn=1:n
            ML = ones(J,1)*-bnd(nn);
            MU = ones(J,1)*bnd(nn);            
            cvx_solver Gurobi
            cvx_precision([eps1 eps1 eps2])
            cvx_begin
                variable z(J,1) binary
                variable W(J,1)
                minimize( norm(G*W - b(:,nn)) )
                subject to
                    {ML.*z <= W, W <= MU.*z, sum(z) <= k};
            cvx_end
            W(log10(abs(W))+logtol<max(log10(abs(W))))=0;
            res(k,nn) = norm(G*W - b(:,nn))/norm(b(:,nn));
            if ~or(isempty(W),sum(M(:))==0)
                W= W.*M;
            end
            Ws(:,nn) = W;
        end
        Wcell{k} = Ws;
        check = max(res(k,:));
        k = k-1;
    end
    Wcell = Wcell(k+1:end);
    res = res(k+1:end,:);

    Ws = zeros(J,n);
    for nn=1:n
        [~,ind] = max(diff(res(:,nn)));
        Ws(:,nn) = Wcell{ind}(:,nn);
    end

end

function support = SPV2(W,b,sparsity)
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
    
    x = (Phi_Lda'*Phi_Lda)\(Phi_Lda' * b);
    r = b - Phi_Lda*x;
    res = norm(r);
    iter = 0;
    
    if (res < 1e-12)
        X = zeros(N,1);
        X(Lda)=x;
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
        x_temp = (Phi_Sga'*Phi_Sga)\(Phi_Sga' * b);        
        [~, x_temp_index] = sort( abs(x_temp) , 'descend' );
        Lda = Sga(x_temp_index(1:sparsity));
        Phi_Lda = W(:,Lda);
        usedlda(Lda)=1;
        
        %%% calculate the residue
        x = (Phi_Lda'*Phi_Lda)\(Phi_Lda' * b);
        r = b - Phi_Lda*x;
        res = norm(r);
    
        X = zeros(N,1);
        X(Lda)=x;
        if ( res/res_old >= 1 || res < 1e-12)
            support = find(X);
            return
        end
    end
    support = find(X);

end

function [err] = computeCrossValidationErrV2_dam(ind, A, b, ntrials)
% script for computing cross validation error
errCrossAccu = zeros(1,ntrials);
for jj = 1:ntrials
    e=computeCVErrV4(ind,A, b);
    errCrossAccu(jj) =  e;
end
err = mean(errCrossAccu)  + std(errCrossAccu);
end


function [err]=computeCVErrV4(support,W,b,ratio)
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


% 
%             proj = R(inds(i),1:inds(i))./vecnorm(G(:,1:inds(i))).^2;
%             [~,a] = max(abs(proj));
%             Wbool(:,i) = Wbool(:,end);
%             Wbool(a,i) = Wbool(a,i) + 1;
%             W(Wbool(:,i)~=0,i) = R(:,Wbool(:,i)~=0) \ w;
