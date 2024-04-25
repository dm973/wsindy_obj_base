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


function [Xi,its,thrs_EL] = sparsifyDynamics_qp(Theta,dXdt,lambda,gamma,M,A,b,excl_inds,const_tol,max_its,disp_opt,max_its_stls)
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
        const_tol=10^-10;
    end
    options = optimoptions('quadprog','Display',disp_opt,'ConstraintTolerance',const_tol,'MaxIterations',max_its);
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
