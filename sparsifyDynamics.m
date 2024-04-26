function [Xi,its] = sparsifyDynamics(Theta,dXdt,lambda,n,gamma,M,maxits,toggle_jointthresh,reg_inds)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% Modified by Daniel A. Messenger, 2020 

[~,nn] =size(Theta);
if  gamma ~= 0
    Theta = [Theta;gamma*eye(nn)];
    dXdt = [dXdt;zeros(nn,n)];
end

if ~isempty(reg_inds) % initial guess: truncated Least-squares
    if isequal(class(reg_inds),'double')
        reg0 = rank(Theta,norm(Theta)*reg_inds);
        reg_inds = abs(dXdt'*Theta)./vecnorm(Theta)/norm(dXdt);
        [~,reg_inds] = sort(reg_inds,'descend');
        reg_inds = reg_inds(1:min(reg0,end));
    elseif isequal(reg_inds,'coltrim')
        reg_inds = coltrim(Theta,lambda,dXdt);
    end
    Xi = zeros(size(Theta,2),1);
    Xi(reg_inds) = Theta(:,reg_inds) \ dXdt;
else % initial guess: Least-squares
    Xi = Theta \ dXdt;
end

if ~isempty(M)
    Xi = M.*Xi;
    % bnds = norm(dXdt)^2./abs(dXdt'*Theta)'.*M;
    if toggle_jointthresh == 1
        bnds = norm(dXdt)./vecnorm(Theta)'.*M;
        LBs = lambda*max(1,bnds);
        UBs = 1/lambda*min(1,bnds);
    elseif toggle_jointthresh == 2
        bnds = norm(dXdt)./vecnorm(Theta)'.*M;
        LBs = lambda*bnds;
        UBs = 1/lambda*bnds;
    elseif toggle_jointthresh == 3
        bnds=norm(dXdt)^2./abs(dXdt'*Theta)'.*M;
        UBs = 1/lambda*bnds;
        LBs = lambda*bnds;
    else
        bnds = norm(dXdt)./vecnorm(Theta)'.*M;
        w0 = abs(dXdt'*Theta);
        nrms = vecnorm(Theta);
        alpha = max(w0./nrms.^2.*M');
        beta = max(w0./nrms/norm(dXdt));
        LBs = lambda*max(alpha,bnds*beta);
        UBs = 1/lambda*min(alpha,bnds*beta);
    end
    thrs_EL = [LBs bnds UBs];
else
    thrs_EL = [];
end

smallinds = 0*Xi;
for j=1:min(nn,maxits)
    if ~isempty(M)
        smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
        if all(smallinds_new(:)==smallinds(:))
            its = j;
            return
        else
            smallinds = smallinds_new;
            Xi(smallinds)=0;    
            for ind=1:n
                Xi(~smallinds,ind) = M(~smallinds).*(Theta(:,~smallinds)\dXdt(:,ind));
            end
        end
    else
        smallinds_new = (abs(Xi)<lambda);
        if all(smallinds_new(:)==smallinds(:))
            its = j;
            return
        else
            smallinds = smallinds_new;
            Xi(smallinds)=0;
            for ind = 1:n        
                biginds = ~smallinds(:,ind);
                Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
            end
        end
    end
end
its = j;
end
