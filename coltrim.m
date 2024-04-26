function inds = coltrim(A,lam,y,exinds)
    if ~exist('exinds','var')
        exinds = [];
    end
    B = abs(((A'*A-diag(diag(A'*A)))./vecnorm(A))./(vecnorm(A)'));
    proj = abs(y'*A./vecnorm(A));
    inds = 1:length(B);
    C = B > lam;
    check = 1;
    while all([any(C(:)) length(inds)>1 check])
        ind2 = inds(sum(C,2)~=0);            
        ind2 = ind2(~ismember(ind2,exinds));
        [~,I] = min(proj(ind2));
        if ~isempty(I)
            a = find(ind2(I)==inds);
            inds = inds([1:a-1 a+1:end]);
            B = B([1:a-1 a+1:end],[1:a-1 a+1:end]);
            C = C([1:a-1 a+1:end],[1:a-1 a+1:end]);
        else
            check = 0;
        end
        % C = B>lam;
        % imagesc(C)
        % drawnow
        % inds
    end
end