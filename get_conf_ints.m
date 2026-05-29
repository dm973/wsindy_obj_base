function [conf_int,c,confbounds,xflip,coverage] = get_conf_ints(w_hat,w_true_hat,CovW,c,c_maxmin)
    xflip = [1:length(w_hat) length(w_hat):-1:1];
    % <(100)c % chance of not containing true val
    check = true;
    stdW = max(sqrt(diag(CovW)),eps);

    if c_maxmin
        while check
            c = c/2;
            conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
            confbounds = [w_hat-conf_int;flipud(w_hat+conf_int)];        
            if any(~and(w_true_hat>=w_hat-conf_int,w_true_hat<=w_hat+conf_int))
                check = true;
            else
                check = false;
            end
        end
    else
        conf_int = arrayfun(@(x)norminv(1 - c/2,0,x),stdW);
        confbounds = [w_hat-conf_int;flipud(w_hat+conf_int)];
        coverage = min(w_true_hat - (w_hat-conf_int), w_hat+conf_int - w_true_hat)./(2*conf_int);
    end

end