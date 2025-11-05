classdef wendyA_model < wsindy_model
    properties
        covfun
    end

    methods
        function obj = wendyA_model(dat,lib,tf,covfun,varargin)
            obj = obj@wsindy_model(dat,lib,tf,varargin{:});
            obj.covfun = covfun;
        end
    end

    methods
        function obj = get_cov(obj,w)
            if exist('w','var')
                obj.weights = w;
            end
            if any(obj.weights)
                obj.cov = obj.covfun(obj);
            else
                obj.cov = speye(sum(cellfun(@(b)length(b),obj.bs{1})));
            end
        end

        function [G,b,RT] = apply_cov(obj,G,b,diag_reg,sparse_inds)
            if isempty(obj.cov)
                obj = obj.get_cov;
            end
            check = 0;
            while check == 0
                try
                    RT = (1-diag_reg)*obj.cov + (diag_reg*mean(diag(obj.cov)))*speye(size(obj.cov,1));
                    RT = chol( RT )';
                    check = 1;
                catch
                    diag_reg = diag_reg*2;
                    fprintf('\nincreasing Cov regularization to %0.2e\n',diag_reg)
                end
            end
            G = RT \ G;
            b = RT \ b;
        end

    end
end