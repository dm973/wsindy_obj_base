classdef testfcn_Vcell < testfcn

    properties
        alphas
        Vcell
    end

    methods
        
        function obj = testfcn_Vcell(dat,Vcell,alphas)
            obj = obj@testfcn(dat);
            
            obj.Vcell = Vcell;
            obj.alphas = alphas;
        end

        function V = get_testmat(obj,difftag)
            if isequal(class(difftag),'diffOp')
                difftag = difftag.difftags;
            elseif length(difftag)~=obj.ndims
                difftag = zeros(obj.ndims,1);
            end
            V = 1;
            for i=1:obj.ndims
                V = kron(obj.Vcell{i}{obj.alphas{i} == difftag(i)},V);
            end
        end

        function vec = test(obj,dat,term)
            if or(isequal(class(term),'diffOp'),any(cellfun(@(x) isequal(x,'diffOp'),superclasses(term)))) % diffop given
                V = obj.get_testmat(term.difftags);
                vec = V*reshape(dat.Uobs{term.stateind},[],1);
            elseif isequal(class(term.linOp),'diffOp') % term has a linear operator
                V = obj.get_testmat(term.linOp.difftags);
                vec = V*reshape(term.evalterm(dat),[],1);
            elseif isempty(term.linOp) % term is a linear operator
                V = obj.get_testmat(zeros(obj.ndims,1));
                vec = V*reshape(term.evalterm(dat),[],1);
            end
        end

    end
end