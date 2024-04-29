
classdef sourceterm < term
    methods
        function obj = sourceterm(varargin)
            obj@term(varargin{:}); 
        end
    end

    methods
        function Y = evalterm(obj,dat) % only fHandle
            grdz = dat.grid;
            grdz{1} = grdz{1}(:);
            for j=2:dat.ndims
                ones_cell = repmat({1},1,j-1);
                grdz{j} = reshape(grdz{j},ones_cell{:},length(grdz{j}));
            end
            Y = obj.fHandle(dat.Uobs{:},grdz{:});
        end
    end

end
