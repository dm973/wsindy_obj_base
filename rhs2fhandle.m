function f = rhs2fhandle(f)
    f = @(varargin) squeeze(reshape(arrayfunvec(cell2mat(arrayfun(@(i)varargin{i}(:),1:length(varargin),'uni',0)),...
        @(z)f(z)',2),[size(varargin{1}) length(varargin)]));
end