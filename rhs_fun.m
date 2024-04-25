function dx = rhs_fun(features,params,x)
    x = num2cell(x);
    dx = zeros(length(params),1);
    for i=1:length(params)
        if ~isempty(features{i})
            dx(i) = cellfun(@(z1) z1(x{:}),features{i})*params{i}(:);
        end
    end
end