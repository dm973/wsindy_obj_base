% w = reshape_cell(WS.weights,cellfun(@(L)length(L.terms),WS.lib))

function y = reshape_cell(x,inds)
    x = x(:);
    y = cell(length(inds),1);
    for i=1:length(inds)
        y{i} = x(1:inds(i));
        x = x(inds(i)+1:end);
    end
end