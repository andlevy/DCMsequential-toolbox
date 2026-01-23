function cX = spm_cell(x,X)



% Get the sizes of each cell in M.x
cellSizes = cellfun(@length, X);

% Verify that the total columns match
totalCols = sum(cellSizes);
if size(x, 2) ~= totalCols
    error('Matrix x has %d columns but M.x expects %d columns', ...
        size(x, 2), totalCols);
end

% Split x into cells using mat2cell
cX = mat2cell(x, size(x, 1), cellSizes);