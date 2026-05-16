function DCM = spm_dcm_seq_collapse(DCM)
% Collapse DCM.options.sequent.event(:).A, .B, .C into standard DCM.A, DCM.B, DCM.C
% C_mode: 'reflect' (default) keeps distinct inputs, aligned by names if present
%         'any' collapses all inputs into a single column
%
% Output
%   DCM.A = {A1, A2, A3[, A4]}
%   DCM.B = {B1, B2, ..., Bq}               % supports multiple B matrices
%   DCM.C = C                                % n by nInputs (or n by 1 if 'any')
% 
% if nargin < 2 || isempty(C_mode),
 C_mode = 'reflect';

% guards
if ~isfield(DCM,'options') || ~isfield(DCM.options,'sequent') ...
        || ~isfield(DCM.options.sequent,'event')
    error('No DCM.options.sequent.event found.');
end
ev = DCM.options.sequent.event;
if ~isstruct(ev) || isempty(ev)
    error('DCM.options.sequent.event must be a nonempty struct array.');
end

% sizes from first event
A0 = ev(1).A;
if ~iscell(A0) || isempty(A0)
    error('event(1).A must be a nonempty cell array.');
end
na = numel(A0);
n  = size(A0{1},1);
for k = 1:na
    if ~isequal(size(A0{k}), [n n])
        error('A{%d} in event(1) is not %d-by-%d.', k, n, n);
    end
end

% number of B matrices across events
nbmax = 0;
for e = 1:numel(ev)
    if isfield(ev(e),'B') && ~isempty(ev(e).B)
        if iscell(ev(e).B), nbmax = max(nbmax, numel(ev(e).B));
        else,               nbmax = max(nbmax, 1);
        end
    end
end
if nbmax == 0, nbmax = 1; end

% collect input names if present and consistent with C width
useNames  = false;
allNames  = {};
for e = 1:numel(ev)
    if isfield(ev(e),'u_names') && ~isempty(ev(e).u_names)
        if isstring(ev(e).u_names), names_e = cellstr(ev(e).u_names(:));
        else,                       names_e = ev(e).u_names(:);
        end
        if numel(names_e) == size(ev(e).C,2)
            allNames = [allNames; names_e]; %#ok<AGROW>
            useNames = true;
        else
            useNames = false;
            break
        end
    end
end
if useNames
    [globalNames, ~, ~] = unique(allNames, 'stable');
    nInputs = numel(globalNames);
else
    nInputs = 0;
    for e = 1:numel(ev)
        nInputs = max(nInputs, size(ev(e).C,2));
    end
    globalNames = {};
end

% accumulators
A = cell(1, na);
for k = 1:na
    A{k} = false(n, n);
end
B = repmat({false(n, n)}, 1, nbmax);
if strcmpi(C_mode, 'any')
    C = false(n, 1);
else
    C = false(n, nInputs);
end

% accumulate with OR
for e = 1:numel(ev)
    % A
    for k = 1:na
        Ak = logical(ev(e).A{k});
        if ~isequal(size(Ak), [n n])
            error('A{%d} in event(%d) is not %d-by-%d.', k, e, n, n);
        end
        A{k} = A{k} | Ak;
    end

    % B
    if isfield(ev(e),'B') && ~isempty(ev(e).B)

        if iscell(ev(e).B)

            this_nb = numel(ev(e).B);

            if this_nb == 1 && nbmax > 1

                Bb = logical(ev(e).B{1});

                if isfield(ev(e),'Bint') && ~isempty(ev(e).Bint)
                    Bb = Bb | diag(~cellfun(@isempty,ev(e).Bint{1}));
                end



                if ~isequal(size(Bb), [n n])
                    error('B{1} in event(%d) is not %d-by-%d.', e, n, n);
                end
                for b = 1:nbmax
                    B{b} = B{b} | Bb;
                end

            elseif this_nb == nbmax

                for b = 1:nbmax
                    Bb = logical(ev(e).B{b});

                    if isfield(ev(e),'Bint') && ~isempty(ev(e).Bint)
                        Bb = Bb | diag(~cellfun(@isempty,ev(e).Bint{1}));
                    end

                    if ~isequal(size(Bb), [n n])
                        error('B{%d} in event(%d) is not %d-by-%d.', b, e, n, n);
                    end
                    B{b} = B{b} | Bb;
                end

            else
                error('Number of B matrices in event(%d) is %d but expected %d.', e, this_nb, nbmax);
            end

        else
            Bb = logical(ev(e).B);

            if isfield(ev(e),'Bint') && ~isempty(ev(e).Bint)
                Bb = Bb | diag(~cellfun(@isempty,ev(e).Bint{1}));
            end


            if ~isequal(size(Bb), [n n])
                error('B in event(%d) is not %d-by-%d.', e, n, n);
            end
            for b = 1:nbmax
                B{b} = B{b} | Bb;
            end
        end
    end

    % C
    Ce = logical(ev(e).C);
    if size(Ce,1) ~= n
        error('C in event(%d) has %d rows. Expected %d.', e, size(Ce,1), n);
    end
    if strcmpi(C_mode, 'any')
        C = C | any(Ce, 2);
    else
        if useNames
            if isstring(ev(e).u_names), names_e = cellstr(ev(e).u_names(:));
            else,                       names_e = ev(e).u_names(:);
            end
            if size(Ce,2) < numel(names_e)
                Ce(:, end+1:numel(names_e)) = false;
            elseif size(Ce,2) > numel(names_e)
                Ce = Ce(:, 1:numel(names_e));
            end
            for j = 1:numel(names_e)
                idx = find(strcmp(globalNames, names_e{j}), 1);
                if ~isempty(idx)
                    C(:, idx) = C(:, idx) | Ce(:, j);
                end
            end
        else
            if size(Ce,2) < nInputs
                Ce(:, end+1:nInputs) = false;
            end
            C = C | Ce(:, 1:nInputs);
        end
    end
end

% cast to double 0 and 1
for k = 1:na
    A{k} = double(A{k});
end
for b = 1:nbmax
    B{b} = double(B{b});
end
C = double(C);

% assign
DCM.A = A;
DCM.B = B;
DCM.C = C;

end
