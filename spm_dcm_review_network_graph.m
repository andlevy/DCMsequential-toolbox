function fig = spm_dcm_review_network_graph(DCM,seq,post)

if isfield(DCM,'Ipos') && ~isempty(DCM.Ipos)
    Ipos  = DCM.Ipos;
else
    Ipos  = DCM.Lpos;
end

Ipos = Ipos;

if nargin < 2
    seq = false;
end

if nargin < 3
    post = false;
end

Sname = DCM.Sname;
fig   = {};
if seq

    event   = DCM.options.sequent.event;
    nmm     = {DCM.options.model.source};
    nEvents = numel(event);

    for e = 1:nEvents

        A = event(e).A;
        B = event(e).B;
        C = event(e).C;

        if isfield(event,'Bint') && ~isempty(event(e).Bint)
            B{1} = B{1} | diag(~cellfun(@isempty,event(e).Bint{1}));
        end

        
        if post

            Q.P = DCM.Pp.B{e};
            Q.B = DCM.Ep.B{e};

            if isfield(DCM.Ep,'int') && isfield(DCM.Ep.int{1},'B') && ~isempty(DCM.Ep.int{1}.B)

                b = []; p =[];
                for i = 1:numel(DCM.Ep.int)
                    b(i) = DCM.Ep.int{i}.B{e}(1);
                    p(i) = DCM.Pp.int{i}.B{e}(1);
                end

                Q.B{1}(logical(eye(size(Q.B{1})))) = b;
                Q.P{1}(logical(eye(size(Q.P{1})))) = p;

            end

        else
            Q   = [];
        end

        if isfield(event,'fu') && ~isempty(event(e).fu)
            fu = event(e).fu;
        else
            fu = [];
        end

        fig{end+1} = draw_network(Ipos,Sname,A,B,C,Q,nmm,fu);

        title(event(e).name)

        if post
            fig{end}.Name = ['network ' sprintf('%02d ',e) lower(event(e).name) 'posterior'];
        else
            fig{end}.Name = ['network ' sprintf('%02d ',e) lower(event(e).name) ];
        end


    end

else

    A     = DCM.A;
    B     = DCM.B;
    C     = DCM.C;

    if post
        Q.A = DCM.Ep.A;
        Q.B = DCM.Ep.B;
    else
        Q = [];
    end

    fig = draw_network(Ipos,Sname,A,B,C,Q);

end


end
%% ========================================================================
%                    Main plotting routines
%==========================================================================
function fig = draw_network(Ipos,Sname,A,B,C,Q,nmm,fu)

if nargin < 6
    Q = [];
end
if nargin < 7
    nmm = [];
end
if nargin < 8
    fu = [];
end

% figure
if exist('figtools_figure','file')
    fig = figtools_figure();
else
    fig = figure;
end
set(fig, 'WindowState', 'maximized');
hold on

N = size(Ipos,2);

% style
Ssep    = 18;                  % separation between forward and backward at midpoint
ah      = 5;                   % arrowhead length (data units)
theta   = deg2rad(18);         % arrowhead half-angle
lw_base = 7;                   % baseline width
lw_mod  = 7;                   % modulation overlay width
r       = 5;

col_mod_fwd  = [0 .7 .4];
col_mod_bck  = [.8 .2 .6];
col_mod_slf  = [.2 .2 .9];
col_base     = [.5 .5 .5];
st_base      = '-';
st_mod       = '--';

% nodes styling
%--------------------------------------------------------------------------
if ~isempty(nmm) && numel(unique(nmm)) > 1
    col_nmm = [];
    for n = 1:numel(nmm)
        switch lower(nmm{n})
            case 'cmc'
                col_nmm(end+1,:) = [0.95 0.80 0.95]; % Lavender/Pink
            case 'bgt'
                col_nmm(end+1,:) = [1.00 0.95 0.80]; % Pale Gold
            case 'mmc'
                col_nmm(end+1,:) = [0.80 1.00 0.80]; % Pale Green
            case 'bcn'
                col_nmm(end+1,:) = [0.80 0.85 1.00]; % Periwinkle Blue
        end
    end

else
    col_nmm = repmat([0.85 0.85 0.85],N,1);
end

ar = []; tx = [];

% ---------- draw baseline extrinsic edges from A{1} and A{2} ----------
if numel(A) >= 1
    Afwd = A{1} ~= 0;
else
    Afwd = false(N);
end
if numel(A) >= 2
    Abwd = A{2} ~= 0;
else
    Abwd = false(N);
end

% forward i->j
for i = 1:N
    for j = 1:N
        if i==j, continue; end
        if Afwd(i,j)
            ar(end+1) = draw_edge_bezier(Ipos, j, i, r, Ssep, ah, theta, col_base, lw_base,st_base);
        end
    end
end
% backward i->j
for i = 1:N
    for j = 1:N
        if i==j, continue; end
        if Abwd(i,j)
            ar(end+1) = draw_edge_bezier(Ipos, j, i, r, Ssep, ah, theta, col_base, lw_base,st_base);
        end
    end
end

% ---------- overlay B modulations on existing baseline edges ----------
for k = 1:numel(B)
    Bk = B{k} ~= 0;
    for i = 1:N
        for j = 1:N
            if i==j, continue; end
            % only overlay where the baseline edge exists in the same direction
            if Afwd(i,j) && Bk(i,j)
                ar(end+1) = draw_edge_bezier(Ipos, j, i, r, Ssep, ah, theta, col_mod_fwd(k,:), lw_mod,st_mod);
                if ~isempty(Q)
                    tx(end+1) = add_param_label_edge(Ipos,j,i,r,Ssep,1,Q,k,col_mod_fwd(k,:));
                end

            elseif Abwd(i,j) && Bk(i,j)
                ar(end+1) = draw_edge_bezier(Ipos, j, i, r, Ssep, ah, theta, col_mod_bck(k,:), lw_mod,st_mod);
                if ~isempty(Q)
                    tx(end+1) = add_param_label_edge(Ipos,j,i,r,Ssep,2,Q,k,col_mod_bck(k,:));
                end

            end
        end
    end
end

% ---------- self-connections ----------
Abase_self = diag(ones(N,1));
for i = 1:N
    if Abase_self(i,i)
        ar(end+1) = draw_selfloop(Ipos, i, r, ah, theta, col_base, lw_base,st_base);
    end
end
% overlay self-loop modulations if any B{k}(i,i) is non-zero
for k = 1:numel(B)
    idx = find(diag(B{k}) ~= 0);
    for ii = reshape(idx,1,[])
        ar(end+1) = draw_selfloop(Ipos, ii, r, ah, theta, col_mod_slf(k,:), lw_mod,st_mod);
        if ~isempty(Q)
            tx(end+1) = add_param_label_self(Ipos, ii, r, Q, k, col_mod_slf(k,:));
        end

    end
end

% ---------- inputs from C ----------
% nodes that receive any input across columns of C
Cin = false(N,1);
if ~isempty(C), Cin = any(C ~= 0, 2); end

% orientation references
cx   = mean(Ipos(1,:));        % network centre x
ymin = min(Ipos(2,:));         % lowest node y
Lin  = 12;                     % input arrow length in data units
col_in = [255, 128, 0]./255;
lw_in  = lw_base;

for i = 1:N
    if Cin(i)
        ar(end+1) = draw_input(Ipos,i,r,Lin,ah,theta,col_in,lw_in,cx,ymin,fu);
    end
end

% ---------- nodes and labels ----------
for i = 1:N
    rectangle('Position',[Ipos(1,i)-r, Ipos(2,i)-r, 2*r, 2*r], ...
        'Curvature',[1 1], 'FaceColor',col_nmm(i,:), 'EdgeColor','k','LineWidth',3);
end
for i = 1:numel(Sname)
    
    if contains(Sname{i},'(silent)','ignorecase',true)
        nodename = { extractBefore(Sname{i},'('), ['(' extractAfter(Sname{i},'(')] };
    else
        nodename = Sname{i};
    end

    text(Ipos(1,i), Ipos(2,i), nodename, ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'fontweight','bold','FontSize',11);

end

for a = ar
    uistack(a,'top')
end

for t = tx
    uistack(t,'top')
end


% ---------- legend ----------
show_nmm = ~isempty(nmm) && numel(unique(lower(nmm))) > 1;
add_network_legend(col_mod_fwd, col_mod_bck, col_mod_slf, st_mod, lw_mod, show_nmm, nmm);


% ---------- axes: equal data units and square view ----------
axis equal
pad  = r + Ssep/2 + ah;
xmin = min(Ipos(1,:)) - pad; xmax = max(Ipos(1,:)) + pad;
ymin = min(Ipos(2,:)) - pad; ymax = max(Ipos(2,:)) + pad;
cx = (xmin+xmax)/2; cy = (ymin+ymax)/2;
half = max(xmax-xmin, ymax-ymin)/2;
xlim([cx-half, cx+half]); ylim([cy-half, cy+half]);
pbaspect([1 1 1])
axis off
end

function tf = diag_any(Acell,idx)
% true for i if any A{idx}(i,i) is non-zero, returns 1xN logical
if numel(Acell) < idx || isempty(Acell{idx}), tf = false(1,0); return; end
a = Acell{idx};
tf = diag(a) ~= 0;
tf = tf(:).';
end

function ar =  draw_edge_bezier(Ipos, i, j, r, Ssep, ah, theta, col, lw,st)
% Draw a directed quadratic Bezier from node i to node j
p = Ipos(1:2,i)';  q = Ipos(1:2,j)';
v = q - p;  L = norm(v);
if L <= eps, return; end
u = v / L;

% trim so curve meets node borders
p2 = p + r*u;
q2 = q - r*u;

% constant midpoint offset gives fixed separation between reciprocals
n  = [-u(2), u(1)];        % unit normal tied to direction
m  = 0.5*(p2 + q2);
h  = Ssep/2;
c  = m + h*n;

% sample Bezier
tt = linspace(0,1,80)';
B  = (1-tt).^2.*p2 + 2*(1-tt).*tt.*c + tt.^2.*q2;
plot(B(:,1), B(:,2), 'Color', col, 'LineWidth', lw,'LineStyle',st)

% arrowhead at end
d  = q2 - c;
ar = draw_arrowhead(q2, d, ah, theta, col);
end

function ar = draw_arrowhead(q, d, ah, theta, col)
% Triangular arrowhead at tip q along direction d
if ~all(isfinite(d)) || norm(d) <= eps, return; end
u = d(:)/norm(d);
w = [-u(2); u(1)];
tip  = q(:);
base = tip - ah*u;
s    = ah*tan(theta);
p1   = base + s*w;
p2   = base - s*w;
ar = patch([tip(1) p1(1) p2(1)], [tip(2) p1(2) p2(2)], col, ...
    'EdgeColor', col, 'FaceColor', col);
end

function ar = draw_selfloop(Ipos, i, r, ah, theta, col, lw,st)
% Arc above the node with arrow at end
rad = 0.8*r;
ctr = [Ipos(1,i) - 1.2*r, Ipos(2,i) + 1.5*rad];
t   = linspace(pi*-0.2, pi*1.5, 80);
x   = ctr(1) + rad*cos(t);
y   = ctr(2) + rad*sin(t);
plot(x, y, 'Color', col, 'LineWidth', lw,'LineStyle',st)
d   = [x(end)-x(end-1); y(end)-y(end-1)];
ar = draw_arrowhead([x(end)+rad*.9; y(end)], d, ah, theta, col);
end

function ar = draw_input(Ipos,i,r,Lin,ah,theta,col,lw,cx,ymin,fu)

% Straight input arrow to node i.
x   = Ipos(1,i);
y   = Ipos(2,i);
tol = 1e-6;
r = r+1;
if y <= ymin + tol
    % lowest node(s) → vertical arrow from below
    p0 = [x, y - r - Lin];
    p1 = [x, y - r];
else
    if x <= cx
        % node on left → horizontal arrow from left
        p0 = [x - r - Lin, y];
        p1 = [x - r,       y];
    else
        % node on right → horizontal arrow from right
        p0 = [x + r + Lin, y];
        p1 = [x + r,       y];
    end
end

dir    = p1 - p0;
u      = dir / norm(dir);
offset = ah/-2;                   % adjust until it looks right
p1_adj = p1 - u * offset;

plot([p0(1) p1(1)], [p0(2) p1(2)], 'Color', col, 'LineWidth', lw)
ar = draw_arrowhead(p1_adj, dir, ah, theta, col);


% Optional: add label for fu
if nargin > 9 && ~isempty(fu)
    % Midpoint of the arrow, offset slightly perpendicular
    mid   = (p0 + p1)/2;
    perp  = [-u(2), u(1)];         % perpendicular to arrow
    offset = 0.6 * r;              % tweak spacing
    pos   = mid + perp*offset;
    pos   = pos - [2 0];

    switch fu
        case 'spm_erp_u'
            label = 'ERP';
        case 'spm_scp_u'
            label = 'SCP';
        otherwise
            label = '';
    end

    if ~isempty(label)
        text(pos(1), pos(2), label, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontSize',11, ...
            'Color',col, ...
            'FontWeight','bold');
    end
end
end

function tx = add_param_label_edge(Ipos,i,j,r,Ssep,dir,Q,k,col)
% Text label "B (A)" near the midpoint of edge i->j.
% dir = 1 for forward (A{1}), 2 for backward (A{2})

if isempty(Q) || ~isfield(Q,'B') || numel(Q.B) < k, return; end
Bij = Q.B{k}(j,i);
Pij = NaN;
if isfield(Q,'P')
    Pij = Q.P{k}(j,i);
end

% compose text: show B; include (A) only if available
if ~isnan(Pij)
    txt = {sprintf('%.2f', Bij) sprintf('(%.2f)', Pij)};
else
    txt = sprintf('%.2f', Bij);
end

% reproduce the geometry from draw_edge_bezier to get the midpoint
p = Ipos(1:2,i)';  q = Ipos(1:2,j)';
v = q - p;  L = norm(v);
if L <= eps, return; end
u = v / L;
p2 = p + r*u;  q2 = q - r*u;
n  = [-u(2), u(1)];
m  = 0.5*(p2 + q2);
h  = Ssep/2;
c  = m + h*n;
tt = 0.35;                                           % midpoint of Bezier
Bmid = (1-tt)^2*p2 + 2*(1-tt)*tt*c + tt^2*q2;

% offset slightly outward along the normal
pos = Bmid;% + 1*r*n;

tx = text(pos(1), pos(2), txt, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
     'FontWeight','bold',...
    'FontSize',10, ...
    'Color', col, ...
    'BackgroundColor',[1 1 1 .8],...
    'Interpreter','none');
end

function tx = add_param_label_self(Ipos,i,r,Q,k,col)
% Text label for self-connection i (uses B; includes A if available)

if isempty(Q) || ~isfield(Q,'B') || numel(Q.B) < k, return; end
Bii = Q.B{k}(i,i);


Pii = NaN;
if isfield(Q,'P')
    Pii = Q.P{k}(i,i);
end

if ~isnan(Pii)
    txt = {sprintf('%.2f', Bii) sprintf('(%.2f)',Pii)};
else
    txt = sprintf('%.2f', Bii);
end

% pick a point just above and to the right of the self-loop
rad = .8*r;
ctr = [Ipos(1,i) - 2.2*r, Ipos(2,i) + 1.5*rad];
t   = pi*1.1;
pos = [ctr(1) + rad*cos(t), ctr(2) + rad*sin(t)] + [0.55*r, 0.35*r];

tx = text(pos(1), pos(2), txt, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontWeight','bold',...
    'FontSize',10, ...
    'Color', col, ...
    'BackgroundColor',[1 1 1 .8],...
    'Interpreter','none');
end

function add_network_legend(colF, colB, colS, st_mod, lw_mod, show_nmm, nmm)
% Legend entries for modulated edges and (optionally) NMM colour mapping.

% Prototype handles for edge modulations
hF = line(nan, nan, 'LineStyle', st_mod, 'LineWidth', lw_mod, 'Color', colF, 'DisplayName','B forward');
hB = line(nan, nan, 'LineStyle', st_mod, 'LineWidth', lw_mod, 'Color', colB, 'DisplayName','B backward');
hS = line(nan, nan, 'LineStyle', st_mod, 'LineWidth', lw_mod, 'Color', colS, 'DisplayName','B self');

H   = [hF, hB, hS];
Lbl = {'Forward','Backward','Self'};

% Optional NMM swatches when multiple types are present
if show_nmm
    % normalise names and get unique list (preserve first-seen order)
    nmml = lower(nmm(:));
    [uniqNames, ia] = unique(nmml, 'stable');

    getCol = @(name) switchlower(name, ...
        'cmc', [0.95 0.80 0.95], ...  % light blue
        'erp', [0.85 1 .85],     ...  % light green
        'bgt', [1.00 0.95 0.80], ...  % light orange
        'mmc',[0.80 1.00 0.80] , ...  % light violet
        'bcn',[0.80 0.85 1.00] , ...
        [0.85 0.85 0.85]);             % default grey

    for t = 1:numel(uniqNames)
        col = getCol(uniqNames{t});
        % round node swatch with black edge
        h = line(nan, nan, ...
            'LineStyle','none', 'Marker','o', 'MarkerSize',10, ...
            'MarkerFaceColor', col, 'MarkerEdgeColor','k', ...
            'DisplayName', upper(uniqNames{t}));
        H(end+1)   = h; %#ok<AGROW>
        Lbl{end+1} = upper(uniqNames{t}); %#ok<AGROW>
    end
end

lgd = legend(H, Lbl, 'Location','northeastoutside', 'AutoUpdate','off');
set(lgd, 'Box','off');
end

function out = switchlower(key, varargin)
% switch on lowercased string keys: switchlower(key,'a',valA,'b',valB,defaultVal)
k = lower(key);
pairs = varargin(1:end-1);
out   = varargin{end};
for i = 1:2:numel(pairs)
    if strcmp(k, pairs{i})
        out = pairs{i+1};
        return
    end
end
end
