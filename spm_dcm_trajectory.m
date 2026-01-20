function Q = spm_dcm_trajectory(P,v,logscale)

if nargin < 3
    logscale = true;
end

if logscale
    Q = exp(P{1});
else
    Q = P{1};
end

for k = 1:size(v,2)

    if logscale
        q1 = exp(P{k+1});
        q2 = exp(P{k});
    else
        q1 = P{k+1};
        q2 = P{k};
    end

    Q  = Q + (q1-q2) *v(:,k);

end