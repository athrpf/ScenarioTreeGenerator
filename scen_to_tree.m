function tr = scen_to_tree(xi, p, nchildren, combo)
% Computes a tree from the scenarios and combinations. If no
% combinations are given, random ones are used.
% USAGE: 
% function tr = scen_to_tree(xi, nchildren, combo)
%
% xi : n_scen times nstages

[nscen, nstages] = size(xi);
assert(length(p)==nscen);

ntreescen = nstages^nchildren;
assert(nscen>=ntreescen);
randomize = false;
if nargin<4
    randomize = true;
    combo = zeros(nscen, 1);
    combo(1:ntreescen) = 1:ntreescen;
    combo(ntreescen+1:end) = randi(nscen, nscen-ntreescen,1);
end
assert(length(combo)==nscen);

tr = tree(nstages+1, nchildren, true);

for stage=1:nstages
    nts  = tr.nodes_this_stage(stage+1); % number of nodes this stage
    fnts = tr.first_node_this_stage(stage+1); % first node this stage
    spn  = ntreescen/nts; % scenarios per node
    for node=1:nts
        ind = find((combo>spn*(node-1))&(combo<=spn*node));
        % node_val = sum(p(ind).*xi(ind,stage))/sum(p(ind))
        node_val = 0;
        if randomize
            node_val = xi(ind(randi(length(ind))),stage);
        else
            c = squareform(pdist(xi(ind,stage)));
            [minc, minidx] = min(sum(c)');
            minidx;
            ind(minidx);
            node_val = xi(ind(minidx),stage);
        end
        tr.node_values(fnts+node-1) = node_val;
    end
end
