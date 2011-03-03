% unit test for tree class
clear all;
close all;

n_children = 3;
n_stages = 4;

tr = tree(n_stages, n_children, true);

% testing functions specific to nodes
for node=1:tr.n_nodes
    stage = tr.stage_of_node(node);
    children = tr.children(node);
    parent   = tr.parent(node);
    disp(['node ',num2str(node),', stage ', num2str(stage), ', parent ',num2str(parent), ', children: ',num2str(children)])
end


% testing functions specific to stages
for stage=1:n_stages
    nts = tr.nodes_this_stage(stage);
    fnts = tr.first_node_this_stage(stage);
    disp(['stage ', num2str(stage), ', nts ',num2str(nts),', fnts ' ...
                        , num2str(fnts)])
end

