% driver script for the geneticDE algorithm
% Hans Pirnay
% 2011-03-03

%close all;
clear all;

n_scenarios = 1000;
n_stages = 4;
n_children = 3;

xi = randn(n_scenarios, n_stages);

p = 1/n_scenarios*ones(n_scenarios,1);
tic
tr = geneticDE(xi, p, n_children);
toc
% $$$ 
% $$$ tr = tree(n_stages+1, n_children, true);
% $$$ nps = n_children;
% $$$ idx = 1;
% $$$ for t=1:n_stages
% $$$     tr.node_values(idx+1:idx+nps) = xi(z(idx:idx+nps-1),t);
% $$$     idx = idx+nps;
% $$$     nps = nps*n_children;
% $$$ end

figure
tr.plot_tree
