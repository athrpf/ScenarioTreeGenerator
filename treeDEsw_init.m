close all;
clear all;

% This script acts as a driver for the function treeDEsw
% 

addpath('/home/hans/MATLAB/externals/gurobi_mex_v1.50');

n_scen = 1000;
n_stages = 3;
n_children = 3;

xi_base = randn(n_stages, n_scen);

xi = ones(n_stages+1, n_scen);

for stage=2:n_stages+1
    xi(stage,:) = xi(stage-1,:).*xi_base(stage-1,:);
end

xi = xi(2:end,:);

tic
[tr,nodes] = treeDEsw(xi, n_children);
toc

figure;
hold on;
colors = 'ymcrgbwk';
for kk=2:tr.n_nodes
    stage = tr.stage_of_node(kk);
    pointtype = [colors(kk-1),'x'];
    plot(stage*ones(1,length(nodes{kk}))'+0.05*(rand-.5),xi(stage-1, nodes{kk})',pointtype);
end

tr.plot_tree
hold off
