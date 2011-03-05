function diff = compareDEtrees(xi, p, n_children)
% Hans Pirnay
% 2011-03-03
% 
% Compares the tree generation algorithms for discrete-event trees
% Usage:
% diff = compareDEtrees(xi, p, n_children)
% xi : scenarios representing the stochastic process. Format
%      scenarios x stages

tic
tr_genetic = geneticDE(xi, p, n_children);
toc

tic
tr_milp    = treeDEsw(xi, p,n_children, 1);
toc

kantoro_genetic = tr_genetic.kantorovich(xi, p, 1);
kantoro_milp = tr_milp.kantorovich(xi, p, 1);

disp('done')
figure;
subplot(2,1,1);
tr_genetic.plot_tree
title(['Genetic algorithm. Error ',num2str(kantoro_genetic)])
subplot(2,1,2);
tr_milp.plot_tree
title(['Stage-wise MILP algorithm. Error ',num2str(kantoro_milp)])
