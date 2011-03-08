function tr = backwardtreeKmediods(xi, p, n_children)
% Very fast computation of a tree using kmeans++ and kmedioids
% works with very large sets of scenarios
% USAGE
% tr = backwardtreeKmediods(xi, p, n_children)
% xi : n_stages x n_scen
% 
% TODO
% - neue Scenariowerte als Mittelwerte ALLER URSPRUENGLICHEN SCENARIEN
% - Wahrscheinlichkeiten ausgeben

[n_stages, n_scen] = size(xi);
K = n_children^n_stages
n_treescen = K;
t = zeros(n_stages, n_treescen);

C = cell(n_stages,1);
L = cell(n_stages,1);

[L{1},C{1}] = kmeans(xi, K);
t(end,:) = C{1}(end,:);
%keyboard
for kk=2:n_stages
    keyboard
    K = K/n_children;
    [L{kk},Cidx] = kmedioids(squareform(pdist(C{kk-1}')), K);
    C{kk} = C{kk-1}(1:n_stages-kk+1,Cidx);
    idx = L{2};
    for jj=3:kk
        idx = L{jj}(idx);
    end
    t(n_stages-kk+1,:) = C{kk}(n_stages-kk+1,idx);
    %keyboard
end

