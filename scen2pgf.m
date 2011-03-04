function scen2pgf(xi, init)
% Takes a stochastic process in scenario form and prints it to stdout
% Usage:
% function scen2pgf(xi, init)
% xi : of size n_stages x n_scen
% init : optional; sets an initial value that is common to all scenarios

[n_stages, n_scen] = size(xi);

if n_stages>10 || n_scen>20
    disp('are you sure you wanna do this??')
end

if nargin>1
    xi = [init*ones(1,n_scen);xi];
end

for scen = xi
    disp('\addplot [color=black] coordinates {')
    for k=1:n_stages
        disp(['(',num2str(k),', ', num2str(scen(k)),')'])
    end
    disp('};')
end
