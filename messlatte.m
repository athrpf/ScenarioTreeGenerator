function m = messlatte(nstages, nchildren, distribution, varargin)
% Computes the value of the Kantorovich Distance that a random
% selection (CE-case) would achieve.
% USAGE
% function m = messlatte(nstages, nchildren, distribution, varargin)


disp('Generating initial scenarios')
nscen = nstages^nchildren;
[xi,p] = generate_scenarios(nscen,nstages, distribution, varargin{1}, varargin{2},varargin{3});

% generate tree from new or old (depending on ce-de) set of
% scenarios and nchildren
disp('Constructing most basic tree')
tr = scen_to_tree(xi,p, nchildren);
tr.p = 1/tr.n_scenarios*ones(tr.n_scenarios,1);

ntestscen = 3000;

disp(['Starting test with ', num2str(ntestscen),' scenarios'])
[nu,q] = generate_scenarios(ntestscen, nstages, distribution,varargin{1}, varargin{2},varargin{3});

m = tr.kantorovich(nu, q, 2);
