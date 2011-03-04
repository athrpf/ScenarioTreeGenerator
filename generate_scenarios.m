function [xi, p] = generate_scenarios(n_scen, n_stages, distribution, ...
                                      varargin)
% Generates a matrix of n_scen scenarios in n_stages stages with
% for the desired distribution.
% 
% Usage:
% [xi, p] = generate_scenarios(n_scen, n_stages, distribution, varargin)
% 
% xi : scenario matrix (nscen times nstages)
% p  : probability of each scenario (nscen,1)
% distribution is a string
% varargin captures all information necessary for the distribution.


switch distribution
  case 'lognormal'
    % varargin for this distribution:
    % starting_value, mu, sigma
    starting_value = varargin{1};
    mu = varargin{2};
    sigma = varargin{3};
    xi = mu+sigma*randn(n_scen, n_stages);
    xi(:,1) = starting_value.*exp(xi(:,1));
    for kk=2:n_stages
        xi(:,kk) = xi(:,kk-1).*exp(xi(:,kk));
    end
    p = 1/n_scen*ones(n_scen,1);
  case 'normal_independent'
    xi = randn(n_scen, n_stages);
    p = 1/n_scen*ones(n_scen,1);
  otherwise
    xi = [];
    p = [];
    disp('Distribution was not recognized')
end

