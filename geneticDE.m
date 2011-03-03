function tr = geneticDE(xi, p, n_children)
% Hans Pirnay
% 2011-03-03
% Computes a tree for a discrete-event stochastic process using
% genetic optimization
% (the algorithm is only implemented for one dimensional
% stochastic processes)

% compute the distance vector c_ijtd

    I = size(xi,1);
    T = size(xi,2);
    n_var = -1 + (1-n_children^(T+1))/(1-n_children);
    
    % create distance matrix
    c = zeros(I,I,T);
    for t=1:T
        for j=1:I
            for i=1:j-1
                c(i,j,t) = abs(xi(i,t) - xi(j,t));
            end
        end
        c(:,:,t) = c(:,:,t) + c(:,:,t)';
    end
    
    
    options = gaoptimset('CreationFcn', @(GenomeLength, FitnessFcn, options) createDE(GenomeLength, FitnessFcn, options, I) , ...
                         'CrossoverFcn', @crossoverscattered, ...
                         'MutationFcn', @(parents, options, nvars, ...
                                          FitnessFcn, state, thisScore, ...
                                          thisPopulation) ...
                         mutateDE(parents, options, nvars, FitnessFcn, state, thisScore, thisPopulation,I), ...
                         'PlotFcns', {@gaplotscores,@gaplotbestf},...
                         'PopulationType', 'Custom');
    z = ga(@(Z) geneticDEfitness(Z,c,p,n_children),n_var,[],[],[],[],[],[],[],options);

    % create tree from integer solution
    tr = tree(T+1, n_children, true);
    nps = n_children;
    idx = 1;
    for t=1:T
        tr.node_values(idx+1:idx+nps) = xi(z(idx:idx+nps-1),t);
        idx = idx+nps;
        nps = nps*n_children;
    end

    opt_diff = tr.compute_optimal_weights(xi, p, 1)
end


function fval = geneticDEfitness(z, c, p, n_children) 
% Computes the Kantorovich Distance for a fixed set of selected scenarios
% fval = geneticDEfitness(z, c, p, n_children, n_stages) 
% 
% z : binary vector of assignments 
    I = size(c,1);
    T = size(c,3);
    J = n_children^T;
    
    fval = 0;
    z_first_idx = 1;
    z_last_idx = n_children;
    z_n = n_children;
    nps = J/n_children; % nodes per stage
    cv = zeros(I,J);
    for t=1:T
        %fval = fval + sum(p.*min(c(:,z(z_first_idx:z_last_idx),t),[], 2));
        for j=1:z_n
            cv(:,1+nps*(j-1):nps*j) = bsxfun(@plus, cv(:,1+nps*(j-1):nps*j), c(:,z(z_first_idx+j-1),t));
        end
        z_first_idx = z_last_idx+1;
        z_n = z_n*n_children;
        z_last_idx = z_first_idx+z_n-1;
        nps = nps/n_children;
    end
    fval = sum(p.*min(cv,[],2));
end

function created = createDE(GenomeLength, FitnessFcn, options, n_scenarios) 
    npop = 100;
    created = ceil(rand(100,GenomeLength)*n_scenarios);
end

function mutated = mutateDE(parents, options, nvars, FitnessFcn, ...
                            state, thisScore, thisPopulation, ...
                            n_scenarios)
    nparents = length(parents);
    
    mutprob = 0.5;
    mutated = zeros(nparents, nvars);
    r = round(randn(nparents, nvars))*4;
    
    mutated = thisPopulation(parents,:) + r;
    mutated(mutated<1) = 1;
    mutated(mutated>n_scenarios) = n_scenarios;
    
end
