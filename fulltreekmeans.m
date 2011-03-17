function tr = fulltreekmeans(xi, p, nc)

    [nscen,ns] = size(xi); % number of scenarios and number of stages
    nl   = nc^ns; % number of leaf nodes
    nxi  = size(xi,2);    
    nts  = nl; % number of tree scenarios = number of leaf nodes
    
    % initialize nu and r
    tr_init = scen_to_tree(xi, p, nc);
    nu = tr_init.node_values;
    r = zeros(1,nscen);
    
    % Convergence criteria
    converged = false;
    maxiter = 100;
    niter = 0;
    
    % initialize movie
    aviobj = avifile('emtree.avi');
    moviefig = figure('visible','off');
    
    while (~converged & niter <maxiter)
        niter = niter+1;
        sprintf('Iteration %d', niter)
        % Expectation step
        [r,convr] = compute_optimal_assignments(r, nu, xi, p, nc);
       
        % Maximization step
        [nu,convnu] = compute_maximization_step(r, nu, xi, p, nc);
        
        % test 
        tr_test = tree(ns+1, nc, true);
        tr_test.node_values(2:end) = nu;
        iter_error = tr_test.compute_optimal_weights(xi, p, 2);
        tr_test.plot_tree
        aviobj=addframe(aviobj,moviefig);
        sprintf('Error in this iteration: %f', iter_error)
        
        % Check for convergence
        converged = convr | convnu;
    end
    disp(['Done after ', num2str(niter),' iterations'])
    disp(['convr = ',num2str(convr),' convnu = ',num2str(convnu)])
    tr = tree(ns+1, nc, true);
    tr.node_values = nu;
    
    % end movie
    aviobj = close(aviobj);
    close(moviefig);
end

function nu = initialize_nodes(xi, ns, nc)
    nl = nc^ns; % number of leaf nodes
    nxi = size(xi,2);    
    
    nn = (1-nc^(ns+1))/(1-nc)
    nu = zeros(nn,1); % tree node values
            
    % seeding
    this_stage_start = 0;
    for kk=1:ns
        for jj=1:nc^kk
            nu(this_stage_start+jj) = xi(kk, randi(nxi));
        end
        this_stage_start = this_stage_start+nc^kk;
    end
    
end


function [r,convr] = compute_optimal_assignments(r, nu, xi, p, nc)
% The expectation step
    [nscen,ns] = size(xi);
    rold = r;
% Compute distances between nu and xi
    tr = tree(ns+1, nc, true);
    tr.node_values(2:end) = nu;
    tscen =tr.tree2scen';
    c = pdist2(tscen(:,2:end),xi);
    [minc, minidx] = min(c);
    r = minidx;
    convr = all(r==rold);
    %r = zeros(size(r));
    %r(minidx, minidx) = 1;
end

function [nu,convnu] = compute_maximization_step(r, nu, xi, p, nc)
    nuold = nu;
    tr = scen_to_tree(xi, p, nc, r);
    nu = tr.node_values;
    convnu = (norm(nu-nuold)<1e-10);
end
