function [tr,node_ptrs] = treeDEsw(xi_in, p, n_children, norm_type)
% Computes a tree in the stagewise fashion. Will consider scenarios
% as fixed (notion of 'discrete events')
% 
% USAGE
% function tr = treeDEsw(xi)
% 
% tr : the tree output as nodewise
% node_ptrs : cell array of all scenarios of xi that belong to each node
% xi : the scenario matrix (n_stages,n_scen)
% n_children : number of children / branchings of each node 
% norm_type  : the norm in which distances are to be computed. will
% be 2 if left empty

    if nargin < 3
        norm_type = 2;
    end
    
    xi = [];
    n_stages = 0;
    n_scen = length(p);
    if n_scen==size(xi_in,2)
        xi = xi_in;
        n_stages = size(xi_in,1);
    elseif n_scen==size(xi_in,1)
        xi = xi_in';
        n_stages = size(xi_in,2);
    else
        disp(['number of scenarios in p does not match any dimension ' ...
              'in xi']);
        return
    end

    %n_nodes = -1 + (1-n_children^(n_stages+1))/(1-n_children);
    tr = tree(n_stages+1, n_children, true);
    tr.node_values(1) = 0;

    n_nodes = tr.n_nodes;
    n_fathering_nodes = n_nodes - tr.n_scenarios;
    node_ptrs = cell(n_nodes,1); % mapping from each node to the
                                 % scenarios that map to it
    node_ptrs{1,1} = 1:n_scen;
    eta_ptrs = cell(n_nodes,1);
    eta_ptrs{1,1} = p;

    for node=1:n_fathering_nodes
        
        %[parent_idx,parent_stage] = tr.parent(node);
        stage = tr.stage_of_node(node);
        indices = node_ptrs{node,1};
        n_fatherscen = length(indices);
        % compute cost matrix
        c = zeros(n_fatherscen);
        for kk=1:n_fatherscen
            for jj=1:kk-1
                c(kk,jj) = norm(xi(stage,indices(kk))-xi(stage,indices(jj)), norm_type);
            end
        end
        c = c+c';
        
        A = treeDEconstraint_matrix(n_fatherscen);
        n_con = size(A,1);
        b = zeros(n_con,1);
        b(1:n_fatherscen) = eta_ptrs{node,1};
        b(n_con) = n_children;
        contypes = [repmat('=',1,n_fatherscen),repmat('<',1, ...
                                                      n_fatherscen^2+1)];
        vartypes = [repmat('C',1,n_fatherscen^2),repmat('B',1,n_fatherscen)];
        
        [xopt,valopt,flag] = gurobi_mex([c(:);zeros(n_fatherscen,1)], ...
                                        1, A, b, contypes, [],[], ...
                                        vartypes);
        if flag~=2
            disp(['Problem encountered while solving MILP with GUROBI: flag',flag]); 
        end
        
        %eta_eps = 1e-11;
        etaopt = zeros(n_fatherscen);
        etaopt(:) = xopt(1:n_fatherscen^2);
        zopt = xopt(n_fatherscen^2+1:end);
        
        tmp = sum(etaopt);
        
        children = tr.children(node)
        selections = find(zopt>0.5);
        if length(selections)<n_children
            selections = [selections;selections(end)* ...
                          ones(n_children-length(selections),1)];
            %keyboard
        end
        
        for child_idx=1:length(selections)
            indices_this_notation = find(etaopt(selections(child_idx),:)>0);
            node_ptrs{children(child_idx)} = node_ptrs{node}(indices_this_notation);
            eta_ptrs{children(child_idx)}  = etaopt(selections(child_idx),indices_this_notation);
        end
        tr.node_values(children) = xi(stage,node_ptrs{node}(selections));
        
        %keyboard
    end
    
    tr.compute_optimal_weights(xi, p, norm_type)
end

