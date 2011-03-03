% Hans Pirnay
% 2011-02-27

classdef tree < handle
    
    properties
        has_father  % if false this tree is actually a forest with
                    % n_children trees
        n_children
        n_stages
        n_scenarios
        n_nodes
        node_values
    end
    
    methods
        function obj = tree(n_stages, n_children, has_father)
            obj.has_father  = has_father;
            obj.n_stages    = n_stages;
            obj.n_children  = n_children;
            obj.n_nodes     = 0;
            if has_father
                obj.n_nodes = (1-n_children^n_stages)/(1- ...
                                                           n_children);
                obj.n_scenarios = n_children^(n_stages-1);
            else
                obj.n_nodes = -1 + (1-n_children^(n_stages+1))/(1- ...
                                                                n_children);
                obj.n_scenarios = n_children^n_stages;
            end
             node_values = zeros(obj.n_nodes,1);
        end
        
        function set_node(obj, stage, scenario, value)
            node_idx = scen_stage_to_node(obj,stage,scenario)
            obj.node_values(node_idx) = value;
        end
        
        function node_idx= stage_scen_to_node(obj, stage, scenario)
            nc = obj.n_children;
            if obj.has_father
                first_node_this_stage = 1+(1-nc^(stage-1))/ ...
                    (1-nc);
            else
                first_node_this_stage = (1-nc^(stage))/(1-nc);
            end
            scenarios_per_node = nc^(n_stages-stage);
            node_idx = first_node_this_stage + floor((scenario-1)/scenarios_per_node);
        end
        
        function stage = stage_of_node(obj, node)
        %returns the stage index of the node
            stage = 1;
            first_node_this_stage = 1;
            nodes_this_stage = obj.n_children;
            if obj.has_father
                nodes_this_stage = 1;
            end
            last_node_this_stage = first_node_this_stage+ ...
                nodes_this_stage-1;
            while node>last_node_this_stage
                stage = stage +1;
                first_node_this_stage = last_node_this_stage+1;
                nodes_this_stage = nodes_this_stage*obj.n_children;
                last_node_this_stage = first_node_this_stage+nodes_this_stage-1;
            end
        end
        
        function [parent_idx,parent_stage] =parent(obj,stage_or_node, scenario)
            node_idx = 0;
            if nargin < 3
                node_idx = stage_or_node;
            else
                node_idx = obj.stage_scen_to_node(stage,scenario);
            end
            
            nc = obj.n_children;
            stage = 1;
            first_node_this_stage = 1;
            nodes_per_stage = nc;
            if obj.has_father
                nodes_per_stage = 1;
            end
            last_node_this_stage = first_node_this_stage+nodes_per_stage-1;
            while last_node_this_stage<node_idx
                stage = stage +1;
                nodes_per_stage = nodes_per_stage*nc;
                first_node_this_stage = last_node_this_stage+1;
                last_node_this_stage = first_node_this_stage+nodes_per_stage-1;
            end
            node_in_stage_idx = node_idx-first_node_this_stage+1;
            parent_stage = stage -1;
            parent_in_stage_idx = floor((node_in_stage_idx-1)/nc);
            if obj.has_father
                first_node_parent_stage = 1+(1-nc^(parent_stage-1))/ ...
                    (1-nc);
            else
                first_node_parent_stage = (1-nc^(parent_stage))/(1-nc); 
            end
            parent_idx = first_node_parent_stage+parent_in_stage_idx;
        end
        
        function child_vec = children(obj, node) 
        % returns a vector of node indices pointing to the children
        % of node
        %nts = obj.nodes_this_stage(node);
            
            nc = obj.n_children;
            st  = obj.stage_of_node(node);
            if st == obj.n_stages
                child_vec = [];
                return
            end
            nts = obj.nodes_this_stage(st);
            fnts = obj.first_node_this_stage(st);
            fnns = obj.first_node_this_stage(st+1);
            child_vec = fnns+(node-fnts)*nc:fnns+(node-fnts+1)*nc-1;
        end
        
        function nts = nodes_this_stage(obj, stage)
            nts = obj.n_children^(stage-1);
            if ~obj.has_father
                nts = nts*obj.n_children;
            end
        end
        
        function fnts = first_node_this_stage(obj, stage)
            nc = obj.n_children;
            if obj.has_father
                fnts = (1-nc^(stage-1))/(1-nc)+1;
            else
                fnts = (1-nc^(stage))/(1-nc);
            end
        end
       
        function scenlist = tree2scen(obj)
            ns = obj.n_stages
            scenlist = zeros(ns, obj.n_scenarios);
            
            node_idx = 1;
            for stage=1:ns
                nts = obj.nodes_this_stage(stage);
                nspn = obj.n_scenarios/nts;
                for node=1:nts
                    scenlist(stage,(node-1)*nspn+1:node*nspn) = ...
                        obj.node_values(node_idx);
                    node_idx = node_idx+1;
                end
            end
        end
        
        function plot_tree(obj)
            scenlist = obj.tree2scen();
            plot(scenlist);            
        end
    end
    
end
