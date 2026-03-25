% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [DOF_u,DOF_p,DISP] = calcDOF(NODES)
% This function calculates the total degrees of freedom consisting of
% classical and enriched DOF.

heavDOF = 2*nnz(NODES(:,2));                                                % Define the number of Heavi DOF
ctipDOF = 8*nnz(NODES(:,4));                                                % Define the number of crack tip DOF
% if heavDOF > 0, heaviNodes;     end                                         % Define heaviside nodal enrichment values
% if ctipDOF > 0, ctipNodes;      end                                         % Define crack tip nodal enrichment values
DOF_u     = 2*max(max(NODES));                                                % Total number of degrees of freedom
DOF_p     = max(max(NODES(:,1)));  
DISP    = sparse(DOF_u,1);                                                    % Initialize displacement vector

end