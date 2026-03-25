% Written By: Matthew Jon Pais, University of Florida (2009)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [enrElem] = enrElem(CONNEC, NODES)
% This function finds the enriched elements.

% Find the elements consisting of the crack tip domain
I = find(NODES(:,4) > 0);                                                   % Nodes containing crack tip enrichment
[a,b,c1]    = intersect(I,CONNEC(:,1)');                                    % Find elements with crack tip enrichment
[a,b,c2]    = intersect(I,CONNEC(:,2)');                                    % Find elements with crack tip enrichment
[a,b,c3]    = intersect(I,CONNEC(:,3)');                                    % Find elements with crack tip enrichment
[a,b,c4]    = intersect(I,CONNEC(:,4)');                                    % Find elements with crack tip enrichment
CTElem      = unique([c1; c2; c3; c4]);                                        % Elements containing crack tip enrichment

                                                
% Find the elements consiting of the crack body domain
I = find(NODES(:,2) > 0);                                               % Nodes containing crack tip enrichment
[a,b,c1]    = intersect(I,CONNEC(:,1)');                                % Find elements with crack tip enrichment
[a,b,c2]    = intersect(I,CONNEC(:,2)');                                % Find elements with crack tip enrichment
[a,b,c3]    = intersect(I,CONNEC(:,3)');                                % Find elements with crack tip enrichment
[a,b,c4]    = intersect(I,CONNEC(:,4)');                                % Find elements with crack tip enrichment
HElem       = unique([c1; c2; c3; c4]);                                    % Elements containing crack tip enrichment
    
% Remove duplicate entries
enrElem     = unique([CTElem; HElem]);

