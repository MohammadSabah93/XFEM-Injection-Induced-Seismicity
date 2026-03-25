function [jump] = Interface2(connectivity, node_coordinates, NODES, CRACK, a, elementType)

counter = 1;
nElements = size(connectivity,1);
nNodes = size (node_coordinates,1);

% Initialze vectors
[~,fullenr] = enrElem(connectivity,NODES);

jump=zeros(fullenr,2);    % Stiffness matrix

for iElem = 1:nElements
        
elementNodes = connectivity(iElem, :)';
nodesPerElement =size(connectivity(iElem, :),2);
Xe = [node_coordinates(connectivity(iElem, :), 1), node_coordinates(connectivity(iElem, :), 2)];
CTN = nnz(NODES(elementNodes,4));        % Number of nodes with crack tip enrichment    
HEN = nnz(NODES(elementNodes,2));        % Number of nodes with Heaviside enrichment
NEN = HEN+CTN;                           % Number of enriched nodes

if (HEN == 4)

[intersections, ~] = findCrackIntersections(CRACK, Xe);
[gp,~] = gaussPointsOnCrack(intersections, Xe, 2, elementType);

     for gpIndex = 1:size(gp, 1)

     [N, ~] = ShapeFunction(gp(gpIndex, :), elementType);

     index = 1;
      Benr = [];
      iLoc = 1;
      for iN = 1:nodesPerElement
        if NODES(elementNodes(iN),2) ~= 0

          Ba = [N(1,iN)     0;
                0    N(1,iN)];
          
          Benr(:,index:(index+1)) = Ba;
          index = index+2;

          DOF(iLoc:(iLoc+1)) = [2*NODES(elementNodes(iN),2)-1 2*NODES(elementNodes(iN),2)] - 2*nNodes;
          iLoc = iLoc+2;

        end
      end
     dis_jump(:,gpIndex) = 2 * Benr * a(DOF);
     end

jump(counter,:) = mean(dis_jump,2)';
counter=counter+1;
     
end
end
end