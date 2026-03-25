function globalF = TractionForces(boundaryConditions, R, E, NODES, elementType)
    

tractionBoundaries = boundaryConditions.tractionBoundaries;
[globalDOF,~,~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * globalDOF;
globalF = spalloc(globalDOF, 1,nnzEst);
nElements = size(E,1);

for iElem = 1:nElements

elementNodes = E(iElem, :);
localF=0;

    for i = 1:size(tractionBoundaries, 1)

        % Check if the current element contains the boundary edge
        [isInBoundary, indexInElement] = ismember(tractionBoundaries(i, 1:2), elementNodes);

        if all(isInBoundary)  % If both nodes are in the current element

            nodesPerElement =size(E(iElem, :),2);
            DOF=reshape([2 * E(iElem, :) - 1; 2 * E(iElem, :)], [], 1);
            node1 = tractionBoundaries(i, 1);
            node2 = tractionBoundaries(i, 2);
            L = sqrt((R(node1, 1) - R(node2, 1))^2 + (R(node1, 2) - R(node2, 2))^2);  
            t = tractionBoundaries(i, 3:4)';  
            boundaryEdge = sort(indexInElement);

            edges = [1, 2; 2, 3; 1, 3];  
            if nodesPerElement == 4
                edges = [1, 2; 2, 3; 3, 4; 1, 4];  
            end
            
            edgeIdx = find(ismember(edges, boundaryEdge, 'rows'), 1);
            [gp, gw] = gaussPoints('1D', 2);

            
            for j = 1:size(gp, 1)
                edgeGaussPoints = getEdgeGaussPoints(edgeIdx, gp(j), nodesPerElement);
                [N, ~] = ShapeFunction(edgeGaussPoints, elementType);

                Ne = zeros(2, 2 * nodesPerElement);
                Ne(1, 1:2:end) = N;  % Shape functions for x-direction
                Ne(2, 2:2:end) = N;  % Shape functions for y-direction

                localF = localF + gw(j) * Ne' * t * L / 2;
            end
            globalF(DOF,1)=globalF(DOF,1)+localF;
        end
    end
end
end

% Helper function to map Gauss points to natural coordinates of an edge
function edgeGaussPoints = getEdgeGaussPoints(edgeIdx, gaussPoint, numNodes)

    if numNodes == 3  % Triangular elements
        switch edgeIdx
            case 1  % Edge 1-2
                edgeGaussPoints = [(gaussPoint + 1) / 2, 0];  % Map to [0, 1], eta = 0
            case 2  % Edge 2-3
                eta = (gaussPoint + 1) / 2;  % Map to [0, 1]
                edgeGaussPoints = [1 - eta, eta];  % xi + eta = 1
            case 3  % Edge 3-1
                edgeGaussPoints = [0, (gaussPoint + 1) / 2];  % Map to [0, 1], xi = 0
        end
    else  % Quadrilateral elements
        switch edgeIdx
            case 1
                edgeGaussPoints = [gaussPoint, -1];  % Edge 1 (xi varies, eta = -1)
            case 2
                edgeGaussPoints = [1, gaussPoint];  % Edge 2 (eta varies, xi = 1)
            case 3
                edgeGaussPoints = [gaussPoint, 1];  % Edge 3 (xi varies, eta = 1)
            case 4
                edgeGaussPoints = [-1, gaussPoint];  % Edge 4 (eta varies, xi = -1)
        end
    end
end
