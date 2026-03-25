function globalD = dampingBoundary(boundaryConditions, XYZ, connectivity, NODES, elementType, params)
% Assembles global boundary damping matrix using absorbing boundary conditions for any element shape

% Material and damping parameters
cs = params.cs;
cp = params.cp;
rho = params.rho_m;
as = 1;
ap = 0.5;

% Damping tensor
V = rho * (as * cs * eye(2) + ap * cp * eye(2));

% Initialization
tractionBoundaries = boundaryConditions.tractionBoundaries;
[globalDOF, ~, ~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * globalDOF;
globalD = spalloc(globalDOF, globalDOF,nnzEst);
nElements = size(connectivity, 1);

for iElem = 1:nElements
    elementNodes = connectivity(iElem, :);
    elementNodes = elementNodes(elementNodes > 0);  % remove zero-padding if any
    nodesPerElement = length(elementNodes);
    DOF = reshape([2 * elementNodes - 1; 2 * elementNodes], [], 1);

    for i = 1:size(tractionBoundaries, 1)
        [isInBoundary, indexInElement] = ismember(tractionBoundaries(i, 1:2), elementNodes);
        if all(isInBoundary)
            node1 = tractionBoundaries(i, 1);
            node2 = tractionBoundaries(i, 2);
            L = norm(XYZ(node1,:) - XYZ(node2,:));

            % Robust direction check based on coordinates
            if abs(XYZ(node1,1) - XYZ(node2,1)) < 1e-8  % horzontal edge (same x)
                V = rho * [as*cs, 0; 0, ap*cp];
            else  % horizontal or general edge
                V = rho * [ap*cp, 0; 0, as*cs];
            end
            edgeNodes = sort(indexInElement);
            edgeIdx = identifyEdgeIndex(nodesPerElement, edgeNodes);

            if isempty(edgeIdx)
                warning('Edge not identified for element %d', iElem);
                continue;
            end

            [gp, gw] = gaussPoints('1D', 2);
            localD = zeros(length(DOF));

            for j = 1:length(gp)
                xi_edge = mapEdgeGaussPoint(elementType, nodesPerElement, edgeIdx, gp(j));
                [N, ~] = ShapeFunction(xi_edge, elementType);

                Ne = zeros(2, 2 * nodesPerElement);
                Ne(1, 1:2:end) = N;
                Ne(2, 2:2:end) = N;

                localD = localD + gw(j) * Ne' * V * Ne * L / 2;
            end

            globalD(DOF, DOF) = globalD(DOF, DOF) + localD;
        end
    end
end
end

function edgeIdx = identifyEdgeIndex(numNodes, edgeNodes)
% General method to find index of edge defined by edgeNodes
if numNodes == 3
    edges = [1 2; 2 3; 1 3];
elseif numNodes == 4
    edges = [1 2; 2 3; 3 4; 1 4];
elseif numNodes == 6
    edges = [1 4; 4 2; 2 5; 5 3; 3 6; 6 1];
else
    edges = nchoosek(1:numNodes, 2);  % fallback to all pairs
end
edgeIdx = find(ismember(edges, edgeNodes, 'rows'), 1);
end

function xi_edge = mapEdgeGaussPoint(elementType, ~, edgeIdx, gaussPoint)
% Maps Gauss point to natural coordinates on an edge for arbitrary shape
switch lower(elementType)
    case 'tri3'
        switch edgeIdx
            case 1, xi_edge = [(gaussPoint + 1) / 2, 0];
            case 2, eta = (gaussPoint + 1) / 2; xi_edge = [1 - eta, eta];
            case 3, xi_edge = [0, (gaussPoint + 1) / 2];
        end
    case 'quad4'
        switch edgeIdx
            case 1, xi_edge = [gaussPoint, -1];
            case 2, xi_edge = [1, gaussPoint];
            case 3, xi_edge = [gaussPoint, 1];
            case 4, xi_edge = [-1, gaussPoint];
        end
    otherwise
        error('Unsupported element type or edge mapping not implemented.');
end
end
