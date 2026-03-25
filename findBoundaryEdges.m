function boundaryEdges = findBoundaryEdges(elementConnectivity,boundaryNodes)
    % This function finds the boundary edges of elements that lie on a specific boundary.
    % It checks which edges of the elements have nodes on the boundary.
    %
    % Inputs:
    %   elementConnectivity - Element connectivity matrix (nElements x nodesPerElement)
    %   nodeCoordinates     - Coordinates of all nodes (nNodes x 2 matrix)
    %   boundaryNodes       - Indices of nodes on the boundary
    %   tolerance           - A small value to determine if a node is on the boundary
    %
    % Output:
    %   boundaryEdges       - Array of boundary edges with format [elementIndex, node1, node2]
    
    nElements = size(elementConnectivity, 1);  % Number of elements
    boundaryEdges = [];  % To store the element index and boundary edge nodes

    % Loop through each element
    for i = 1:nElements
        elementNodes = elementConnectivity(i, :);  % Nodes of the current element
        nNodesPerElement = length(elementNodes);   % Number of nodes in the element

        % Get all edges of the element
        if nNodesPerElement == 3  % For triangular elements
            edges = [elementNodes(1), elementNodes(2);  % Edge 1-2
                     elementNodes(2), elementNodes(3);  % Edge 2-3
                     elementNodes(3), elementNodes(1)]; % Edge 3-1
        elseif nNodesPerElement == 4  % For quadrilateral elements
            edges = [elementNodes(1), elementNodes(2);  % Edge 1-2
                     elementNodes(2), elementNodes(3);  % Edge 2-3
                     elementNodes(3), elementNodes(4);  % Edge 3-4
                     elementNodes(4), elementNodes(1)]; % Edge 4-1
        else
            error('Unsupported element type: must be TRI3 or QUAD4.');
        end

        % Check each edge to see if it lies on the boundary
        for j = 1:size(edges, 1)
            edgeNodes = edges(j, :);
            if all(ismember(edgeNodes, boundaryNodes))
                % If both nodes of the edge are in boundaryNodes, it's a boundary edge
                boundaryEdges = [boundaryEdges; edgeNodes];  % Store element index and edge nodes
            end
        end
    end
end
