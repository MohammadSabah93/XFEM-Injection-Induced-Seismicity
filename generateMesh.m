function [connectivity_global, node_coordinates_global, nNodes, nElements, nodesPerElement, elementType] = generateMesh()
   
   % Generate 2D mesh
    domain_length = 10;
    domain_width =0.5;
    nex = 80;
    ney = 20;
    
    [connectivity_global, node_coordinates_global, nNodes, nElements] = MeshGeneration_2D(domain_length, domain_width, nex, ney);

    nodesPerElement = size(connectivity_global, 2);

    if nodesPerElement == 4
        elementType = 'QUAD4';
    else
        elementType = 'TRI3';
    end
end