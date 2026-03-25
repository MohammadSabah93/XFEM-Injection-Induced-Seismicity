function fe_Insitu = InsituForce(Xe, params, elementType, nodesPerElement)
    % Computes the in-situ force vector for a 2D element
    % Inputs:
    %   Xe              - Coordinates of the element nodes (matrix)
    %   params          - Structure containing stress and pressure components:
    %                     Sxxr: In-situ horizontal stress in x-direction
    %                     Syyr: In-situ horizontal stress in y-direction
    %                     Sxyr: In-situ shear stress
    %                     P0: Initial pore pressure
    %   elementType     - Type of the element (e.g., 'QUAD4', 'TRI3')
    %   nodesPerElement - Number of nodes in the element
    % Output:
    %   fe_Insitu       - Element in-situ force vector (2*nodesPerElement x 1)

    % Extract in-situ stresses and pore pressure from parameters
    Sxxr = params.Sxxr;  % In-situ horizontal stress in x-direction
    Syyr = params.Syyr;  % In-situ horizontal stress in y-direction
    Sxyr = params.Sxyr;  % In-situ shear stress
    P0 = params.P0;      % Initial pore pressure

    % Initialize the in-situ force vector
    fe_Insitu = zeros(2 * nodesPerElement, 1);

    % Get Gauss points and weights for numerical integration
    [gp, weights] = gaussPoints(elementType, 2);

    % Loop over Gauss points to compute the in-situ force vector
    for gpIndex = 1:size(gp, 1)
        % Compute shape function derivatives in the local coordinate system
        [~, DNL] = ShapeFunction(gp(gpIndex, :), elementType);

        % Compute the Jacobian matrix and its determinant
        J = DNL * Xe;
        detJ = det(J);

        % Compute the shape function derivatives in the global coordinate system
        DNG = J \ DNL;

        % Construct the strain-displacement matrix (B matrix)
        Be = zeros(3, 2 * size(Xe, 1));
        Be(1, 1:2:end) = DNG(1, :);  % ∂N/∂x for x-direction
        Be(2, 2:2:end) = DNG(2, :);  % ∂N/∂y for y-direction
        Be(3, 1:2:end) = DNG(2, :);  % ∂N/∂y for x-direction (shear)
        Be(3, 2:2:end) = DNG(1, :);  % ∂N/∂x for y-direction (shear)

        % Assemble the in-situ force vector
        fe_Insitu = fe_Insitu + weights(gpIndex) * Be' * [Sxxr + P0; Syyr + P0; Sxyr] * detJ;
    end
end
