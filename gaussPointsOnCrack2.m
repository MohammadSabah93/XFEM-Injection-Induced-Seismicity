function [gaussPoint, w] = gaussPointsOnCrack2(intersections, Xe, order, elementType)
    % Compute Gauss points on a crack given by two physical intersection points.
    %
    % Inputs:
    %   intersections: [2 x 2] array with each row [x, y] for a crack endpoint.
    %   Xe: [4 x 2] nodal coordinates of the Quad4 element.
    %   order: Order of the 1D Gauss quadrature.
    %
    % Outputs:
    %   gaussPoint: Gauss points along the crack in (ksi, eta) parametric space.
    %   w: Quadrature weights.
    
    % Map each physical intersection point to parametric coordinates using a robust
    % inverse isoparametric mapping (works for points anywhere inside the element).
    parametricCoord = zeros(2, 2);
    for i = 1:2
       parametricCoord(i,:) = mapToParametricGeneral(intersections(i,:), Xe, elementType);
    end
    
    % Compute Gauss points along the crack in parametric space.
    [gaussPoint, w] = findGaussPointsOnCrack(parametricCoord, order);
end

function [xi_eta] = mapToParametricGeneral(point, Xe, elementType)
    % Robustly compute the parametric coordinates (xi, eta) for a given physical point
    % in a Quad4 element using Newton–Raphson iteration.
    %
    % Inputs:
    %   point: [1 x 2] physical coordinates [x, y].
    %   Xe: [4 x 2] nodal coordinates of the element.
    %
    % Outputs:
    %   xi_eta: [1 x 2] parametric coordinates [xi, eta].
    
    tol     = 1e-8;
    maxIter = 50;
    % Initial guess: if the point is inside the element, starting at the center ([0,0])
    xi     = 0;
    eta    = 0;
    
    for iter = 1:maxIter
        % Compute shape functions and their derivatives at (xi, eta)
        [N, DNL]=ShapeFunction([xi,eta],elementType);
        
        % Map to physical coordinates
        x_val = N * Xe(:,1);
        y_val = N * Xe(:,2);
        
        % Residual vector
        R = [x_val - point(1); y_val - point(2)];
        
        % Check for convergence
        if norm(R) < tol
            break;
        end
        
        % Compute the Jacobian matrix of the mapping
        J = [DNL(1,:) * Xe(:,1), DNL(2,:) * Xe(:,1);
             DNL(1,:) * Xe(:,2), DNL(2,:) * Xe(:,2)];
         
        % Update (xi, eta) by solving J * delta = -R
        delta = -J \ R;
        xi  = xi  + delta(1);
        eta = eta + delta(2);
    end
    
    if iter == maxIter
        warning('Newton-Raphson did not converge in mapToParametricGeneral');
    end
    
    xi_eta = [xi, eta];
end


function [gaussPoint, w] = findGaussPointsOnCrack(parametricCoord, order)
    % Compute Gauss points in the (xi, eta) parametric space along a crack.
    %
    % Inputs:
    %   parametricCoord: [2 x 2] array with rows [xi, eta] for the two crack endpoints.
    %   order: Order of the 1D Gauss quadrature.
    %
    % Outputs:
    %   gaussPoint: Gauss points along the crack in parametric coordinates.
    %   w: Gauss quadrature weights.
    
    % Retrieve 1D Gauss points and weights (assuming gaussPoints('1D', order) is defined)
    [gp, w] = gaussPoints('1D', order);
    
    % Extract parametric endpoints.
    xi1  = parametricCoord(1, 1); 
    eta1 = parametricCoord(1, 2);
    xi2  = parametricCoord(2, 1); 
    eta2 = parametricCoord(2, 2);
    
    numPoints  = length(gp);
    gaussPoint = zeros(numPoints, 2);
    
    % Parameterize the crack line in parametric space.
    for i = 1:numPoints
        % Correctly map Gauss points from [-1, 1] to the crack segment
        gaussPoint(i, 1) = (1 - gp(i))/2 * xi1 + (1 + gp(i))/2 * xi2;
        gaussPoint(i, 2) = (1 - gp(i))/2 * eta1 + (1 + gp(i))/2 * eta2;
    end
end
