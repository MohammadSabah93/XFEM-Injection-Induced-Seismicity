function [gp, w] = gaussPoints(elementType, order)
    % Function to compute Gauss points and weights for 1D, 2D triangular,
    % and 2D quadrilateral elements.
    %
    % Inputs:
    %   elementType - '1D', 'TRI3', 'QUAD4'
    %   order - order of quadrature (1, 2, or higher depending on requirement)
    %
    % Outputs:
    %   gp - Gauss points (matrix where each row represents a Gauss point)
    %   w  - Weights corresponding to each Gauss point
    
    switch elementType
        case '1D'
            [gp, w] = gauss1D(order);
            
        case 'TRI3'
            [gp, w] = gauss2DTriangle(order);
            
        case 'QUAD4'
            [gp, w] = gauss2DQuadrilateral(order);
            
        otherwise
            error('Invalid element type. Choose from ''1D'', ''TRI3'', or ''QUAD4''.');
    end
end

%% Gauss points for 1D element
function [gp, w] = gauss1D(order)
    % Compute Gauss points and weights for 1D elements.
    
    % Use built-in function for standard Gauss-Legendre quadrature in 1D
    [gp, w] = gaussLegendre1D(order);
end

%% Gauss points for 2D triangular element
function [gp, w] = gauss2DTriangle(order)
    % Gauss points and weights for triangular element (on the reference triangle)
    switch order
        case 1
            gp = [1/3, 1/3]; % Barycentric coordinates for a 1-point rule
            w  = 1/2;        % Weight for the 1-point rule (area = 1/2)
            
        case 3
            gp = [1/6, 1/6; 2/3, 1/6; 1/6, 2/3];
            w  = [1/3; 1/3; 1/3];
            
        case 6
            gp = [0.0915762135, 0.0915762135;
                  0.8168475730, 0.0915762135;
                  0.0915762135, 0.8168475730;
                  0.4459484909, 0.1081030182;
                  0.1081030182, 0.4459484909;
                  0.4459484909, 0.4459484909];
            w  = [0.1099517437; 0.1099517437; 0.1099517437;
                 0.2233815897; 0.2233815897; 0.2233815897];
        case 7
            gp(1,:) = [0.101286507323456 0.101286507323456];
            gp(2,:) = [0.470142064105115 0.059715871789770];
            gp(3,:) = [0.797426985353087 0.101286507323456];
            gp(4,:) = [0.333333333333333 0.333333333333333];
            gp(5,:) = [0.059715871789770 0.470142064105115];
            gp(6,:) = [0.470142064105115 0.470142064105115];
            gp(7,:) = [0.101286507323456 0.797426985353087];

            w(1,:) = 0.125939180544827;
            w(2,:) = 0.132394152788506;
            w(3,:) = 0.125939180544827;
            w(4,:) = 0.225030000300000;
            w(5,:) = 0.132394152788506;
            w(6,:) = 0.132394152788506;
            w(7,:) = 0.125939180544827;
             case 9
            gp(1,:) = [0.165409927389841 0.037477420750088];
            gp(2,:) = [0.797112651860071 0.037477420750088];
            gp(3,:) = [0.437525248383384 0.124949503233232];
            gp(4,:) = [0.037477420750088 0.165409927389841];
            gp(5,:) = [0.797112651860071 0.165409927389841];
            gp(6,:) = [0.124949503233232 0.437527248383384];
            gp(7,:) = [0.437525248383384 0.437525248383384];
            gp(8,:) = [0.037477420750088 0.797112651860071];
            gp(9,:) = [0.165409927389841 0.797112651860071];

            w(1,:) = 0.063691414286223;
            w(2,:) = 0.063691414286223;
            w(3,:) = 0.205950504760887;
            w(4,:) = 0.063691414286223;
            w(5,:) = 0.063691414286223;
            w(6,:) = 0.205950504760887;
            w(7,:) = 0.205950504760887;
            w(8,:) = 0.063691414286223;
            w(9,:) = 0.063691414286223;
        case 12
            gp(1,:)  = [0.063089014491502 0.063089014491502];
            gp(2,:)  = [0.310352451033785 0.053145049844816];
            gp(3,:)  = [0.636502499121399 0.053145049844816];
            gp(4,:)  = [0.873821971016996 0.063089014491502];
            gp(5,:)  = [0.249286745170910 0.249286745170910];
            gp(6,:)  = [0.501426509658179 0.249286745170910];
            gp(7,:)  = [0.053145049844816 0.310352451033785];
            gp(8,:)  = [0.636502499121399 0.310352451033785];
            gp(9,:)  = [0.249286745170910 0.501426509658179];
            gp(10,:) = [0.053145049844816 0.636502499121399];
            gp(11,:) = [0.310352451033785 0.636502499121399];
            gp(12,:) = [0.063089014491502 0.873821971016996];

            w(1,:)  = 0.050844906370207;
            w(2,:)  = 0.082851075618374;
            w(3,:)  = 0.082851075618374;
            w(4,:)  = 0.050844906370207;
            w(5,:)  = 0.116786275726379;
            w(6,:)  = 0.116786275726379;
            w(7,:)  = 0.082851075618374;
            w(8,:)  = 0.082851075618374;
            w(9,:)  = 0.116786275726379;
            w(10,:) = 0.082851075618374;
            w(11,:) = 0.082851075618374;
            w(12,:) = 0.050844906370207;
        case 13
            gp(1,:)  = [0.065130102902216 0.065130102902216];
            gp(2,:)  = [0.312865496004875 0.048690315425316];
            gp(3,:)  = [0.638444188569809 0.048690315425316];
            gp(4,:)  = [0.869739794195568 0.065130102902216];
            gp(5,:)  = [0.048690315425316 0.312865496004875];
            gp(6,:)  = [0.260345966079038 0.260345966079038];
            gp(7,:)  = [0.479308067841923 0.260345966079038];
            gp(8,:)  = [0.638444188569809 0.312865496004875];
            gp(9,:)  = [0.333333333333333 0.333333333333333];
            gp(10,:) = [0.260345966079038 0.479308067841923];
            gp(11,:) = [0.048690315425316 0.638444188569809];
            gp(12,:) = [0.312865496004875 0.638444188569809];
            gp(13,:) = [0.065130102902216 0.869739794195568];

            w(1,:)  =  0.053347235608839;
            w(2,:)  =  0.077113760890257;
            w(3,:)  =  0.077113760890257;
            w(4,:)  =  0.053347235608839;
            w(5,:)  =  0.077113760890257;
            w(6,:)  =  0.175615257433204;
            w(7,:)  =  0.175615257433204;
            w(8,:)  =  0.077113760890257;
            w(9,:)  = -0.149570044467670;
            w(10,:) =  0.175615257433204;
            w(11,:) =  0.077113760890257;
            w(12,:) =  0.077113760890257;
            w(13,:) =  0.053347235608839;
            
        otherwise
            error('Quadrature order not supported for triangular elements.');
    end
end

%% Gauss points for 2D quadrilateral element
function [gp, w] = gauss2DQuadrilateral(order)
    % Gauss points and weights for quadrilateral element
    [gp1D, w1D] = gaussLegendre1D(order);
    
    % Compute 2D Gauss points by taking Cartesian product of 1D points
    [x, y] = meshgrid(gp1D, gp1D);
    gp = [x(:), y(:)];
    
    % Compute the weights as the product of 1D weights
    [wx, wy] = meshgrid(w1D, w1D);
    w = wx(:) .* wy(:);
end

%% Gauss-Legendre quadrature in 1D
function [gp, w] = gaussLegendre1D(order)
    % Compute Gauss-Legendre points and weights for a given order
    switch order
        case 1
            gp = 0;
            w = 2;
        case 2
            gp = [-1/sqrt(3); 1/sqrt(3)];
            w = [1; 1];
        case 3
            gp = [-sqrt(3/5); 0; sqrt(3/5)];
            w = [5/9; 8/9; 5/9];
        case 4
            gp = [-sqrt((3 + 2*sqrt(6/5))/7); -sqrt((3 - 2*sqrt(6/5))/7);
                   sqrt((3 - 2*sqrt(6/5))/7);  sqrt((3 + 2*sqrt(6/5))/7)];
            w = [(18 - sqrt(30))/36; (18 + sqrt(30))/36; (18 + sqrt(30))/36; (18 - sqrt(30))/36];
        otherwise
            error('Order not supported for 1D Gauss quadrature.');
    end
end