function [globalF_total] = FluidForce (params, connectivity, XYZ, PSI, NODES, CRACK, t, W, elementType, xs, qs, P, nt, dt)


% Extract material properties
mu=params.mu;
rho_f=params.rho_f;
com_f=params.com_f;
g = params.g;        


[~,DOF_p,~] = calcDOF(NODES);
globalH = zeros(DOF_p,1); 
globalS = zeros(DOF_p,1); 
globalF = zeros(DOF_p,1); 
nElements = size(connectivity,1);
c = 0;
% Flag to ensure injection flux is added only once
injectionAdded = false;

for iElem = 1:nElements
        
nodesPerElement =size(connectivity(iElem, :),2);
Xe = [XYZ(connectivity(iElem, :), 1), XYZ(connectivity(iElem, :), 2)];
N1  = connectivity(iElem,1);                                                  % Node 1 for current element
N2  = connectivity(iElem,2);                                                  % Node 2 for current element
N3  = connectivity(iElem,3);                                                  % Node 3 for current element
N4  = connectivity(iElem,4);                                                  % Node 4 for current element
NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element    
HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment

local = connectivity(iElem, :);                                                                                  % Initialize stiffness for current element
iLoc = 5;
localH = 0;
localS = 0;
localF = 0;

if (HEN ~= 0)

[intersections, jacobian] = findCrackIntersections(CRACK, Xe);

 if isempty(intersections)
        continue;  % Skip this element if no crack intersection is found
 end
 if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [PSI(N1);  PSI(N2);  PSI(N3);  PSI(N4)];                 % Nodal crack level set values
 end
[gp,gw] = gaussPointsOnCrack2(intersections, Xe, 2, elementType);
c = c + 1;

     for gpIndex = 1:size(gp, 1)

     [N, dNLocal] = ShapeFunction(gp(gpIndex, :), elementType);
     J= (dNLocal * Xe);
     dNGlobal = J \ dNLocal;

     enr = find(NN(:, 2) ~= 0)';
     psi_modi = N(1,enr) * abs(PN(enr,1))- abs(N(1,enr) * PN(enr,1));
     dpsi_dx = dNGlobal(1,enr) * abs(PN(enr,1)) - dNGlobal(1,enr) * PN(enr,1) * sign (N(1,enr) * PN(enr,1));
     dpsi_dy = dNGlobal(2,enr) * abs(PN(enr,1)) - dNGlobal(2,enr) * PN(enr,1) * sign (N(1,enr) * PN(enr,1));
     dpsi = [dpsi_dx;dpsi_dy];  

     index = 1;
     Benr = [];
     Nenr = [];
     for iN = 1:nodesPerElement
        if NODES(elementNodes(iN),2) ~= 0
       
           Benr(:,index) = dNGlobal(:,iN) * psi_modi + dpsi * N(iN);
           Nenr(:,index) = N(iN) * psi_modi;
           index = index+1;

          if (gpIndex == length(gp))
               local(iLoc) = NN(iN,2);
               iLoc = iLoc+1;
          end
        end
     end

     B1=[dNGlobal, Benr];
     B2=[N, Nenr];
     localH = localH - gw(gpIndex) * (t'*B1)' * W(c)^3/(12*mu) * (t'* B1) * det(jacobian);
     localS = localS - gw(gpIndex) * B2' * W(c) * com_f * B2 * det(jacobian);
     localF = localF + gw(gpIndex) * (t'*B1)' * W(c)^3/(12*mu) * rho_f * (t' * g) * det(jacobian);
    
     end

    globalH(local,1) = globalH(local,1) + localH * P(local,nt+1); 
    globalS(local,1) = globalS(local,1) + localS * (P(local,nt+1)-P(local,nt))/dt; 
    globalF(local,1) = globalF(local,1) + localF; 

 % Injection flux: Add injection flux from qs only once if xs is within this element
            if  ~isempty(xs) && isBetweenIntersections(xs, intersections, 1e-6) && ~injectionAdded 
                [xi] = global2LocalCoord2D(xs, Xe, elementType);
                [Nstd, ~] = ShapeFunction(xi, elementType);
                psi_modi_xs = Nstd * abs(PSI(elementNodes)) - abs(Nstd * PSI(elementNodes));
                index  = 1;
                Benr_xs = [];
                for iN = 1:nodesPerElement
                    if NODES(elementNodes(iN),2) ~= 0
                        N_enr_xs = Nstd(iN) * psi_modi_xs;
                        Benr_xs(:, index) = N_enr_xs;
                        index = index + 1;
                    end
                end
                B_xs = [Nstd, Benr_xs];
                localFs = B_xs' * W(c) * qs;
                globalF(local, 1) = globalF(local, 1) + localFs;
                injectionAdded = true;
            end
end
end

globalF_total = globalH +  globalS + globalF;
end


function [xi] = global2LocalCoord2D(x_s, Xe, elementType)
    
    % Initialize
    tol = 1e-8;
    maxIter = 10;
    xi = [0, 0];  % guess

    for iter = 1:maxIter
        [N, dNLocal] = ShapeFunction(xi, elementType);
        x_current = N * Xe;  % approximate global coords
        R = x_s - x_current;
        if norm(R) < tol
            break;
        end
        % Compute the Jacobian dx/dxi
        J = dNLocal * Xe;   % 2x2 for a quad
        dxi = J \ R';
        xi  = xi + dxi';
    end
end

function between = isBetweenIntersections(xs, intersections, tol)
    % Checks if xs lies between the intersection points of the crack with the element edges.
    % Assumes intersections contains at least two points.
    if size(intersections, 1) < 2
        between = false;
        return;
    end
    % For simplicity, take the first two intersection points.
    A = intersections(1, :);
    B = intersections(2, :);
    % Compute the vector from A to B.
    AB = B - A;
    % Compute the projection parameter of xs onto AB.
    t = dot(xs - A, AB) / dot(AB, AB);
    % Compute the projection of xs on the line AB.
    projection = A + t * AB;
    dist = norm(xs - projection);
    % xs is between if t is between 0 and 1 and its distance to the line is less than tol.
    between = (t >= 0 && t <= 1 && dist < tol);
end