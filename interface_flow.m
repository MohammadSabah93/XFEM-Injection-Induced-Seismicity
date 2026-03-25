function [globalH, globalS, globalQ, globalF] = interface_flow (params, connectivity, XYZ, PSI, NODES, CRACK, t, W, elementType)


% Extract material properties
mu=params.mu;
com_f=params.com_f;       
g = params.g;
rho_f = params.rho_f;

m = [1;1;0];

tx = t(1);
ty = t(2);
pt = [tx^2, ty^2, 2*tx*ty];

% Initialize the FE matrixes
[DOF_u,DOF_p,~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * (DOF_p+DOF_u);
nnzEst_p = avgBW * DOF_p;
globalQ = spalloc(DOF_p,DOF_u,nnzEst);
globalH = spalloc(DOF_p,DOF_p,nnzEst_p); 
globalS = spalloc(DOF_p,DOF_p,nnzEst_p); 
globalF = spalloc(DOF_p,1,nnzEst_p); 
nElements = size(connectivity,1);
c = 0;

for iElem = 1:nElements
        
nodesPerElement =size(connectivity(iElem, :),2);
Xe = [XYZ(connectivity(iElem, :), 1), XYZ(connectivity(iElem, :), 2)];
N1  = connectivity(iElem,1);                                                  % Node 1 for current element
N2  = connectivity(iElem,2);                                                  % Node 2 for current element
N3  = connectivity(iElem,3);                                                  % Node 3 for current element
N4  = connectivity(iElem,4);                                                  % Node 4 for current element
NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element    
HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
PEN = nnz(NN(:,1)); 
NEN = HEN+PEN;  

local_u = reshape([2 * connectivity(iElem, :) - 1; 2 * connectivity(iElem, :)], [], 1);
local_p = connectivity(iElem, :); 
iLoc_u = 9;                                                                                 
iLoc_p = 5;
localC = 0;
localH = 0;
localS = 0;
localF = 0;

if (NEN ~= 0)

[intersections, jacobian] = findCrackIntersections2(CRACK, Xe);

 if isempty(intersections)
        continue;  % Skip this element if no crack intersection is found
 end

 if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [ PSI(N1);  PSI(N2);  PSI(N3);  PSI(N4)];                 % Nodal crack level set values
 end

 [gp,gw] = gaussPointsOnCrack2(intersections, Xe, 2, elementType);
 c = c + 1;

     for gpIndex = 1:size(gp, 1)

     [N, dNLocal] = ShapeFunction(gp(gpIndex, :), elementType);
     J= (dNLocal * Xe);
     dNGlobal = J \ dNLocal;

     B_u = zeros(3, 2 * nodesPerElement);
     B_u(1, 1:2:end) = dNGlobal(1, :);  % ∂N/∂x for x-direction
     B_u(2, 2:2:end) = dNGlobal(2, :);  % ∂N/∂y for y-direction
     B_u(3, 1:2:end) = dNGlobal(2, :);  
     B_u(3, 2:2:end) = dNGlobal(1, :);  

     enr = find(NN(:, 1) ~= 0)';
     psi_modi = N(1,enr) * abs(PN(enr,1))- abs(N(1,enr) * PN(enr,1));
     dpsi_dx = dNGlobal(1,enr) * abs(PN(enr,1)) - dNGlobal(1,enr) * PN(enr,1) * sign (N(1,enr) * PN(enr,1));
     dpsi_dy = dNGlobal(2,enr) * abs(PN(enr,1)) - dNGlobal(2,enr) * PN(enr,1) * sign (N(1,enr) * PN(enr,1));
     dpsi = [dpsi_dx;dpsi_dy];  

     index_u = 1;
     index_p = 1;
     Benr_u = [];
     Benr_p = [];
     Nenr = [];
     for iN = 1:nodesPerElement
         if NN(iN,2) ~= 0
           
           Hi  = NN(iN,3);                                 
           H   = -Hi;          

           Ba = [dNGlobal(1,iN)*H   0;
                 0    dNGlobal(2,iN)*H
                 dNGlobal(2,iN)*H dNGlobal(1,iN)*H];

           Benr_u(:,index_u:(index_u+1)) = Ba;
           index_u = index_u+2;
                    
           if (gpIndex == length(gp))
               local_u(iLoc_u:(iLoc_u+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
               iLoc_u = iLoc_u+2;
           end
         end

        if NN(iN,1) ~= 0
       
           Benr_p(:,index_p) = dNGlobal(:,iN) * psi_modi + dpsi * N(iN);
           Nenr(:,index_p) = N(iN) * psi_modi;
           index_p = index_p+1;

          if (gpIndex == length(gp))
               local_p(iLoc_p) = NN(iN,1);
               iLoc_p = iLoc_p+1;
          end
        end
     end

     B1=[B_u, Benr_u];
     B2=[dNGlobal, Benr_p];
     B3=[N, Nenr];
     localC = localC + gw(gpIndex) * B3' * W(c) * (pt * B1) * det(jacobian);
     localH = localH + gw(gpIndex) * (t'*B2)' * W(c)^3/(12*mu) * (t'* B2) * det(jacobian);
     localS = localS + gw(gpIndex) * B3' * W(c) * com_f * B3 * det(jacobian);
     localF = localF + gw(gpIndex) * (t'*B2)' * W(c)^3/(12*mu) * rho_f * (t'* g) * det(jacobian);
    
     end
    globalQ(local_p,local_u) = globalQ(local_p,local_u) + localC;  
    globalH(local_p,local_p) = globalH(local_p,local_p) + localH; 
    globalS(local_p,local_p) = globalS(local_p,local_p) + localS; 
    globalF(local_p,1) = globalF(local_p,1) + localF; 
end
end
end

