function globalQ = CouplingMatrix(params, connectivity, XYZ, PSI, PHI, NODES, CRACK, elementType, omega, domain_length, nex)
   
% Extract material properties
alpha=params.alpha;
m = [1;1;0];
nCT = size(PHI,2); 
lXElem = domain_length/nex; 

% Initialize the FE coupling matrix
[DOF_u,DOF_p,~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * (DOF_u+DOF_p);
globalQ = spalloc(DOF_u,DOF_p,nnzEst); 

mm = size(CRACK,1);                                                          % Determine number of data points defining crack
if mm > 0
    if nCT == 1
        xCT = CRACK(mm,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(mm,2);                                                   % Y-coordinate of crack tip
    elseif nCT == 2
        xCT = [CRACK(mm,1) CRACK(1,1)];                                      % X-coordinates of crack tips
        yCT = [CRACK(mm,2) CRACK(1,2)];                                      % Y-coordinates of crack tips
    end
end

% Initialize the variables which will create the traditional sparse matrix
nIndexT  = 0;                                                 % Initialize index
nElements = size(connectivity,1);
[enrelem] = enrElem(connectivity, NODES);
uenrElem = nElements-length(enrelem);                         % Number of unenriched elements
allRowT  = ones(uenrElem*32,1);                               % Row indices
allColT  = ones(uenrElem*32,1);                               % Column indices
allValT  = zeros(uenrElem*32,1);      

for iElem = 1:nElements
        
Xe = [XYZ(connectivity(iElem, :), 1), XYZ(connectivity(iElem, :), 2)];
nodesPerElement =size(connectivity(iElem, :),2);
N1  = connectivity(iElem,1);                                                  % Node 1 for current element
N2  = connectivity(iElem,2);                                                  % Node 2 for current element
N3  = connectivity(iElem,3);                                                  % Node 3 for current element
N4  = connectivity(iElem,4);                                                  % Node 4 for current element
NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element    
CTN = nnz(NN(:,4));                                                     % Number of nodes with crack tip enrichment    
HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
PEN = nnz(NN(:,1));    
NEN = HEN+CTN+PEN;                                                    % Number of enriched nodes

localQ = 0;                                                                                  % Initialize stiffness for current element
local_u  = reshape([2 * connectivity(iElem, :) - 1; 2 * connectivity(iElem, :)], [], 1);             % Traditional index locations
local_p  = connectivity(iElem, :); 
iLoc_u=9;
iLoc_p=5;

if (NEN == 0)

    [gp, gw] = gaussPoints(elementType, 2);
    for gpIndex = 1:length(gp)

    [N, dNLocal] = ShapeFunction(gp(gpIndex,:), elementType);
    jacobian = dNLocal * Xe;
    dNGlobal = jacobian \ dNLocal;
         
    Bstd = zeros(3, 2 * nodesPerElement);
    Bstd(1, 1:2:end) = dNGlobal(1, :);  % ∂N/∂x for x-direction
    Bstd(2, 2:2:end) = dNGlobal(2, :);  % ∂N/∂y for y-direction
    Bstd(3, 1:2:end) = dNGlobal(2, :);  % ∂N/∂y for x-direction (shear)
    Bstd(3, 2:2:end) = dNGlobal(1, :);  % ∂N/∂x for y-direction (shear)
    localQ = localQ + gw(gpIndex) * Bstd' * m * alpha* N * det(jacobian);

    end

else

    if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [PSI(N1);  PSI(N2);  PSI(N3);  PSI(N4)];                 % Nodal crack level set values
    end
    crackTipInside = false;
     for ct = 1:nCT
     % Crack tip coordinates
     ctX = xCT(ct);
     ctY = yCT(ct);
     in = inpolygon(ctX, ctY, Xe(:,1), Xe(:,2));
     if in
        crackTipInside = true;
        break;
     end
     end

    if  CTN == 4 || crackTipInside
        [gp, gw, J] = subDomain(13, PN, Xe, 1, 0, 0,CRACK);
    elseif HEN == 4 || PEN ==4
        [gp, gw, J] = subDomain(7, PN, Xe, 0, 0, 0,CRACK);
    else
        [gp, gw, J] = subDomain(7, PN, Xe, 0, 0, 0,CRACK);
    end

for gpIndex = 1:length(gp)

    jacobian = [J(gpIndex,1) J(gpIndex,2);J(gpIndex,3) J(gpIndex,4)];                                                                                               
    [N, dNLocal] = ShapeFunction(gp(gpIndex,:), elementType);
    dNGlobal = (dNLocal * Xe) \ dNLocal;

    Xgp = N * Xe(:,1);                          % The global X for the current gauss point
    Ygp = N * Xe(:,2);                          % The global Y for the current gauss point

    Bstd = zeros(3, 2 * nodesPerElement);
    Bstd(1, 1:2:end) = dNGlobal(1, :);  % ∂N/∂x for x-direction
    Bstd(2, 2:2:end) = dNGlobal(2, :);  % ∂N/∂y for y-direction
    Bstd(3, 1:2:end) = dNGlobal(2, :);  % ∂N/∂y for x-direction (shear)
    Bstd(3, 2:2:end) = dNGlobal(1, :);  % ∂N/∂x for y-direction (shear)
      
    enr = find(NN(:, 1) ~= 0)';
    if ~isempty(enr)
        psi_modi = N(1,enr) * abs(PN(enr,1))- abs(N(1,enr) * PN(enr,1)); 
    end

    index_u = 1;
    index_p = 1;
    Benr = [];
    Nenr = [];
    for iN = 1:nodesPerElement
        if NN(iN,2) ~= 0

           psi1 = PSI(N1);                                         % Psi level set value at node 1
           psi2 = PSI(N2);                                         % Psi level set value at node 2
           psi3 = PSI(N3);                                         % Psi level set value at node 3
           psi4 = PSI(N4);                                         % Psi level set value at node 4
           psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Psi level set value at current gauss point
           Hgp = sign(psi);                                        % Heaviside value at current gauss point
           Hi  = NN(iN,3);                                         % Nodal Heaviside value
           H   = Hgp-Hi;                                           % Shifted Heaviside value

           Ba = [dNGlobal(1,iN)*H   0;
                 0    dNGlobal(2,iN)*H;
                 dNGlobal(2,iN)*H dNGlobal(1,iN)*H];

           Benr(:,index_u:(index_u+1)) = Ba;
           index_u = index_u+2;
                    
           if (gpIndex == length(gp))
               local_u(iLoc_u:(iLoc_u+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
               iLoc_u = iLoc_u+2;
           end
        elseif NN(iN,4) ~= 0

            if nCT == 1
               X     = Xgp-xCT;                                    % Horizontal distance from crack tip to gauss point
               Y     = Ygp-yCT;                                    % Vertical distance from crack tip to gauss point
               CCS   = [cos(omega) sin(omega);-sin(omega) cos(omega)];
               XYloc = CCS*[X Y]';                                 % Change to crack tip coordinates
               r     = sqrt(XYloc(1)^2+XYloc(2)^2);                % Radius from crack tip to current gauss point
               if r < 0.001*lXElem; r = 0.05*lXElem; end
               theta = atan2(XYloc(2),XYloc(1));                   % Angle from crack tip to current gauss point

            elseif nCT == 2
              X1  = (Xgp-xCT(1));
              Y1  = (Ygp-yCT(1));
              X2  = (Xgp-xCT(2));
              Y2  = (Ygp-yCT(2));
              CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
              XY1 = CCS*[X1 Y1]';
              CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];
              XY2 = CCS*[X2 Y2]';
              r1  = sqrt(XY1(1)^2+XY1(2)^2);                      % Radius from crack tip to current gauss point
              r2  = sqrt(XY2(1)^2+XY2(2)^2);
              if r1 > r2
                 r = r2; theta = atan2(XY2(2),XY2(1));
                 CCS = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];                            
              elseif r2 > r1
                 r = r1; theta = atan2(XY1(2),XY1(1));
                 CCS = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];
              end
              if r < 0.001*lXElem; r = 0.05*lXElem; end     
            end
            c = 1/(2*sqrt(r)); ct = CCS(1,1); st = CCS(1,2);         % Constants
           
            a1gp = sqrt(r)*sin(theta/2);                        % Node 1 crack tip enrichment value
            a2gp = sqrt(r)*cos(theta/2);                        % Node 2 crack tip enrichment value
            a3gp = sqrt(r)*sin(theta)*sin(theta/2);             % Node 3 crack tip enrichment value
            a4gp = sqrt(r)*sin(theta)*cos(theta/2);             % Node 4 crack tip enrichment value

            a1 = a1gp-NN(iN,5);                                 % Shifted alpha 1 enrichment value
            a2 = a2gp-NN(iN,7);                                 % Shifted alpha 2 enrichment value
            a3 = a3gp-NN(iN,9);                                 % Shifted alpha 3 enrichment value
            a4 = a4gp-NN(iN,11);                                % Shifted alpha 4 enrichment value

             % Derivative of crack tip enrichment functions with respect to X
             Px = c*[-sin(theta/2)*ct              + cos(theta/2)*-st;...
                     cos(theta/2)*ct               + sin(theta/2)*-st;...
                     -sin(3*theta/2)*sin(theta)*ct + (sin(theta/2)+sin(3*theta/2)*cos(theta))*-st;...
                     -cos(3*theta/2)*sin(theta)*ct + (cos(theta/2)+cos(3*theta/2)*cos(theta))*-st];

                        % Derivative of crack tip enrichment functions with respect to Y
             Py = c*[-sin(theta/2)*st              + cos(theta/2)*ct;...
                     cos(theta/2)*st               + sin(theta/2)*ct;...
                     -sin(3*theta/2)*sin(theta)*st + (sin(theta/2)+sin(3*theta/2)*cos(theta))*ct;...
                     -cos(3*theta/2)*sin(theta)*st + (cos(theta/2)+cos(3*theta/2)*cos(theta))*ct];

              B1 = [dNGlobal(1,iN)*a1+N(iN)*Px(1)          0;
                                       0            dNGlobal(2,iN)*a1+N(iN)*Py(1);
                              dNGlobal(2,iN)*a1+N(iN)*Py(1) dNGlobal(1,iN)*a1+N(iN)*Px(1)];

              B2 = [dNGlobal(1,iN)*a2+N(iN)*Px(2)          0;
                                       0            dNGlobal(2,iN)*a2+N(iN)*Py(2);
                              dNGlobal(2,iN)*a2+N(iN)*Py(2) dNGlobal(1,iN)*a2+N(iN)*Px(2)];

              B3 = [dNGlobal(1,iN)*a3+N(iN)*Px(3)          0;
                                       0            dNGlobal(2,iN)*a3+N(iN)*Py(3);
                              dNGlobal(2,iN)*a3+N(iN)*Py(3) dNGlobal(1,iN)*a3+N(iN)*Px(3)];

              B4 = [dNGlobal(1,iN)*a4+N(iN)*Px(4)          0;
                                       0            dNGlobal(2,iN)*a4+N(iN)*Py(4);
                              dNGlobal(2,iN)*a4+N(iN)*Py(4) dNGlobal(1,iN)*a4+N(iN)*Px(4)];
              
             Bb = [B1 B2 B3 B4];
             Benr(:,index_u:(index_u+7)) = Bb;
             index_u = index_u+8;

             if (gpIndex == length(gp))
                local_u(iLoc_u:(iLoc_u+7)) = [2*NN(iN,4)-1 2*NN(iN,4) 2*NN(iN,6)-1  2*NN(iN,6)...
                2*NN(iN,8)-1 2*NN(iN,8) 2*NN(iN,10)-1 2*NN(iN,10)];
                iLoc_u = iLoc_u+8;
             end
        end
       
        if NN(iN,1) ~= 0
        
        Nenr(:,index_p) = N(iN) * psi_modi;
        index_p = index_p+1;
                    
        if (gpIndex == length(gp))
            local_p(iLoc_p) = NN(iN,1);
            iLoc_p = iLoc_p+1;
        end
        end
         
    end
    Bu=[Bstd, Benr];
    Bp=[N, Nenr];
    localQ = localQ + gw(gpIndex) * Bu' * m * alpha * Bp * det(jacobian);
end
end

if size(localQ,1) == 8                                                  % Unenriched element
        for row = 1:8
            for col = 1:4
                nIndexT = nIndexT+1;
                allRowT(nIndexT) = local_u(row);
                allColT(nIndexT) = local_p(col);
                allValT(nIndexT) = localQ(row,col);
            end
        end
    else
        globalQ(local_u,local_p) = globalQ(local_u,local_p) + localQ;            % Assemble the global stiffness
end
end

globalQ = globalQ + sparse(allRowT,allColT,allValT,DOF_u,DOF_p);
end