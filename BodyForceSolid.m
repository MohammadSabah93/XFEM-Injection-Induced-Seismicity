function globalF = BodyForceSolid(params, connectivity, XYZ, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType)
    
% Extract material properties
rho = params.rho_s;  % Density of the material
g = params.g;        % Gravitational acceleration
S = rho * g;

lXElem = domain_length/nex; 
nCT = size(PHI,2);                                                      % Number of crack tips 

% Initialize the FE stiffness matrix
[globalDOF,~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * globalDOF;
globalF = spalloc(globalDOF,1,nnzEst);    

m = size(CRACK,1);                                                          % Determine number of data points defining crack
if m > 0
    if nCT == 1
        xCT = CRACK(m,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(m,2);                                                   % Y-coordinate of crack tip
    elseif nCT == 2
        xCT = [CRACK(m,1) CRACK(1,1)];                                      % X-coordinates of crack tips
        yCT = [CRACK(m,2) CRACK(1,2)];                                      % Y-coordinates of crack tips
    end
end

% Initialize the variables which will create the traditional sparse matrix
nIndexT  = 0;                                                 % Initialize index
nElements = size(connectivity,1);
[enrelem] = enrElem(connectivity, NODES);
uenrElem = nElements-length(enrelem);                         % Number of unenriched elements
allRowT  = ones(uenrElem*64,1);                               % Row indices
allValT  = zeros(uenrElem*64,1);                              % Stiffness matrix values

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
NEN = HEN+CTN;  

localF = 0;                                                                                  % Initialize stiffness for current element
local  = reshape([2 * connectivity(iElem, :) - 1; 2 * connectivity(iElem, :)], [], 1);             % Traditional index locations
iLoc=9;

if (NEN == 0)

    [gp, gw] = gaussPoints(elementType, 2);
    for gpIndex = 1:length(gp)

    [N, dNLocal] = ShapeFunction(gp(gpIndex,:), elementType);
    jacobian = dNLocal * Xe;
    Nstd = zeros (2, 2*nodesPerElement);
    Nstd (1,1:2:end)= N;
    Nstd (1,2:2:end)= N;
    localF = localF + gw(gpIndex) * Nstd' * S * det(jacobian);

    end

elseif NEN > 0

     if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                 % Nodal crack level set values
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
    elseif HEN == 4
        [gp, gw, J] = subDomain(7, PN, Xe, 0, 0, 0,CRACK);
    else
        [gp, gw, J] = subDomain(7, PN, Xe, 0, 0, 0,CRACK);
    end
    

for gpIndex = 1:length(gp)
    
    jacobian = [J(gpIndex,1) J(gpIndex,2);J(gpIndex,3) J(gpIndex,4)];                       
                                                                               
    [N, ~] = ShapeFunction(gp(gpIndex,:), elementType);
    Nstd = zeros (2, 2*nodesPerElement);
    Nstd (1,1:2:end)= N;
    Nstd (1,2:2:end)= N;
    
    Xgp = N * Xe(:,1);                          % The global X for the current gauss point
    Ygp = N * Xe(:,2);                          % The global Y for the current gauss point
            
    index = 1;
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

           Ba = [N(1,iN)*H   0;
                 0    N(1,iN)*H];

           Nenr(:,index:(index+1)) = Ba;
           index = index+2;
                    
           if (gpIndex == length(gp))
               local(iLoc:(iLoc+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
               iLoc = iLoc+2;
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
              elseif r2 > r1
                 r = r1; theta = atan2(XY1(2),XY1(1)); 
              end
              if r < 0.001*lXElem; r = 0.05*lXElem; end     
            end
            
           
            a1gp = sqrt(r)*sin(theta/2);                        % Node 1 crack tip enrichment value
            a2gp = sqrt(r)*cos(theta/2);                        % Node 2 crack tip enrichment value
            a3gp = sqrt(r)*sin(theta)*sin(theta/2);             % Node 3 crack tip enrichment value
            a4gp = sqrt(r)*sin(theta)*cos(theta/2);             % Node 4 crack tip enrichment value

            a1 = a1gp-NN(iN,5);                                 % Shifted alpha 1 enrichment value
            a2 = a2gp-NN(iN,7);                                 % Shifted alpha 2 enrichment value
            a3 = a3gp-NN(iN,9);                                 % Shifted alpha 3 enrichment value
            a4 = a4gp-NN(iN,11);                                % Shifted alpha 4 enrichment value

             B1 = [N(iN)*a1    0;
                   0        N(iN)*a1];
                  

             B2 = [N(iN)*a2    0;
                    0       N(iN)*a2];
                    

             B3 = [N(iN)*a3    0;
                    0       N(iN)*a3];
                   
             B4 = [N(iN)*a4    0;
                   0        N(iN)*a4];
                   

             Bb = [B1 B2 B3 B4];
             Nenr(:,index:(index+7)) = Bb;
             index = index+8;

             if (gpIndex == length(gp))
                local(iLoc:(iLoc+7)) = [2*NN(iN,4)-1 2*NN(iN,4) 2*NN(iN,6)-1  2*NN(iN,6)...
                2*NN(iN,8)-1 2*NN(iN,8) 2*NN(iN,10)-1 2*NN(iN,10)];
                iLoc = iLoc+8;
             end
        end
   end
    B=[Nstd, Nenr];
    localF = localF + gw(gpIndex) * B' * S * det(jacobian);

end
end

if length(localF) == 8                                                  % Unenriched element
        for row = 1:8
            nIndexT = nIndexT+1;
            allRowT(nIndexT) = local(row);
            allValT(nIndexT) = localF(row);
        end
    else
        globalF(local,1) = globalF(local,1) + localF;            % Assemble the global force vector
end

end

globalF = globalF + sparse(allRowT,1,allValT,globalDOF,1);
end
    




   