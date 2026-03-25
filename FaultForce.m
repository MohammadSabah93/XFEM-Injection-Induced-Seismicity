function globalF = FaultForce(connectivity, XYZ, NODES, CRACK, Sn, Ss, t, n, elementType, PHI, omega, domain_length, nex)



lXElem = domain_length/nex; 
nCT     = size(PHI,2); 

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


[DOF_u,~,~] = calcDOF(NODES);
nElements = size(connectivity,1);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * DOF_u;
globalF = spalloc (DOF_u,1,nnzEst);
counter = 1;

for iElem = 1:nElements
        
Xe = [XYZ(connectivity(iElem, :), 1), XYZ(connectivity(iElem, :), 2)];
nodesPerElement =size(connectivity(iElem, :),2);       
N1  = connectivity(iElem,1);                                                  % Node 1 for current element
N2  = connectivity(iElem,2);                                                  % Node 2 for current element
N3  = connectivity(iElem,3);                                                  % Node 3 for current element
N4  = connectivity(iElem,4);                                                  % Node 4 for current element
NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element    
HEN = nnz(NN(:,2));                                                     % Number of nodes with Heaviside enrichment
CTN = nnz(NN(:,4));
NEN = HEN+CTN;  

localF = 0;                                                            % Initialize stiffness for current element
local = 0;
iLoc = 1;


if (NEN ~= 0)

[intersections, jacobian] = findCrackIntersections2(CRACK, Xe);

   if isempty(intersections)
        continue;  % Skip this element if no crack intersection is found
   end
   
[gp,gw] = gaussPointsOnCrack2(intersections, Xe, 2, elementType);
P = Sn(counter,1)*n + Ss(counter,1)*t; % Normal and shear stress on fault

     for gpIndex = 1:size(gp, 1)

     [N, ~] = ShapeFunction(gp(gpIndex, :), elementType);
     Xgp = N * Xe(:,1);                          
     Ygp = N * Xe(:,2); 
     index = 1;
     Nstd = [];
     for iN = 1:nodesPerElement
        if NN(iN,2) ~= 0

          Ba = [N(1,iN)     0;
                0    N(1,iN)];
          
          Nstd(:,index:(index+1)) = Ba;
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

             if r1 > r2
             a1gp = sign(omega(1))*sqrt(r)*abs(sin(theta/2));                
             else
             a1gp = sign(omega(2))*sqrt(r)*abs(sin(theta/2));                  
             end

             B1 = [N(iN)*a1gp    0;
                   0        N(iN)*a1gp];

             Bb = B1;
             Nstd(:,index:(index+1)) = Bb;
             index = index+2;

             if (gpIndex == length(gp))
                local(iLoc:(iLoc+1)) = [2*NN(iN,4)-1 2*NN(iN,4)];
                iLoc = iLoc+2;
             end
        end
     end
    localF = localF + 2 * gw(gpIndex) * Nstd' * P * det(jacobian);
     end
    globalF(local,1) = globalF(local,1) + localF; 
    counter = counter +1;
end
end
end
