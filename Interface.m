function [avg_pn, avg_pt, avg_gn, avg_gt, LM_avg] = Interface(params, connectivity, XYZ, NODES, CRACK, U, n, t, interface_elements, LM, elementType, omega, domain_length, nex, PHI)


% Extract material properties
Ks=params.Ks;
Kn=params.Kn;

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


counter = 1;
nElements = size(connectivity,1);
nInterfaceElem = size(interface_elements,1);

pn=zeros(nInterfaceElem,1);    % Stiffness matrix
pt=zeros(nInterfaceElem,1); 
gn=zeros(nInterfaceElem,1);
gt=zeros(nInterfaceElem,1);
avg_pn = zeros(nInterfaceElem, 1);
avg_pt = zeros(nInterfaceElem, 1);
avg_gn = zeros(nInterfaceElem, 1);
avg_gt = zeros(nInterfaceElem, 1);


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
DOF=[];
if (NEN ~= 0)

[intersections, ~] = findCrackIntersections2(CRACK, Xe);

   if isempty(intersections)
        continue;  % Skip this element if no crack intersection is found
   end
   
[gp,~] = gaussPointsOnCrack2(intersections, Xe, 2, elementType);

     for gpIndex = 1:size(gp, 1)

      [N, ~] = ShapeFunction(gp(gpIndex, :), elementType);
      Xgp = N * Xe(:,1);                          
      Ygp = N * Xe(:,2);     
      index = 1;
      Benr = [];
      iLoc = 1;
      for iN = 1:nodesPerElement
        if NN(iN,2) ~= 0

          Ba = [N(1,iN)     0;
                0    N(1,iN)];
          
          Benr(:,index:(index+1)) = Ba;
          index = index+2;

          DOF(iLoc:(iLoc+1)) = [2*NN(iN,2)-1 2*NN(iN,2)];
          iLoc = iLoc+2;

       elseif NN(iN,4) ~= 0
            if nCT == 1
               XYZ     = Xgp-xCT;                                    % Horizontal distance from crack tip to gauss point
               Y     = Ygp-yCT;                                    % Vertical distance from crack tip to gauss point
               CCS   = [cos(omega) sin(omega);-sin(omega) cos(omega)];
               XYloc = CCS*[XYZ Y]';                                 % Change to crack tip coordinates
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
             Benr(:,index:(index+1)) = Bb;
             index = index+2;

             DOF(iLoc:(iLoc+1)) = [2*NN(iN,4)-1 2*NN(iN,4)];
             iLoc = iLoc+2;
             
        end
      end
    
     dis_jump = 2 * Benr * U(DOF);
     gn(counter,gpIndex)= n'*dis_jump;
     gt(counter,gpIndex)= t'*dis_jump;
     pn(counter,gpIndex) = Kn*gn(counter,gpIndex);      % Stiffness matrix
     pt(counter,gpIndex) = Ks*gt(counter,gpIndex); 
     end
     % Compute the average stress (normal and tangential) for this interface element
     avg_pn(counter) = mean(pn(counter, :));
     avg_pt(counter) = mean(pt(counter, :));
     avg_gn(counter) = mean(gn(counter, :));
     avg_gt(counter) = mean(gt(counter, :));
     counter = counter +1;
end
end


LM_avg = zeros(nInterfaceElem, 1);
for elem = 1:nInterfaceElem
    nodes = interface_elements(elem, :); % two nodes for 1D interface element
    LM_avg(elem) = 0.5 * (LM(nodes(1)) + LM(nodes(2))); % simple average
end

end







