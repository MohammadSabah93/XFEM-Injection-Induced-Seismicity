function globalS = StorageMatrix(params, connectivity, XYZ, PSI, PHI, NODES, CRACK, elementType)
   
% Extract material properties
S=params.Sc_m;
nCT = size(PHI,2);     
% Initialize the FE storage matrix
[~,DOF_p,~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * DOF_p;
globalS = spalloc(DOF_p,DOF_p,nnzEst);     

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
allRowT  = ones(uenrElem*16,1);                               % Row indices
allColT  = ones(uenrElem*16,1);                               % Column indices
allValT  = zeros(uenrElem*16,1);                              % Stiffness matrix values

for iElem = 1:nElements
        
nodesPerElement =size(connectivity(iElem, :),2);
Xe = [XYZ(connectivity(iElem, :), 1), XYZ(connectivity(iElem, :), 2)];
N1  = connectivity(iElem,1);                                                  % Node 1 for current element
N2  = connectivity(iElem,2);                                                  % Node 2 for current element
N3  = connectivity(iElem,3);                                                  % Node 3 for current element
N4  = connectivity(iElem,4);                                                  % Node 4 for current element
NN  = NODES([N1 N2 N3 N4]',:);                                          % Nodal data for current element    
CTN = nnz(NN(:,4));                                                     % Number of nodes with crack tip enrichment    
HEN = nnz(NN(:,1));                                                     % Number of nodes with pressure enrichment

localS = 0;                                  % Initialize stiffness for current element
local  = connectivity(iElem, :);             % Traditional index locations
iLoc=5;

if (HEN == 0)

    [gp, gw] = gaussPoints(elementType, 2);
    for gpIndex = 1:length(gp)

    [N, dNLocal] = ShapeFunction(gp(gpIndex,:), elementType);
    jacobian = dNLocal * Xe;
    localS = localS + gw(gpIndex) * N' * S * N * det(jacobian);

    end
else

     if numel(PSI) == 0, PN = [0 0 0 0]; else
        PN = [ PSI(N1);  PSI(N2) ; PSI(N3) ; PSI(N4)];                 % Nodal crack level set values
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

    enr = find(NN(:, 1) ~= 0)';
    if ~isempty(enr)
    psi_modi = N(1,enr) * abs(PN(enr,1))- abs(N(1,enr) * PN(enr,1)); 
    end

    index = 1;
    Nenr = [];
    for iN = 1:nodesPerElement
        if NN(iN,1) ~= 0
           
           Nenr(:,index) = N(iN) * psi_modi;
           index = index+1;
                    
           if (gpIndex == length(gp))
               local(iLoc) = NN(iN,1);
               iLoc = iLoc+1;
           end
        end
    end

    B=[N,Nenr];
    localS = localS + gw(gpIndex) * B' * S * B * det(jacobian);

end
end
if length(localS) == 4                                                  % Unenriched element
        for row = 1:4
            for col = 1:4
                nIndexT = nIndexT+1;
                allRowT(nIndexT) = local(row);
                allColT(nIndexT) = local(col);
                allValT(nIndexT) = localS(row,col);
            end
        end
    else
        globalS(local,local) = globalS(local,local) + localS;            % Assemble the global stiffness
end

end

globalS = globalS + sparse(allRowT,allColT,allValT,DOF_p,DOF_p);
end
