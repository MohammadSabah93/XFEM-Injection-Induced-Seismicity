function [omega, PSI, PHI, NODES, interface_elements, interface_nodes]=levelSet(connectivity, node_coordinates, nNodes, domain_length, nex, ney, CRACK)

NODES = zeros(nNodes,1);
lXElem = domain_length/nex;                               % Length of elements in the x-direction
nBand = (sqrt(2)+0.0005)*lXElem;                                % Radius of narrow band about crack

nPt = size(CRACK,1);                                      % Number of data points defining crack                                              

% Create the level set functions (PHI and PSI) defining the crack
if isempty(CRACK) == 0
    nCT = 2;                                                                % Default number of crack tips
    %% Check if crack is at the domain boundary (edge crack)
    if     (CRACK(1,1)   == 0) || (CRACK(1,1)   == nex*lXElem)           % Check for edge crack
        nCT = nCT-1;
    elseif (CRACK(nPt,1) == 0) || (CRACK(nPt,1) == nex*lXElem)
        nCT = nCT-1;
    elseif (CRACK(1,2)   == 0) || (CRACK(1,2)   == ney*lXElem)
        nCT = nCT-1;
    elseif (CRACK(nPt,2) == 0) || (CRACK(nPt,2) == ney*lXElem)
        nCT = nCT-1;
    end

    %% Compute Crack Orientation (Angle)**
    if nCT == 2
        % Vector from start to end
        vec = CRACK(2,:) - CRACK(1,:);
        omega(1) = atan2(vec(2), vec(1));
        omega(2) = atan2(-vec(2), -vec(1));  % reverse vector
    elseif nCT == 1
        % Compute orientation for a single crack tip
        disc  = -CRACK(nPt,:)+CRACK(nPt-1,:);                                % Horizontal and vertical distances for current crack segment
        omega = atan2(disc(2),disc(1));                                     % Crack angle with respect to horizontal
    end

    %% Initialize Level Set Functions**
    % Find the elements to search for new level set values
    PHI = sparse(nNodes,nCT);                                                % Initialize phi
    PSI = sparse(nNodes, 1); % PSI defines perpendicular distance from crack
  
    %%  Compute PHI and PSI for Each Crack Tip**
    for iCT = 1:nCT
        
        if iCT == 1
            dFSeg  = sqrt((CRACK(nPt-1,1)-CRACK(nPt,1))^2+(CRACK(nPt-1,2)-CRACK(nPt,2))^2);
            radius = dFSeg+nBand;                                               % Define radius of nodal search for defining level set
            xCT    = CRACK(nPt,1);                                              % X-coordinate of crack tip
            yCT    = CRACK(nPt,2);                                              % Y-coordinate of crack tip
            CCS    = [cos(omega(1)) sin(omega(1));-sin(omega(1)) cos(omega(1))];% Convert from global to crack tip coordinate system
        elseif iCT == 2
            dFSeg  = sqrt((CRACK(1,1)-CRACK(2,1))^2+(CRACK(1,2)-CRACK(2,2))^2);
            radius = dFSeg+nBand;                                               % Define radius of nodal search for defining level set
            xCT    = CRACK(1,1);                                                % X-coordinate of crack tip
            yCT    = CRACK(1,2);                                                % Y-coordinate of crack tip
            CCS    = [cos(omega(2)) sin(omega(2));-sin(omega(2)) cos(omega(2))];% Convert from global to crack tip coordinate system
        end
        
        dist = zeros(1,nNodes);                                              % Initialize distance vector
        for iN = 1:nNodes
            Xn       = node_coordinates(iN,1);                                           % X-coordinate for the current node
            Yn       = node_coordinates(iN,2);                                           % Y-coordinate for the current node
            X        = Xn-xCT;                                              % Horizontal distance from crack tip to current node
            Y        = Yn-yCT;                                              % Vertical distance from crack tip to current node
            XYloc    = CCS*[X Y]';                                          % Change to crack tip coordinates
            r        = sqrt(XYloc(1)^2+XYloc(2)^2);                         % Radius from crack tip to current gauss point
            dist(iN) = r;                                                   % Store radius value
        end
        
        temp   = dist-radius;                                               % Determine whether or not the node is outside the search radius
        domain = find(temp <= 0);                                           % Find nodes within search radius
        
        % Compute the PHI level set functions for the main crack tip(s)
        for iNode = 1:length(domain)
            cNode = domain(iNode);                                          % Current node within search radius
            x     = node_coordinates(cNode,1);                                           % X-coordinate for the current node
            y     = node_coordinates(cNode,2);                                           % Y-coordinate for the current node
            
            % Define phi for first crack tip
            disc = CRACK(nPt,:)-CRACK(nPt-1,:);                             % Horizontal and vertical distances for current crack segment
            t    = 1/norm(disc)*disc;                                       % Tangent to current crack segment
            phi  = ([x y]-CRACK(nPt,:))*t';
            if phi == 0, phi = 1e-6; end
            
            if nCT == 2
                disc = CRACK(1,:)-CRACK(2,:);                               % Horizontal and vertical distances for current crack segment
                t    = 1/norm(disc)*disc;                                   % Tangent to current crack segment
                phi(2)  = ([x y]-CRACK(1,:))*t';
                if phi(2) == 0, phi(2) = 1e-6; end
            end
            
            % Define phi and psi at nodes within narrow band
            if iCT == 1
                if abs(phi(1)) < nBand                                      % Check if phi is within narrow band
                    
                    % Define psi
                    dist = zeros(nPt-1,1); sine = dist;
                    for iSeg = 1:(nPt-1)
                        x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                        x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);

                        if x1 == x2
                        dist(iSeg) = abs(x - x1);  % Distance is horizontal distance
                        sine(iSeg) = sign(x - x1); % Left side (-), Right side (+)
                        else
                        m = (y2-y1)/(x2-x1);
                        if isinf(m) == 1, m = 1e6; end
                        b = y1-m*x1;
                        
                        xo = (m*y+x-m*b)/(m^2+1);
                        yo = (m^2*y+m*x+b)/(m^2+1);
                        
                        if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                        if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                        dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                        sine(iSeg) = -sign(x2-x1)*sign(y-yo);
                        end
                    end
                    
                    dMin = min(abs(dist));
                    ind  = find(abs(dist) == dMin);
                    if length(ind) == 2, ind(2) = []; end
                    psi = dist(ind)*sine(ind);
                    if psi == 0, psi = 1e-6; end
                    
                    % Check if psi is within defined narrow band
                    if abs(psi) < nBand
                        PHI(cNode,1) = phi(1);
                        PSI(cNode,1) = psi;
                    end
                elseif phi(1) < 0                                           % Check if phi is negative
                    
                    % Define psi
                    dist = zeros(nPt-1,1); sine = dist;
                    for iSeg = 1:(nPt-1)
                        x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                        x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);
                        if x1 == x2
                        dist(iSeg) = abs(x - x1);  % Distance is horizontal distance
                        sine(iSeg) = sign(x - x1); % Left side (-), Right side (+)
                        else
                        m = (y2-y1)/(x2-x1);
                        if isinf(m) == 1, m = 1e6; end
                        b = y1-m*x1;
                        
                        xo = (m*y+x-m*b)/(m^2+1);
                        yo = (m^2*y+m*x+b)/(m^2+1);
                        
                        if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                        if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                        dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                        sine(iSeg) = -sign(x2-x1)*sign(y-yo);
                        
                        end
                    end
                    
                    dMin = min(abs(dist));
                    ind  = find(abs(dist) == dMin);
                    if length(ind) == 2, ind(2) = []; end
                    psi = dist(ind)*sine(ind);
                    if psi == 0, psi = 1e-6; end
                    
                    % Check if psi is within narrow band
                    if abs(psi) < nBand
                        PSI(cNode,1) = psi;
                    end
                end
            end
            
            % Define phi for second crack tip
            if iCT == 2
                if abs(phi(2)) < nBand                                      % Check if phi is within narrow band
                    
                    % Define psi
                    dist = zeros(nPt-1,1); sine = dist;
                    for iSeg = 1:(nPt-1)
                        x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                        x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);

                        if x1 == x2
                        dist(iSeg) = abs(x - x1);  % Distance is horizontal distance
                        sine(iSeg) = sign(x - x1); % Left side (-), Right side (+)
                        else

                        m = (y2-y1)/(x2-x1);
                        if isinf(m) == 1, m = 1e6; end
                        b = y1-m*x1;
                        
                        xo = (m*y+x-m*b)/(m^2+1);
                        yo = (m^2*y+m*x+b)/(m^2+1);
                        
                        if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                        if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                        
                        dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                        sine(iSeg) = -sign(x2-x1)*sign(y-yo);
                         
                        end
                    end
                    
                    dMin = min(abs(dist));
                    ind  = find(abs(dist) == dMin);
                    if length(ind) == 2, ind(2) = []; end
                    psi = dist(ind)*sine(ind);
                    if psi == 0, psi = 1e-6; end
                    
                    % Check if psi is within defined narrow band
                    if abs(psi) < nBand
                        PHI(cNode,1) = phi(2);
                        PSI(cNode,1) = psi;
                    end
                elseif phi(2) < 0                                           % Check if phi is negative
                    if phi(1) < 0
                        % Define psi
                        dist = zeros(nPt-1,1); sine = dist;
                        for iSeg = 1:(nPt-1)
                            x1 = CRACK(iSeg,1);            y1 = CRACK(iSeg,2);
                            x2 = CRACK(iSeg+1,1);          y2 = CRACK(iSeg+1,2);

                            if x1 == x2
                            dist(iSeg) = abs(x - x1);  % Distance is horizontal distance
                            sine(iSeg) = sign(x - x1); % Left side (-), Right side (+)
                            else

                            m = (y2-y1)/(x2-x1);
                            if isinf(m) == 1, m = 1e6; end
                            b = y1-m*x1;
                            
                            xo = (m*y+x-m*b)/(m^2+1);
                            yo = (m^2*y+m*x+b)/(m^2+1);
                            
                            if iSeg ~= 1, if xo < x1, xo = x1; yo = y1; end, end
                            if iSeg ~= (nPt-1), if xo > x2, xo = x2; yo = y2; end, end
                            
                            dist(iSeg) = sqrt((xo-x)^2+(yo-y)^2);
                            sine(iSeg) = -sign(x2-x1)*sign(y-yo);
                            end
                        end
                        
                        dMin = min(abs(dist));
                        ind  = find(abs(dist) == dMin);
                        if length(ind) == 2, ind(2) = []; end
                        psi = dist(ind)*sine(ind);
                        if psi == 0, psi = 1e-6; end
                        
                        % Check if psi is within narrow band
                        if abs(psi) < nBand
                            PSI(cNode,1) = psi;
                        end
                    end
                end
            end
        end
    end
    
    % Define the crack tip enriched nodes
    ctNodes  = [];
    I        = find(PHI ~= 0);                                              % Nodes with defined phi
    [~,~,c1] = intersect(I,connectivity(:,1)');                             % Find elements with defined phi
    [~,~,c2] = intersect(I,connectivity(:,2)');                             % Find elements with defined phi
    [~,~,c3] = intersect(I,connectivity(:,3)');                             % Find elements with defined phi
    [~,~,c4] = intersect(I,connectivity(:,4)');                             % Find elements with defined phi
    ctElem   = unique([c1(:); c2(:); c3(:); c4(:)]);                        % Candidate elements for crack tip enrichment
    
    for iElem = 1:length(ctElem)
        cElem = ctElem(iElem);
        phiE  = PHI(connectivity(cElem,1:4));
        if nnz(phiE) == 4
            psiE = PSI(connectivity(cElem,1:4));
            if max(psiE)*min(psiE) < 1e-6
                if max(phiE)*min(phiE) < 1e-6
                    for iN = 1:4
                        ctNodes = [ctNodes connectivity(cElem,iN)];
                    end
                end
            end
        end
    end
    % ctNodes=[];
    NODES(ctNodes,4) = NaN;
   
    % Define Heaviside enriched nodes
    hNodes   = [];
    I        = find(PSI ~= 0);                                              % Nodes with defined psi
    [~,~,c1] = intersect(I,connectivity(:,1)');                             % Find elements with defined psi
    [~,~,c2] = intersect(I,connectivity(:,2)');                             % Find elements with defined psi
    [~,~,c3] = intersect(I,connectivity(:,3)');                             % Find elements with defined psi
    [~,~,c4] = intersect(I,connectivity(:,4)');                             % Find elements with defined psi
    hElem    = unique([c1(:); c2(:); c3(:); c4(:)]);                        % Candidate elements for Heaviside enrichment

    for iElem = 1:length(hElem)
        cElem = hElem(iElem);
        psiE  = PSI(connectivity(cElem,1:4));
        if nnz(psiE) == 4
            if (max(psiE) == 1e-6) || (min(psiE) == 1e-6)
                for iN = 1:4
                    gN = connectivity(cElem,iN);
                    if NODES(gN,4) == 0
                        if psiE(iN) == 1e-6
                            if PHI(gN) <= 0
                                if NODES(gN,2) == 0
                                    hNodes = [hNodes gN];
                                end
                            end
                        end
                    end
                end
            elseif max(psiE)*min(psiE) < 0
                for iN = 1:4
                    gN = connectivity(cElem,iN);
                    if PHI(gN) <= 0
                        if NODES(gN,4) == 0
                            if NODES(gN,2) == 0
                                hNodes = [hNodes gN];
                            end
                        end
                    end
                end
            end
        end
    end

 % for iElem = 1:size(connectivity,1)
 % 
 %        Xe = [node_coordinates(connectivity(iElem, :), 1), node_coordinates(connectivity(iElem, :), 2)];
 % 
 %        [intersections, ~] = findCrackIntersections2(CRACK, Xe);
 %        if ~isempty(intersections)
 %        for iN = 1:4    
 %        gN = connectivity(iElem,iN);
 %        hNodes = [hNodes gN];
 %        end
 %        end
 % 
 % end
 % Define pressure modified enriched nodes
pNodes   = [];
I        = find(PSI ~= 0);                                              % Nodes with defined psi
[~,~,c1] = intersect(I,connectivity(:,1)');                             % Find elements with defined psi
[~,~,c2] = intersect(I,connectivity(:,2)');                             % Find elements with defined psi
[~,~,c3] = intersect(I,connectivity(:,3)');                             % Find elements with defined psi
[~,~,c4] = intersect(I,connectivity(:,4)');                             % Find elements with defined psi
pElem    = unique([c1(:); c2(:); c3(:); c4(:)]);                        % Candidate elements for Heaviside enrichment

    for iElem = 1:length(pElem)
        cElem = pElem(iElem);
        psiE  = PSI(connectivity(cElem,1:4));
        if nnz(psiE) == 4
            if (max(psiE) == 1e-6) || (min(psiE) == 1e-6)
                for iN = 1:4
                    gN = connectivity(cElem,iN);
                        if psiE(iN) == 1e-6
                            if PHI(gN) <= 0
                               pNodes = [pNodes gN];
                            end
                        end
                end
            elseif max(psiE)*min(psiE) < 0
                for iN = 1:4
                    gN = connectivity(cElem,iN);
                    if PHI(gN) <= 0
                       if NODES(gN,1) == 0
                              pNodes = [pNodes gN];
                       end
                    end
                end
            end
        end
    end
% Preliminaries to numbering enriched nodes
nNodes   = nNodes+1;
nNodes2=nNodes;
pNodes = unique(pNodes);
ctNodes = unique(ctNodes);
hNodes  = unique(hNodes);
    

% % find std nodes in blending elements around crack tip    
blendNode = [];    
% for iElem = 1:size(connectivity,1)
%     elementNodes = connectivity(iElem, :);           % Nodes of the current element
%     % Find which nodes are already tip enriched (from ctNodes)
%     enriched = intersect(elementNodes, ctNodes);
% 
%     % Check if this element is a blending element:
%     % It must have at least one enriched node, but not all nodes enriched.
%     if ~isempty(enriched) && (length(enriched) < length(elementNodes))
%         % This element is a blending element
%         % Loop through each node in the element
%         for j = 1:length(elementNodes)
%             nodeId = elementNodes(j);
%             % If the node is not tip enriched (i.e. not in ctNodes)
%             if ~ismember(nodeId, ctNodes)
%                 % And if this node has not yet been assigned a blending number.
%                 % Here we check column 11 of NODES (adjust if you use more columns)
%                 blendNode=[blendNode nodeId];
%             end
%         end
%     end
% end

blendNode = unique(blendNode);   
% blenNode_heavi=[299,300,1027,1028];
% allEnriched_hevi = union(hNodes, blenNode_heavi);  % merge and remove duplicates

% Number the pressure modified nodes
for i = 1:length(pNodes)
    NODES(pNodes(i),1) = nNodes2;
    nNodes2 = nNodes2+1;
end
% Number the Heaviside nodes
for i = 1:length(hNodes)
    NODES(hNodes(i),2) = nNodes;
    nNodes = nNodes+1;
end
    
allEnriched = union(ctNodes, blendNode);  % merge and remove duplicates
% ---------- Assign DOFs ----------
for i = 1:length(allEnriched)
    nodeId = allEnriched(i);
    NODES(nodeId,4)  = nNodes; nNodes = nNodes + 1;
    NODES(nodeId,6)  = nNodes; nNodes = nNodes + 1;
    NODES(nodeId,8)  = nNodes; nNodes = nNodes + 1;
    NODES(nodeId,10) = nNodes; nNodes = nNodes + 1;
end  
    
else
    omega  = [];
end

for iNode = 1:size(NODES,1)
    if NODES(iNode,2) ~= 0
       NODES(iNode,3) = sign(PSI(iNode));
    end
end

% This function assigns nodes enriched with the crack tip enrichment
% function values associated with the near tip asymptotic displacement
% field.
iSeg = size(CRACK,1);                                                       % Number of crack segments
nCT  = size(PHI,2);                                                         % Number of crack tips

% Define coordinates of crack tip(s)
if nCT == 1
    xCT = CRACK(iSeg,1);                                                    % X-coordinate of crack tip
    yCT = CRACK(iSeg,2);                                                    % Y-coordinate of crack tip
elseif nCT == 2
    xCT = [CRACK(iSeg,1) CRACK(1,1)];                                       % X-coordinate of crack tips
    yCT = [CRACK(iSeg,2) CRACK(1,2)];                                       % Y-coordinate of crack tips
end

for iNode = 1:size(NODES,1)
    if NODES(iNode,4) ~= 0
    XN = node_coordinates(iNode,1);                                                  % Nodal x-coordinate
    YN = node_coordinates(iNode,2);                                                  % Nodal y-coordinate
    X1  = (XN-xCT(1));
    Y1  = (YN-yCT(1));
    X2  = (XN-xCT(2));
    Y2  = (YN-yCT(2));
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
        NODES(iNode,5)  = sqrt(r)*sin(theta/2);                             % Alpha 1 crack tip enrichment value
        NODES(iNode,7)  = sqrt(r)*cos(theta/2);                             % Alpha 2 crack tip enrichment value
        NODES(iNode,9)  = sqrt(r)*sin(theta)*sin(theta/2);                  % Alpha 3 crack tip enrichment value
        NODES(iNode,11) = sqrt(r)*sin(theta)*cos(theta/2);                  % Alpha 4 crack tip enrichment value  
    end
end

% Determine number of interface elements and nodes

interface_elements = [];
interface_nodes = [];
nodeID = 0;               
elemID = 0;  
nElements = nex*ney;
for iElem = 1:nElements
 
    elementNodes = connectivity(iElem, :)';
    Xe = node_coordinates(elementNodes, :);  
    HEN = nnz(NODES(elementNodes, 2));
    CTN = nnz(NODES(elementNodes, 4));
    NEN = HEN+CTN; 
    
    if NEN ~= 0
        
        [intersections, ~] = findCrackIntersections2(CRACK, Xe);  
        
        if ~isempty(intersections)  

            nodeID1 = nodeID + 1;
            nodeID2 = nodeID + 2;
            nodeID = nodeID + 1;           
            elemID = elemID + 1;
            interface_elements(elemID, :) = [nodeID1, nodeID2];
            interface_nodes(nodeID1,:) = intersections (1,:);
            interface_nodes(nodeID2,:) = intersections (2,:);
        end
    end
end
end

