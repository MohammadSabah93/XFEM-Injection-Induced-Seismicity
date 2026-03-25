function [gp, gw, J] = subDomain(npt, psi, xyz, ctflag, jiflag, CCS, CRACK)
% This function subdivides enriched elements and determines the Gauss 
% points and weights to be used in the integration during the assembly 
% of the stiffness matrix.
%
% Inputs:
%   - npt: Number of Gauss points per sub-triangle
%   - psi: Level set function values at element nodes
%   - xyz: Coordinates of the element nodes (4x2 for quadrilateral)
%   - ctflag: Crack tip flag (1 if crack tip enrichment is needed, 0 otherwise)
%   - jiflag: Jacobian flag (1 if transformation to crack tip coordinate system is required)
%   - CCS: Crack tip coordinate system transformation matrix
%
% Outputs:
%   - gp: Gauss points in parametric space
%   - gw: Gauss weights for integration
%   - J: Jacobian matrix for each Gauss point

% Initialize variables
ntip = []; 
xCT = []; 
yCT = [];

%% Identify Crack Tip in the Element (if needed)
if ctflag == 1
    m = size(CRACK,1);  % Number of points defining the crack
    xCTip = [CRACK(m,1), CRACK(1,1)];  % X-coordinates of crack tip
    yCTip = [CRACK(m,2), CRACK(1,2)];  % Y-coordinates of crack tip

    % Find bounding box of the element
    xMn = min(xyz(:,1)); xMx = max(xyz(:,1)); xMd = (xMn + xMx) / 2;
    yMn = min(xyz(:,2)); yMx = max(xyz(:,2)); yMd = (yMn + yMx) / 2;
    
    % Locate crack tip within the element's bounding box
    for i = 1:length(xCTip)
        if (xCTip(i) >= xMn) && (xCTip(i) <= xMx) && (yCTip(i) >= yMn) && (yCTip(i) <= yMx)
            xCT = xCTip(i); 
            yCT = yCTip(i);
            break;
        end
    end

    % If crack tip is not found, select the closest point to the element center
    if isempty(xCT)
        d = sqrt((xCTip - xMd).^2 + (yCTip - yMd).^2);
        [~, ind] = min(d);
        xCT = xCTip(ind); 
        yCT = yCTip(ind);
    end

    % Add crack tip as an additional node
    xyz = [xyz; xCT, yCT];

    % Map crack tip to parametric space
    le = xyz(2,1) - xyz(1,1);
    xm = mean(xyz(1:2,1));
    ym = mean(xyz(2:3,2));
    xi = 2 * (xCT - xm) / le;
    eta = 2 * (yCT - ym) / le;
    ntip = [xi, eta];
end

%% Construct Subdomain Nodes
corner = [1, 2, 3, 4, 1];  % Element corner connectivity
node = [-1, -1; 1, -1; 1, 1; -1, 1];  % Parametric space coordinates

% Include crack tip node if applicable
node = [node; ntip];

% Find intersection points of the level set function (if any)
if ~isempty(psi)
    for i = 1:4
        n1 = corner(i);
        n2 = corner(i+1);
        if psi(n1) * psi(n2) < 0  % Sign change indicates intersection
            r = psi(n1) / (psi(n1) - psi(n2));
            pnt = (1 - r) * node(n1, :) + r * node(n2, :);
            xi = pnt(1); eta = pnt(2);
            N = 1/4 * [(1-xi)*(1-eta); (1+xi)*(1-eta); (1+xi)*(1+eta); (1-xi)*(1+eta)];
            xpnt = dot(N, xyz(1:4,1)');
            ypnt = dot(N, xyz(1:4,2)');
            xyz = [xyz; xpnt, ypnt];  % Store intersection point
            node = [node; pnt];  % Add to node set
        end
    end
end

%% Triangulate the Subdomain
tol = 1e-6;  % Set an appropriate tolerance
node_rounded = round(node/tol)*tol;
[node_unique, ~, ~] = unique(node_rounded, 'rows', 'stable');
[xyz_unique, ~, ~] = unique(xyz, 'rows', 'stable');
tri = delaunay(node_unique(:,1), node_unique(:,2));      % Perform Delaunay triangulation

%% Compute Gauss Points and Weights for Each Sub-Triangle
[q, w] = gaussPoints('TRI3', npt);  % Get Gauss points & weights for triangular elements

%% Transform to Crack Tip Coordinate System (if needed)
if jiflag == 1
    m = size(CRACK,1);
    xCT = CRACK(m,1); 
    yCT = CRACK(m,2);

    % Translate to crack tip as origin
    xyz(:,1) = xyz(:,1) - xCT;  
    xyz(:,2) = xyz(:,2) - yCT;

    % Rotate to crack tip coordinate system
    xyz = (CCS * xyz')';
end

%% Compute Gauss Points and Jacobian for Each Subdomain
pt = 1;
gp = []; gw = []; J = [];

for e = 1:size(tri,1)
    coord = node_unique(tri(e,:),:);  % Nodes of the current triangle
    xyzl = xyz_unique(tri(e,:),:);  % Physical coordinates of triangle nodes
    
    for i = 1:length(w)
        xi = q(i,1); eta = q(i,2);
        N = [1 - xi - eta; xi; eta];  % Shape functions
        gp(pt,:) = N' * coord;  % Compute Gauss point in parametric space
        gw(pt,1) = w(i) / 2;    % Adjusted weight for area integration

        % Compute Jacobian
        J(pt,:) = [-xyzl(1,1) + xyzl(2,1), -xyzl(1,2) + xyzl(2,2), ...
                   -xyzl(1,1) + xyzl(3,1), -xyzl(1,2) + xyzl(3,2)];

        pt = pt + 1;
    end
end
