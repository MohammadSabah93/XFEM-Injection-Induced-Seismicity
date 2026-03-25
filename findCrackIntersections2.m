function [intersections, jacobian] = findCrackIntersections2(CRACK, Xe)
% findCrackIntersections  Computes the intersection between a crack and an element edge.
%   [intersections, jacobian] = findCrackIntersections(CRACK, Xe)
%
%   Inputs:
%     CRACK - A 2x2 matrix where each row is a point [x y] of the crack.
%             Row 1 is the crack start and row 2 is the crack end.
%     Xe    - A nx2 matrix of nodal coordinates for a polygonal element.
%             n can be 3 (triangle), 4 (quadrilateral), or more.
%
%   Outputs:
%     intersections - A 2x2 matrix containing two unique intersection points
%                     [x, y] defining the portion of the crack inside the element.
%     jacobian      - A scalar computed as half the length of the crack segment 
%                     inside the element. If no valid intersection is found, [] is returned.
%
%   The function:
%     1. Checks each edge of the polygon (element) for an intersection with the crack.
%     2. Removes duplicate intersection points (using a tolerance to overcome round‐off).
%     3. If no boundary intersections are found but both crack endpoints fall
%        inside the element, the crack endpoints are used.
%     4. If only one intersection is found, the function checks whether one
%        of the endpoints is inside the element (and not coincident with the found intersection)
%        and adds it.
%     5. If more than two intersection points are computed (due to round‐off or complex
%        intersection scenarios), they are sorted along the crack line, and the extreme points
%        are taken.

    % Check for a valid polygonal element (at least 3 nodes)
    if size(Xe,1) < 3
        error('Element must have at least 3 nodes.');
    end
    % Check that CRACK is a 2x2 matrix.
    if ~isequal(size(CRACK), [2,2])
        error('CRACK must be a 2x2 matrix, where each row is a [x y] coordinate.');
    end

    % Tolerance for comparing floating point numbers.
    tol = 1e-10;
    
    % Initialize intersections array.
    intersections = [];
    
    % Extract crack endpoints for clarity.
    crackStart = CRACK(1,:);
    crackEnd   = CRACK(2,:);
    
    % Number of nodes in the element (polygon)
    num_nodes = size(Xe, 1);
    
    % Loop over each edge of the element.
    for i = 1:num_nodes
        % Determine the current edge as a segment from node i to node j.
        j = mod(i, num_nodes) + 1;  % cyclic index: after last node, go to first
        [is_int, x_int, y_int] = edgeIntersection(crackStart, crackEnd, Xe(i,:), Xe(j,:));
        if is_int
            intersections = [intersections; [x_int, y_int]];  %#ok<AGROW>
        end
    end
    
    % Remove duplicate intersection points using tolerance.
    % uniquetol (introduced in R2018b) removes rows that differ by less than tol.
    intersections = uniquetol(intersections, tol, 'ByRows', true);
    
    % If no intersections were found, check if the crack lies completely inside the element.
    if isempty(intersections)
        if isPointInsideElement(crackStart, Xe) && isPointInsideElement(crackEnd, Xe)
            intersections = [crackStart; crackEnd];
            jacobian = norm(crackEnd - crackStart) / 2;
            return;
        else
            jacobian = [];
            return;
        end
    end
    
    % If only one intersection is found, check if one of the crack endpoints is inside and not coincident.
    if size(intersections,1) == 1
        if isPointInsideElement(crackStart, Xe) && ~isApproxEqual(intersections, crackStart, tol)
            intersections = [intersections; crackStart];
        elseif isPointInsideElement(crackEnd, Xe) && ~isApproxEqual(intersections, crackEnd, tol)
            intersections = [intersections; crackEnd];
        else
            intersections = [];
            jacobian = [];
            return;
        end
    end
    
    % If more than two intersections are found, sort them by the parametric coordinate along the crack.
    if size(intersections,1) > 2
        % Compute the parametric coordinate t for each intersection on the crack:
        t_vals = zeros(size(intersections,1), 1);
        crackDir = crackEnd - crackStart;
        for k = 1:size(intersections,1)
            t_vals(k) = dot(intersections(k,:) - crackStart, crackDir) / (norm(crackDir)^2);
        end
        [~, sortIdx] = sort(t_vals);
        intersections = intersections(sortIdx,:);
        % Use the first and last points as the extreme intersections.
        intersections = [intersections(1,:); intersections(end,:)];
    end
    
    % If exactly two intersections are now available, compute the Jacobian.
    if size(intersections,1) == 2
        jacobian = norm(intersections(2,:) - intersections(1,:)) / 2;
    else
        % If something unexpected happens, return empty.
        intersections = [];
        jacobian = [];
    end
end

%% Helper Function: Edge Intersection
function [is_intersect, x_int, y_int] = edgeIntersection(crackStart, crackEnd, edgeStart, edgeEnd)
% edgeIntersection  Determines the intersection between two line segments.
%
%   [is_intersect, x_int, y_int] = edgeIntersection(crackStart, crackEnd, edgeStart, edgeEnd)
%
%   Inputs:
%     crackStart, crackEnd - 1x2 arrays for the crack endpoints.
%     edgeStart, edgeEnd   - 1x2 arrays for the edge endpoints.
%
%   Outputs:
%     is_intersect - A boolean flag that is true if the crack and edge intersect.
%     x_int, y_int - The coordinates of the intersection point (NaN if none).
    
    % Unpack coordinates for readability.
    x1_c = crackStart(1); y1_c = crackStart(2);
    x2_c = crackEnd(1);   y2_c = crackEnd(2);
    x1_e = edgeStart(1);  y1_e = edgeStart(2);
    x2_e = edgeEnd(1);    y2_e = edgeEnd(2);
    
    % Set up the 2x2 linear system
    A = [x2_c - x1_c, -(x2_e - x1_e);
         y2_c - y1_c, -(y2_e - y1_e)];
    b = [x1_e - x1_c; y1_e - y1_c];
    
    % Check for near singularity (parallel or collinear lines)
    if abs(det(A)) < 1e-10
        is_intersect = false;
        x_int = NaN; 
        y_int = NaN;
        return;
    end
    
    % Solve for the parameters t (along crack) and s (along edge)
    st = A \ b;
    t = st(1);  % parameter along the crack segment
    s = st(2);  % parameter along the edge segment
    
    % Determine if the intersection is within both segments (0 ≤ t,s ≤ 1)
    if t >= 0 && t <= 1 && s >= 0 && s <= 1
        is_intersect = true;
        x_int = x1_c + t * (x2_c - x1_c);
        y_int = y1_c + t * (y2_c - y1_c);
    else
        is_intersect = false;
        x_int = NaN; 
        y_int = NaN;
    end
end

%% Helper Function: Point-In-Polygon Test
function inside = isPointInsideElement(pt, Xe)
% isPointInsideElement  Determines if a point lies within a polygon.
%
%   inside = isPointInsideElement(pt, Xe)
%
%   Inputs:
%     pt - A 1x2 array [x y] representing the point.
%     Xe - An nx2 matrix of polygon vertices.
%
%   Output:
%     inside - A logical value that is true if the point is inside the polygon.
%
%   The polygon is explicitly closed by appending the first node to the end.
    
    polygonX = [Xe(:,1); Xe(1,1)];
    polygonY = [Xe(:,2); Xe(1,2)];
    inside = inpolygon(pt(1), pt(2), polygonX, polygonY);
end

%% Helper Function: Approximate Equality Check
function eq = isApproxEqual(pt1, pt2, tol)
% isApproxEqual  Checks if two 1x2 points are approximately equal within a tolerance.
%
%   eq = isApproxEqual(pt1, pt2, tol)
%
%   Inputs:
%     pt1, pt2 - 1x2 arrays containing [x y] coordinates.
%     tol      - A tolerance value (default is 1e-10 if not provided).
%
%   Output:
%     eq - True if the Euclidean distance between the points is less than tol.
    
    if nargin < 3
        tol = 1e-10;
    end
    eq = norm(pt1 - pt2) < tol;
end
