function [connectivity_global, node_coordinates_global, nodes, elements] = MeshGeneration_2D(domain_length, domain_width, nex, ney, crackStart, crackEnd)

% User data
show_nodes = get_user_input('Do you want to display node numbers in the mesh?\nPress 1 for yes or 0 for no: ');
show_elements = get_user_input('Do you want to display element numbers in the mesh?\nPress 1 for yes or 0 for no: ');

%% Mesh generation
x = linspace(0, domain_length, nex + 1);
y = linspace(0, domain_width, ney + 1);

% Generate grid points
[X, Y] = meshgrid(x, y);
X = X'; Y = Y';
X = X(:); Y = Y(:); % Reshape to column vectors

% Node locations and total number of nodes/elements
node_coordinates = [(1:(nex+1)*(ney+1))', X, Y]; % Node number, x, y
nodes = size(node_coordinates, 1); % Total number of nodes
elements = nex * ney; % Total number of elements

% Plot nodes
figure(1)
scatter(X, Y, 10, 'filled'); % Plot nodes
axis equal; axis off;
hold on;

%% Element-node connectivity matrix
connectivity = zeros(elements, 4); % Preallocate for speed

el_id = 1;
for j = 1:ney
    for i = 1:nex
        n1 = (j-1)*(nex+1) + i;
        n2 = n1 + 1;
        n3 = n2 + nex + 1;
        n4 = n1 + nex + 1;
        connectivity(el_id, :) = [n1, n2, n3, n4]; % Without element ID
        % Plot element edges
        fill(node_coordinates([n1 n2 n3 n4 n1], 2), node_coordinates([n1 n2 n3 n4 n1], 3), 'w', 'EdgeColor', 'b');
        el_id = el_id + 1;
    end
end

%% Display node numbers if requested
if show_nodes
    text(X, Y, num2str((1:nodes)'), 'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', 'r');
end

%% Display element numbers if requested
if show_elements
    % Compute element centroids for labeling
    centroids_x = mean([X(connectivity(:,1)), X(connectivity(:,2)), X(connectivity(:,3)), X(connectivity(:,4))], 2);
    centroids_y = mean([Y(connectivity(:,1)), Y(connectivity(:,2)), Y(connectivity(:,3)), Y(connectivity(:,4))], 2);
    text(centroids_x, centroids_y, num2str((1:elements)'), 'FontSize', 8, 'Color', 'k', 'HorizontalAlignment', 'center');
end

hold on;

%% Plot the crack as a black line
plot([crackStart(1), crackEnd(1)], [crackStart(2), crackEnd(2)], 'k-', 'LineWidth', 2);

hold off;

%% Output global connectivity and node coordinate matrices
connectivity_global = connectivity; % Return without the element ID column
node_coordinates_global = node_coordinates(:, 2:3); % Return only the coordinates (x, y)

end

%% Helper function to get user input and validate
function result = get_user_input(prompt)
    result = input(prompt);
    while ~ismember(result, [0, 1])
        disp('Incorrect input! Please enter 1 for yes or 0 for no.');
        result = input(prompt);
    end
end
