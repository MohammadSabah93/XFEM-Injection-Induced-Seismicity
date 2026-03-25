function postProcess(node_coordinates_global, U, P, connectivity_global, sigma_x_nodal, sigma_y_nodal, tau_xy_nodal, elementType)
    % Post-Processing and Visualization
    % Inputs:
    %   node_coordinates_global - Global node coordinates (nNodes x 2)
    %   U                       - Global displacement vector (2*nNodes x 1)
    %   P                       - Pressure values at nodes (nNodes x 1)
    %   connectivity_global     - Element connectivity matrix (nElements x nodesPerElement)
    %   sigma_x_nodal           - \sigma_x stresses at nodes (nNodes x 1)
    %   sigma_y_nodal           - \sigma_y stresses at nodes (nNodes x 1)
    %   tau_xy_nodal            - \tau_{xy} stresses at nodes (nNodes x 1)
    %   elementType             - Type of element ('TRI3' or 'QUAD4')

    % Extract nodal displacements
    Ux = U(1:2:end);
    Uy = U(2:2:end);

    % Compute deformed node positions
    deformedNodes = node_coordinates_global + [Ux, Uy];

    % Plot original and deformed mesh
    figure;
    hold on;
    if strcmp(elementType, 'TRI3')
        triplot(connectivity_global, node_coordinates_global(:, 1), node_coordinates_global(:, 2), 'k-');
        triplot(connectivity_global, deformedNodes(:, 1), deformedNodes(:, 2), 'r--');
    elseif strcmp(elementType, 'QUAD4')
        patch('Vertices', node_coordinates_global, 'Faces', connectivity_global, 'FaceColor', 'none', 'EdgeColor', 'k');
        patch('Vertices', deformedNodes, 'Faces', connectivity_global, 'FaceColor', 'none', 'EdgeColor', 'r', 'LineStyle', '--');
    end
    xlabel('x');
    ylabel('y');
    title('Original and Deformed Mesh');
    legend('Original Mesh', 'Deformed Mesh');
    hold off;

   plotFieldOnDeformedMesh(deformedNodes, connectivity_global, sigma_x_nodal, 'Sxx [Pa]');
   plotFieldOnDeformedMesh(deformedNodes, connectivity_global, sigma_y_nodal, 'Syy [Pa]');
   plotFieldOnDeformedMesh(deformedNodes, connectivity_global, tau_xy_nodal, 'Sxy [Pa]');
   plotFieldOnDeformedMesh(deformedNodes, connectivity_global, Ux, 'Displacement in x [m]');
   plotFieldOnDeformedMesh(deformedNodes, connectivity_global, Uy, 'Displacement in y [m]');
   plotFieldOnDeformedMesh(deformedNodes, connectivity_global, P, 'Pressure [Pa]');

   
end

function plotFieldOnDeformedMesh(deformedNodes, connectivity_global, nodalField, fieldLabel)
    % Plots a scalar field on the deformed mesh using nodal values
    % fieldLabel: e.g., '\sigma_x [MPa]', 'Displacement [m]'

    figure('Color', 'w'); 
    hold on;

    patch('Vertices', deformedNodes, ...
          'Faces', connectivity_global, ...
          'FaceVertexCData', nodalField, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    colormap(jet);
    hcb = colorbar;

    % Label the colorbar with field name and units
    ylabel(hcb, fieldLabel, ...
        'FontName', 'Calibri', ...
        'FontSize', 12, ...
        'Rotation', 270, ...
        'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold');

    clim([min(nodalField), max(nodalField)]);
    xlabel('x [m]');
    ylabel('y [m]');
    set(gca, 'XColor', 'none', 'YColor', 'none');  % hides axis lines
    axis equal;
    hold off;
end
