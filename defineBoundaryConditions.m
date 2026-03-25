function boundaryConditions = defineBoundaryConditions(R,E)
    % Initialize boundary condition structure
    boundaryConditions = struct();

    tol = 1e-6;

    % Identify boundary nodes
    boundaryConditions.leftEdgeNodes = find(abs(R(:, 1) - min(R(:, 1))) < tol);
    boundaryConditions.rightEdgeNodes = find(abs(R(:, 1) - max(R(:, 1))) < tol);
    boundaryConditions.bottomEdgeNodes = find(abs(R(:, 2) - min(R(:, 2))) < tol);
    boundaryConditions.topEdgeNodes = find(abs(R(:, 2) - max(R(:, 2))) < tol);

    % Fully fixed nodes (zero displacement)
    fixedNodes = [];   % ensure column vector

    % Map node IDs to displacement DOFs (Ux = 2*node-1, Uy = 2*node)
    fullyFixedDOFs = [2*fixedNodes - 1; 2*fixedNodes];    % stack into one column
    fullyFixedDOFs = unique(fullyFixedDOFs, 'stable');    % remove duplicates (if any)

    % Zero displacement values
    fullyFixedDisp = zeros(numel(fullyFixedDOFs), 1);

    % Roller boundary conditions (e.g., only u_x fixed)
    rollerXNodes = [boundaryConditions.leftEdgeNodes; boundaryConditions.rightEdgeNodes];
    rollerXDOFs = 2*rollerXNodes - 1; % fix only u_x
    rollerXDisp = zeros(length(rollerXDOFs),1);

    % RollerYNodes 
    rollerYNodes = [boundaryConditions.bottomEdgeNodes;boundaryConditions.topEdgeNodes];
    rollerYDOFs = 2*rollerYNodes; % fix only u_y
    rollerYDisp = zeros(length(rollerYDOFs),1);

    % Nonzero prescribed displacement 
    EdgeNodes = [];
    EdgeDOFs = 2*EdgeNodes; 
    EdgeDisp = -0.1 * ones(length(EdgeDOFs),1); 

   % Traction boundaries
    boundaryConditions.tractionBoundaries = [];
    tractionConfig = struct(); 
    tractionConfig.leftEdge=[];
    tractionConfig.rightEdge=[];
    tractionConfig.topEdge=[];
    tractionConfig.bottomEdge=[];
    

    if  ~isempty(tractionConfig.leftEdge)
        Tx = tractionConfig.leftEdge(1);
        Ty = tractionConfig.leftEdge(2);
        leftTractionBoundaries = findBoundaryEdges(E,boundaryConditions.leftEdgeNodes);
        boundaryConditions.tractionBoundaries = [
            boundaryConditions.tractionBoundaries;
            leftTractionBoundaries, repmat([Tx, Ty], size(leftTractionBoundaries, 1), 1)
        ];
    end
    if ~isempty(tractionConfig.rightEdge)
        Tx = tractionConfig.rightEdge(1);
        Ty = tractionConfig.rightEdge(2);
        rightTractionBoundaries = findBoundaryEdges(E,boundaryConditions.rightEdgeNodes);
        boundaryConditions.tractionBoundaries = [
            boundaryConditions.tractionBoundaries;
            rightTractionBoundaries, repmat([Tx, Ty], size(rightTractionBoundaries,1), 1)
        ];
    end
    if ~isempty(tractionConfig.bottomEdge)
        Tx = tractionConfig.bottomEdge(1);
        Ty = tractionConfig.bottomEdge(2);
        bottomTractionBoundaries = findBoundaryEdges(E,boundaryConditions.bottomEdgeNodes);
        boundaryConditions.tractionBoundaries = [
            boundaryConditions.tractionBoundaries;
            bottomTractionBoundaries, repmat([Tx, Ty], size(bottomTractionBoundaries,1), 1)
        ];
    end
    if ~isempty(tractionConfig.topEdge)
        Tx = tractionConfig.topEdge(1);
        Ty = tractionConfig.topEdge(2);
        topTractionBoundaries = findBoundaryEdges(E,boundaryConditions.topEdgeNodes);
        boundaryConditions.tractionBoundaries = [
            boundaryConditions.tractionBoundaries;
            topTractionBoundaries, repmat([Tx, Ty], size(topTractionBoundaries,1), 1)
        ];
    end


% Flux boundaries
    boundaryConditions.fluxBoundaries = [];
    fluxConfig = struct(); 
    fluxConfig.leftEdge=[];
    fluxConfig.rightEdge=[];
    fluxConfig.topEdge=[];
    fluxConfig.bottomEdge=[];
   

    if  ~isempty(fluxConfig.leftEdge)
        q = fluxConfig.leftEdge;
        leftFluxBoundaries = findBoundaryEdges(E,boundaryConditions.leftEdgeNodes);
        boundaryConditions.fluxBoundaries = [
            boundaryConditions.fluxBoundaries;
            leftFluxBoundaries, repmat(q, size(leftFluxBoundaries, 1), 1)
        ];
    end
    if ~isempty(fluxConfig.rightEdge)
        q = fluxConfig.rightEdge;
        
        rightFluxBoundaries = findBoundaryEdges(E,boundaryConditions.rightEdgeNodes);
        boundaryConditions.fluxBoundaries = [
            boundaryConditions.fluxBoundaries;
            rightFluxBoundaries, repmat(q, size(rightFluxBoundaries,1), 1)
        ];
    end
    if ~isempty(fluxConfig.bottomEdge)
        q = fluxConfig.bottomEdge;
        bottomFluxBoundaries = findBoundaryEdges(E,boundaryConditions.bottomEdgeNodes);
        boundaryConditions.fluxBoundaries = [
            boundaryConditions.fluxBoundaries;
            bottomFluxBoundaries, repmat(q, size(bottomFluxBoundaries,1), 1)
        ];
    end
    if ~isempty(fluxConfig.topEdge)
        q = fluxConfig.topEdge;
        topFluxBoundaries = findBoundaryEdges(E,boundaryConditions.topEdgeNodes);
        boundaryConditions.fluxBoundaries = [
            boundaryConditions.fluxBoundaries;
            topFluxBoundaries, repmat(q, size(topFluxBoundaries,1), 1)
        ];
    end


    % Fixed pressure DOFs and values
    boundaryConditions.fixedPressureDOF = [];
    pressureConfig = struct();
    pressureConfig.leftEdge =[];  % Fixed pressure on left edge
    pressureConfig.rightEdge = 0; % Fixed pressure on right edge
    pressureConfig.topEdge = [];   % fixed pressure on top edge
    pressureConfig.bottomEdge = []; % Fixed pressure on bottom edge

    % Assign fixed pressure DOFs for left edge
    if ~isempty(pressureConfig.leftEdge)
        leftNodes = boundaryConditions.leftEdgeNodes;
        boundaryConditions.fixedPressureDOF = [
            boundaryConditions.fixedPressureDOF;
            leftNodes, repmat(pressureConfig.leftEdge, size(leftNodes, 1), 1)
        ];
    end

    % Assign fixed pressure DOFs for right edge
    if ~isempty(pressureConfig.rightEdge)
        rightNodes = boundaryConditions.rightEdgeNodes;
        boundaryConditions.fixedPressureDOF = [
            boundaryConditions.fixedPressureDOF;
            rightNodes, repmat(pressureConfig.rightEdge, size(rightNodes, 1), 1)
        ];
    end

    % Assign fixed pressure DOFs for bottom edge
    if ~isempty(pressureConfig.bottomEdge)
        bottomNodes = boundaryConditions.bottomEdgeNodes;
        boundaryConditions.fixedPressureDOF = [
            boundaryConditions.fixedPressureDOF;
            bottomNodes, repmat(pressureConfig.bottomEdge, size(bottomNodes, 1), 1)
        ];
    end

    % Assign fixed pressure DOFs for top edge
    if ~isempty(pressureConfig.topEdge)
        topNodes = boundaryConditions.topEdgeNodes;
        boundaryConditions.fixedPressureDOF = [
            boundaryConditions.fixedPressureDOF;
            topNodes, repmat(pressureConfig.topEdge, size(topNodes, 1), 1)
        ];
    end

    % Concentrated sources
    boundaryConditions.concentratedSources = [];
    boundaryConditions.concentratedSources = [5304,1e-4];
    % Displacement DOFs
    boundaryConditions.fixedDOF_disp = [
    fullyFixedDOFs , fullyFixedDisp;
    rollerXDOFs    , rollerXDisp;
    rollerYDOFs    , rollerYDisp;
    EdgeDOFs       , EdgeDisp];

    boundaryConditions.freeDOF_disp = setdiff(1:2 * size(R, 1), boundaryConditions.fixedDOF_disp);
end