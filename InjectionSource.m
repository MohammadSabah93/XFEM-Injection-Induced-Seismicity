function globalF = InjectionSource(boundaryConditions, NODES)
    
[~,globalDOF,~] = calcDOF(NODES);
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * globalDOF;
globalF = spalloc(globalDOF, 1,nnzEst);

concentratedSources = boundaryConditions.concentratedSources;
    for i = 1:size(concentratedSources, 1)
        nodeIndex = concentratedSources(i, 1);    % Global node index
        sourceValue = concentratedSources(i, 2);    % Source term value
        
        % Directly add the source value to the corresponding global DOF entry.
        globalF(nodeIndex) = globalF(nodeIndex) + sourceValue;
    end
    
end


