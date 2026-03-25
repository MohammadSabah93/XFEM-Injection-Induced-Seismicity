function globalG = StabLagrange(params,interface_elements,interface_nodes)

M = params.Yu;
t=100;
nElements = size(interface_elements,1);
nInterfaceNodes = max(interface_elements(:));
avgBW = 20;  % tune this to your mesh
nnzEst = avgBW * nInterfaceNodes;
globalG = spalloc(nInterfaceNodes,nInterfaceNodes,nnzEst);    

for iElem = 1:nElements

localG = 0;                              % Initialize stiffness for current element

InterNodes = interface_elements(iElem, :);  % Node indices for this element
x1 = interface_nodes(InterNodes(1), :);     % [x, y] of node 1
x2 = interface_nodes(InterNodes(2), :);     % [x, y] of node 2
Le = norm(x2 - x1);                         % Element length in 2D
jacobian = Le / 2;                          % Jacobian for 1D is half the length
[gp, gw] = gaussPoints('1D', 2);

     for gpIndex = 1:size(gp, 1)

     [N, ~] = ShapeFunction(gp(gpIndex, :), '1D');
     
     Nstd = [N(1) - 0.5, N(2)-0.5];
   
     localG = localG + gw(gpIndex) * (t/(2*M))*(Nstd')* Nstd * det(jacobian);
    
     end
    globalG(InterNodes,InterNodes) = globalG(InterNodes,InterNodes) + localG; 
end
end

