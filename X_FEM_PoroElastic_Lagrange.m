clc;
clear;
close all;

% 2D Mesh Generation 

% Generate 2D mesh
domain_length = 100;
domain_width = 100;
nex = 100;
ney = 100;

crackStart = [20 ,20.5];
crackEnd = [80, 80.8];

CRACK=[crackStart;crackEnd];

[connectivity, node_coordinates, nNodes, nElements] = MeshGeneration_2D(domain_length, domain_width, nex, ney, crackStart, crackEnd);

[omega, PSI, PHI, NODES, interface_elements, interface_nodes]=levelSet(connectivity, node_coordinates, nNodes, domain_length, nex, ney, CRACK);

nodesPerElement = size(connectivity, 2);
if nodesPerElement == 4
   elementType = 'QUAD4';
else
   elementType = 'TRI3';
end

%% Calculate unit normal and tangential vectors on fault

  % Extract coordinates of the crack
    x1 = crackStart(1);
    y1 = crackStart(2);
    x2 = crackEnd(1);
    y2 = crackEnd(2);

    % Compute the direction vector and length of the crack
    crackVector = [x2 - x1; y2 - y1];
    crackLengthSquared = dot(crackVector, crackVector); % Avoid recomputing square

    % Normalize the normal vector to the crack
    n = [crackVector(2); -crackVector(1)] / sqrt(crackLengthSquared);

    % Normalize the tangential vector to the crack
    t = crackVector / sqrt(crackLengthSquared);

    crackLength = sqrt(crackLengthSquared);

%% Model Parameters

params = defineModelParameters();

%% Boundary Conditions

boundaryConditions = defineBoundaryConditions(node_coordinates,connectivity);


%% Pre-allocate solution vectors
[DOF_u,DOF_p,DISP] = calcDOF(NODES);
nInterfaceNodes = max(interface_elements(:));
nInterfaceElem = size(interface_elements,1);

params.a = params.a*ones(nInterfaceElem,1);
params.b = params.b*ones(nInterfaceElem,1);
params.Dc = params.Dc*ones(nInterfaceElem,1);

a = zeros(2*nnz(NODES(:,2)),1);
Ux = zeros(nNodes, 1); 
Uy = zeros(nNodes, 1); 
U_ddotx = zeros(nNodes, 1); 
U_ddoty = zeros(nNodes, 1);
U_ddot_total = zeros(nNodes, 1);
U_dotx = zeros(nNodes, 1); 
U_doty = zeros(nNodes, 1);
U_dot_total = zeros(nNodes, 1);
P = zeros(DOF_p, 1);
Pp = zeros(nNodes, 1);
Dn = zeros(nInterfaceElem,1); 
Ds = zeros(nInterfaceElem,1); 
Sn = zeros(nInterfaceElem,1); 
Ss = zeros(nInterfaceElem,1);
Pn = zeros(nInterfaceElem,1); 
W = zeros(nInterfaceElem,1);
LM = zeros(nInterfaceNodes,1);
slip = zeros(nInterfaceElem,1);
uf = zeros(nInterfaceElem,1);
theta = zeros(nInterfaceElem,1);
Vs = zeros(nInterfaceElem,1);
stickflag = zeros(nInterfaceElem,1);
duf = zeros(nInterfaceElem,1);

%% Initial conditions

Nt = 5000;
dt_c=3600;
dt = dt_c;
W (:,1) = params.W0;
uf(:,1) = params.uf0;
theta(:,1) = params.theta0;
P(:,1)=params.P0;
Time(1,1)=0;
dampingFactor = 1;
Ss(:,1)=-2e+7;
Sn(:,1)=-4e+7;
Ds(:,1)=-0.000602895532919860;
landa=0.6;
beta=0.3025;

%% Assembling Process for PoroElasticity
globalM = MassMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType, omega, domain_length, nex);
globalC = InertialMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType, omega, domain_length, nex);

globalK = StiffnessMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType);
globalH = ConductanceMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType);
globalS = StorageMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType);
globalQ = CouplingMatrix(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, elementType,omega, domain_length, nex);
globalDa = dampingBoundary(boundaryConditions, node_coordinates, connectivity, NODES, elementType, params);


globalQ_inter = CouplingInterface(connectivity, node_coordinates, PSI, NODES, CRACK, n, elementType, PHI, omega, domain_length, nex);
globalFu_ext1 = TractionForces(boundaryConditions, node_coordinates, connectivity, NODES, elementType);
globalFu_ext2 = BodyForceSolid(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType);
globalFu_ext3 = InsituStress(params, connectivity, node_coordinates, PSI, PHI, NODES, CRACK, omega, domain_length, nex, elementType);
globalFp_ext1 = BoundaryFluidSource(boundaryConditions, node_coordinates, connectivity, NODES, elementType);
globalFp_ext2 = BodyForceFluid(params, node_coordinates, connectivity, PSI, PHI, NODES, CRACK, elementType);
globalFp_ext3 = InjectionSource(boundaryConditions, NODES);
globalG_inter = Lagrange(connectivity, node_coordinates, NODES, CRACK, n, interface_elements, elementType, PHI, omega, domain_length, nex);
globalG_inter2 = StabLagrange(params,interface_elements,interface_nodes);
globalQ = globalQ + globalQ_inter;
globalFu = globalFu_ext1 + globalFu_ext2 - globalFu_ext3;


omega1=pi*(2*1-1)*params.cs/(2*domain_width);
omega3=pi*(2*3-1)*params.cs/(2*domain_width);
a_ray=2*omega1*omega3*0.04/(omega1+omega3);
b_ray=2*0.04/(omega1+omega3);
globalD = a_ray*globalM+b_ray*globalK + globalDa;
% globalD = 0*globalM+0*globalK + globalDa;
% Create full lagrange included stiffness matrix
DOF_lag = size(globalG_inter, 2);
globalK = [globalK, sparse(DOF_u,DOF_lag);sparse(DOF_lag, DOF_u),sparse(DOF_lag,DOF_lag)];
globalK(1:DOF_u, DOF_u+1:end)     = globalG_inter;
globalK(DOF_u+1:end, 1:DOF_u)     = globalG_inter';
globalK(DOF_u+1:end, DOF_u+1:end) = -globalG_inter2;
globalFu = [globalFu; sparse(DOF_lag, 1)];
globalM = [globalM, sparse(DOF_u,DOF_lag);sparse(DOF_lag, DOF_u),sparse(DOF_lag,DOF_lag)];
globalD = [globalD, sparse(DOF_u,DOF_lag);sparse(DOF_lag, DOF_u),sparse(DOF_lag,DOF_lag)];
globalQ = [globalQ;sparse(DOF_lag,size(globalQ,2))];
globalC = [globalC sparse(size(globalC,1),DOF_lag)];

U = zeros(DOF_u+DOF_lag, 1); 
U_ddot = zeros(DOF_u+DOF_lag, 1); 
U_dot = zeros(DOF_u+DOF_lag, 1); 
ResidualHistory = cell(Nt+1,1);

%% Apply Boundary conditions

fixedDOF_disp = boundaryConditions.fixedDOF_disp;
dispDOFs = fixedDOF_disp(:,1);
fixedPressureDOF = boundaryConditions.fixedPressureDOF;

if ~isempty(fixedPressureDOF)
    pressureDOFs = fixedPressureDOF(:,1) + DOF_u + DOF_lag;
    fixedPressures = fixedPressureDOF(:,2);
else
    pressureDOFs = [];
    fixedPressures = [];
    diagPressIdx = [];
end

%% Newton-Raphson iteration

tolerance = 1e-6;
maxIter   = 10;

for nt = 1:Nt

% Initial guess for next time step
    U(:,nt+1)  = U(:,nt);
    P(:,nt+1)  = P(:,nt);
    W(:,nt+1)  = W(:,nt);
    Sn(:,nt+1) = Sn(:,nt);
    Pn(:,nt+1) = Pn(:,nt);
    Ss(:,nt+1) = Ss(:,nt);
    uf(:,nt+1) = uf(:,nt);
    stickflag(:,nt+1) = stickflag(:,nt);
    Vs(:,nt+1) = Vs(:,nt);
    

    if all(Vs(:,nt+1) == 1e-9)
        dt_ev = inf;
    else
        dt_ev = min(0.1 * params.Dc ./ Vs(:,nt+1));
    end
    dt_candidate = min([max(1e-4, dt_ev), dt_c]);
    dt = min(dt_candidate, 1.2 * dt);

    iter = 0;
    Norm = Inf;

while Norm > tolerance && iter < maxIter

        iter = iter + 1;

        dK_dU = StiffnessInterface_Lagrange(params, connectivity, node_coordinates, NODES, CRACK, n, t, Sn(:,nt+1), duf, dt, stickflag(:,nt+1), elementType, PHI, omega, domain_length, nex);
        globalF_inter = FaultForce(connectivity, node_coordinates, NODES, CRACK, Pn(:,nt+1), Ss(:,nt+1), t, n, elementType, PHI, omega, domain_length, nex);
        [H_inter, S_inter, Q_inter, F_inter] = interface_flow(params, connectivity, node_coordinates, PSI, NODES, CRACK, t, W(:,nt+1), elementType);
        H = globalH + H_inter;
        S = globalS + S_inter;

        % Form residuals 
        Gu = globalFu + globalM * (1/(beta*dt^2)*U(:,nt) + 1/(beta*dt)*U_dot(:,nt) + (1/(2*beta)-1)*U_ddot(:,nt)) + ...
                        globalD * (landa/(beta*dt) * U(:,nt) + (landa/beta-1) * U_dot(:,nt) + dt*(landa/(2*beta)-1)*U_ddot(:,nt));
        Gp = (globalFp_ext1+globalFp_ext2+globalFp_ext3) + F_inter + globalC * (1/(beta*dt^2)*U(:,nt) + 1/(beta*dt)*U_dot(:,nt) + (1/(2*beta)-1)*U_ddot(:,nt)) + ...
                        (globalQ'+[Q_inter,sparse(size(Q_inter,1),DOF_lag)])* (landa/(beta*dt)*U(:,nt) + (landa/beta-1)*U_dot(:,nt) + dt * (landa/(2*beta)-1)*U_ddot(:,nt)) + S/dt * P(:,nt);
                      

        Residual  = [-landa/(beta*dt) * ((1/(beta*dt^2) * globalM + landa/(beta*dt) * globalD + globalK) * U(:,nt+1) - globalQ * P(:,nt+1) + [globalF_inter;zeros(DOF_lag, 1)] - Gu);...
                    (1/(beta*dt^2) * globalC + landa/(beta*dt) * (globalQ'+[Q_inter,sparse(size(Q_inter,1),DOF_lag)])) * U(:,nt+1) + (S/dt + H) * P(:,nt+1) - Gp];


        % Form jacobian
        K = globalK;
        K(1:DOF_u, 1:DOF_u) = K(1:DOF_u, 1:DOF_u) + dK_dU;
        Jacobian = [-landa/(beta*dt) * ((1/(beta*dt^2) * globalM + landa/(beta*dt) * globalD + K)), landa/(beta*dt) * globalQ; ...
                    (1/(beta*dt^2) * globalC + landa/(beta*dt) * (globalQ'+[Q_inter,sparse(size(Q_inter,1),DOF_lag)])), (1/dt*S + H)];


        % Apply essential boundary conditions (displacement)
        Jacobian(dispDOFs, :) = 0;
        Jacobian(:, dispDOFs) = 0;
        diagDispIdx = sub2ind(size(Jacobian), dispDOFs, dispDOFs);
        Jacobian(diagDispIdx) = 1;
        Residual(dispDOFs) = 0;

        % Apply essential boundary conditions (pressure)
        if ~isempty(pressureDOFs)
        Jacobian(pressureDOFs, :) = 0;
        Jacobian(:, pressureDOFs) = 0;
        diagPressIdx = sub2ind(size(Jacobian), pressureDOFs, pressureDOFs);
        Jacobian(diagPressIdx) = 1;
        Residual(pressureDOFs) = P(fixedPressureDOF(:,1), nt+1) - fixedPressures;
        end


%--- 1) Equilibration (once) ---------------------------------------------
% if iter == 1
%     [~, R0, C0] = equilibrate(Jacobian);  
% end

Jt = sparse(Jacobian);  
rt = - Residual;
deltaX = Jt \ rt;

% Normalize Displacement and Pressure Corrections**
        ResidualU = norm(Residual(1:(DOF_u+DOF_lag),1)) / norm(globalFu);
        ResidualP = norm(Residual((DOF_u+DOF_lag)+1:end,1)) / (norm(globalFp_ext3)+1e-7);

        % Compute L2 Norm with Normalized Values**
        Norm = sqrt(norm(ResidualU, 2)^2 + norm(ResidualP, 2)^2);
        % ResidualHistory{nt+1,1}=[ResidualHistory{nt+1}; Norm];
        fprintf('Time step %d, Iteration %d: Norm = %e\n', nt, iter, Norm);

        % Update displacement and pressure for time step nt+1
        U(:, nt+1) = U(:, nt+1) + dampingFactor*deltaX(1:(DOF_u+DOF_lag));  % update displacement
        P(:, nt+1) = P(:, nt+1) + dampingFactor*deltaX((DOF_u+DOF_lag)+1:end);  % Update pressure

        Ux(:,nt+1) = U(1:2:2*nNodes,nt+1);
        Uy(:,nt+1) = U(2:2:2*nNodes,nt+1);
        a(:,nt+1) = U(2*nNodes+1:2*nNodes+2*nnz(NODES(:,2)),nt+1);
        LM(:,nt+1) = U(DOF_u+1:end,nt+1);
        Pp(:,nt+1) = P(1:nNodes,nt+1);

        [pn,pt,gn,gt,LM_elem] = Interface(params, connectivity, node_coordinates, NODES, CRACK, U(:,nt+1), n, t, interface_elements, LM(:,nt+1), elementType, omega, domain_length, nex, PHI);
        Dn(:,nt+1) = gn;
        Ds(:,nt+1) = gt;
        Sn(:,nt+1) = LM_elem;
        Ss(:,nt+1) = pt;
        Pn(:,nt+1) = 0;
        W(:,nt+1) = params.W0 + min(0, LM_elem)*params.Dn_max./((params.Kn)*params.Dn_max + min(0, LM_elem));
        % W(:,nt+1) = params.W0 + Sn(:,nt+1)/params.Kn;

        [Ss, slip, Vs, uf, duf, theta, stickflag] = updateInterface(params, Ds, Ss, Sn, Vs, theta, uf, slip, stickflag, nt, dt, iter);

        if Norm < tolerance
            fprintf('Convergence achieved at time step %d, iteration %d\n', nt, iter);
            break;
        end
end

U_ddot(1:DOF_u,nt+1) = 1/(beta*dt^2) * (U(1:DOF_u, nt+1)-U(1:DOF_u, nt)) - 1/(beta*dt) * U_dot(1:DOF_u,nt) - (1/(2*beta)-1) * U_ddot(1:DOF_u,nt);
U_ddotx(:,nt+1) = U_ddot(1:2:2*nNodes,nt+1);
U_ddoty(:,nt+1) = U_ddot(2:2:2*nNodes,nt+1);
U_ddot_total(:,nt+1) = sqrt(U_ddotx(:,nt+1).^2 + U_ddoty(:,nt+1).^2);
U_dot(1:DOF_u,nt+1) = landa/(beta*dt)*(U(1:DOF_u,nt+1)-U(1:DOF_u,nt)) - (landa/beta-1)*U_dot(1:DOF_u,nt) - dt * (landa/(2*beta)-1)*U_ddot(1:DOF_u,nt);
U_dotx(:,nt+1) = U_dot(1:2:2*nNodes,nt+1);
U_doty(:,nt+1) = U_dot(2:2:2*nNodes,nt+1);
U_dot_total(:,nt+1)  = sqrt(U_dotx(:,nt+1).^2  + U_doty(:,nt+1).^2);

Time(nt+1,1) = dt + Time(nt,1);
current_time_hours = Time(end) / 3600;
if current_time_hours >= 360
        fprintf('Reached total simulation time of 360 hours. Quitting simulation.\n');
        break;
end
if iter==maxIter
    break
end

if mod(nt, 500)==0

fprintf('save data at time step %d...\n', nt);
save('simulation_WW_SH.mat','-v7.3');

end
end

Qinj = boundaryConditions.concentratedSources(1,2);
InjectedVolume =  Qinj*Time;  % time is Nt × 1 in seconds

%% Calculate seismicity parameters

[AS, M0, Mw, CFS, SSD, ST] = SeismicityParameters(interface_elements, interface_nodes, Vs, slip, Ss, Sn, uf, params);

%% Calculate stresses 

[Sxx,Sxy,Syy,Svm,S1,S2,theta_p] = elemStress2(params, U, connectivity, node_coordinates, PSI, NODES, CRACK, omega, domain_length, nex, elementType, Nt);

%% Post-Processing and Visualization

% postProcess(node_coordinates, U(1:2*nNodes,end), Pp(:,end), connectivity, Sxx(:,end), Syy(:,end), Sxy(:,end), elementType);
save('simulation_WW_SH.mat','-v7.3');



