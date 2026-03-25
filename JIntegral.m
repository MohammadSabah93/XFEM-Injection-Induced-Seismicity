% Original version of J-Integral code written by Nguyen Vinh Phu (2006)
% Modified By: Matthew Jon Pais, University of Florida (2010)
% Website: http://sites.google.com/site/matthewjpais/Home
% Email: mpais@ufl.edu, matthewjpais@gmail.com

function [KI, KII] = JIntegral(params, connectivity, node_coordinates, omega, DISPLACEMENT, PSI, NODES, CRACK, domain_length, nex, elementType)
% This function calculates the mixed-mode stress intensity factors for a
% cracked body using the domain form of the interaction integrals based on
% the original formulation by Rice.

nElem  = size(connectivity,1);                                               % Total number of elements
nNode  = max(connectivity(:));                                       % Total number of nodes
lXElem = domain_length/nex;                                                         % Length of elements in the x-direction
Em     = params.Yu;                                                            % Young's modulus for the matrix
vm     = params.nu;                                                            % Poisson's ratio for the matrix
plane  = 2;                                                            % Plane stress or plane strain
Gm     = Em/2/(1+vm);                                                       % Shear modulus for the matrix
XYZ = node_coordinates;

% Create elastic constant matrix
if plane == 1                                                               % Plane stress
    t = MAT(6);                                                             % Plane stress thickness
    
    C1 = Em/(1-vm^2);                                                       % Constant for elastic constant matrix
    C2 = Em*vm/(1-vm^2);                                                    % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm = t*[C1 C2  0;...
            C2 C1  0;...
             0  0 C3];

    km = (3-vm)/(1+vm);
 
elseif plane == 2                                                           % Plane strain
    C1 = Em*(1-vm)/(1+vm)/(1-2*vm);                                         % Constant for elastic constant matrix
    C2 = Em*vm/(1+vm)/(1-2*vm);                                             % Constant for elastic constant matrix
    C3 = Em/2/(1+vm);                                                       % Constant for elastic constant matrix
    Cm = [C1 C2  0;...
          C2 C1  0;...
           0  0 C3];
    km = 3-4*vm;
end

nCT = length(omega);                                                        % Number of crack tips in domain

for iJ = 1:nCT

    if iJ == 1
        nPt = size(CRACK,1);                                                % Determine number of data points defining crack
        xCT = CRACK(nPt,1);                                                 % X-coordinate of crack tip
        yCT = CRACK(nPt,2);                                                 % Y-coordinate of crack tip
    else
        xCT = CRACK(1,1);                                                   % X-coordinate of crack tip
        yCT = CRACK(1,2);                                                   % Y-coordinate of crack tip
    end
    
    CCS = [cos(omega(iJ)) sin(omega(iJ));-sin(omega(iJ)) cos(omega(iJ))];
    
    % Geometric predicates for crack geometry
    area  = lXElem*lXElem;                                                  % Elemental area
    
    % Define elements to be used in J-integral
    c      = 1;                                                             % Magnification factor
    radius = c*sqrt(area);                                                  % Search radius
    
    dist = zeros(1,nNode);
    for iN = 1:nNode
        Xn       = XYZ(iN,1);                                               % X-coordinate for the current node
        Yn       = XYZ(iN,2);                                               % Y-coordinate for the current node
        X        = Xn-xCT;                                                  % Horizontal distance from crack tip to current node
        Y        = Yn-yCT;                                                  % Vertical distance from crack tip to current node
        XYloc    = CCS*[X Y]';                                              % Change to crack tip coordinates
        r        = sqrt(XYloc(1)^2+XYloc(2)^2);                             % Radius from crack tip to current gauss point
        dist(iN) = r;                                                       % Store radius value
    end
    
    % Determine elements in J-integral and assign nodal q values
    temp    = dist-radius;                                                  % Determine whether or not the node is outside the search radius
    temp    = temp(connectivity(:,1:4))';                                         % Build elemental distance vector
    Jdomain = NaN(1,nElem);                                                 % Initialize Jdomain to NaN
    index   = 1;                                                            % Index to track Jdomain
    for i = 1:nElem
        if (min(temp(:,i)) < 0)
            Jdomain(index) = i;
            index = index+1;
        end
    end
    
    Jdomain(isnan(Jdomain)) = [];                                           % Remove the unused locations in Jdomain
    temp  = dist-radius;                                                    % Determine whether or not the node is outside the search radius
    temp  = temp(connectivity(Jdomain,1:4))';                                     % Build elemental distance vector for elements in J-integral domain
    temp  = (temp<=0);                                                      % Calculate nodes inside/outside search radius
    qNode = temp';                                                          % Store the nodal q values used in J-integral calculations
    
    I = zeros(2,1);                                                         % Interaction integral for combined states 1 and 2
    for iElem = 1:length(Jdomain)
        nElem = Jdomain(iElem);                                             % The global number for the current element
        N1  = connectivity(nElem,1);                                              % Node 1 for current element
        N2  = connectivity(nElem,2);                                              % Node 2 for current element
        N3  = connectivity(nElem,3);                                              % Node 3 for current element
        N4  = connectivity(nElem,4);                                              % Node 4 for current element
        NN  = NODES([N1 N2 N3 N4]',:);                                      % Nodal data for current element
        CTN = nnz(NN(:,4));                                                 % Number of nodes with crack tip enrichment
        HEN = nnz(NN(:,2));                                                 % Number of nodes with Heaviside enrichment
        NEN = HEN+CTN;                                                  % Number of enriched nodes
        
        X1 = XYZ(N1,1); X2 = XYZ(N2,1); X3 = XYZ(N3,1); X4 = XYZ(N4,1);     % X-coordinates of nodes
        Y1 = XYZ(N1,2); Y2 = XYZ(N2,2); Y3 = XYZ(N3,2); Y4 = XYZ(N4,2);     % Y-coordinates of nodes
        
        xyz = [X1 Y1;X2 Y2;X3 Y3;X4 Y4];                                    % Nodal coordinate matrix
        
        % Define the gauss points and the gauss weights
        if NEN == 4                                                         % Fully enriched element
            PN = [ PSI(N1)  PSI(N2)  PSI(N3)  PSI(N4)];                     % Nodal crack level set values
                      
             if  CTN == 4 
             [gp, gw, J] = subDomain(13, PN, xyz, 1, 1, CCS,CRACK);
             elseif HEN == 4
             [gp, gw, J] = subDomain(7, PN, xyz, 0, 1, CCS,CRACK);
             else
             [gp, gw, J] = subDomain(7, PN, xyz, 0, 1, CCS,CRACK);
             end

        else
            [gp, gw] = gaussPoints(elementType, 2);
            J = [];
        end
        
        % Loop through gauss points in current element to solve for J-integral
        Uenr = []; iGP = 1; iLoc = 1;
        for i = 1:length(gp)
            xi = gp(i,1); eta = gp(i,2);
            W  = gw(i);
            
            if isempty(J) == 0
                Ji   = [J(i,1) J(i,2);J(i,3) J(i,4)];                       % Jacobian for current subdivision
                detJ = det(Ji);                                             % Determinant of Jacobian for current subdivision
            else
                xyz1 = CCS*([X1 Y1]'-[xCT yCT]');                           % Crack tip coordinates of node 1
                xyz2 = CCS*([X2 Y2]'-[xCT yCT]');                           % Crack tip coordinates of node 2
                xyz3 = CCS*([X3 Y3]'-[xCT yCT]');                           % Crack tip coordinates of node 3
                xyz4 = CCS*([X4 Y4]'-[xCT yCT]');                           % Crack tip coordinates of node 4
                
                dxpdxi  = 1/4*(-(1-eta)*xyz1(1)+(1-eta)*xyz2(1)+(1+eta)*xyz3(1)-(1+eta)*xyz4(1));
                dxpdeta = 1/4*(-(1- xi)*xyz1(1)-(1+ xi)*xyz2(1)+(1+ xi)*xyz3(1)+(1- xi)*xyz4(1));
                dypdxi  = 1/4*(-(1-eta)*xyz1(2)+(1-eta)*xyz2(2)+(1+eta)*xyz3(2)-(1+eta)*xyz4(2));
                dypdeta = 1/4*(-(1- xi)*xyz1(2)-(1+ xi)*xyz2(2)+(1+ xi)*xyz3(2)+(1- xi)*xyz4(2));
                
                Je   = [dxpdxi dypdxi;dxpdeta dypdeta];                     % Jacobian for current element
                detJ = det(Je);                                             % Determinant of the Jacobian
            end
            
            % Define quadrilateral shape functions and derivatives
            N  = 1/4*[(1-xi)*(1-eta);(1+xi)*(1-eta);(1+xi)*(1+eta);(1-xi)*(1+eta)];
            Nx = 2/lXElem*1/4*[-(1-eta);1-eta;1+eta;-(1+eta)];
            Ny = 2/lXElem*1/4*[-(1-xi);-(1+xi);1+xi;1-xi];
            
            Xgp = N(1)*X1+N(2)*X2+N(3)*X3+N(4)*X4;                          % The global X for the current gauss point
            Ygp = N(1)*Y1+N(2)*Y2+N(3)*Y3+N(4)*Y4;                          % The global Y for the current gauss point
            
            X     = Xgp-xCT;                                                % Horizontal distance from crack tip to gauss point
            Y     = Ygp-yCT;                                                % Vertical distance from crack tip to gauss point
            XYloc = CCS*[X Y]';                                             % Change to crack tip coordinate system
            r     = sqrt(XYloc(1)^2+XYloc(2)^2);                            % Radius from crack tip to current gauss point
            theta = atan2(XYloc(2),XYloc(1));                               % Angle from crack tip to current gauss point
            
            U = [DISPLACEMENT(2*N1-1) DISPLACEMENT(2*N1) DISPLACEMENT(2*N2-1) DISPLACEMENT(2*N2)...
                DISPLACEMENT(2*N3-1) DISPLACEMENT(2*N3) DISPLACEMENT(2*N4-1) DISPLACEMENT(2*N4)];
            
            
            Benr = [];
            Bu = [Nx(1)   0   Nx(2)   0   Nx(3)   0   Nx(4)   0;...
                    0   Ny(1)   0   Ny(2)   0   Ny(3)   0   Ny(4);...
                  Ny(1) Nx(1) Ny(2) Nx(2) Ny(3) Nx(3) Ny(4) Nx(4)];
            
            iB = 1;
            for iN = 1:4
                if NN(iN,2) ~= 0                                            % Heaviside node
                    psi1 = PSI(N1);                                         % Value of phi at Node 1
                    psi2 = PSI(N2);                                         % Value of phi at Node 2
                    psi3 = PSI(N3);                                         % Value of phi at Node 3
                    psi4 = PSI(N4);                                         % Value of phi at Node 4
                    psi  = N(1)*psi1+N(2)*psi2+N(3)*psi3+N(4)*psi4;         % Calculate phi at the current gauss point
                    Hgp  = sign(psi);                                       % Find the value of H at the current gauss point
                    
                    Hi = NN(iN,3);                                          % Find the value of H at the current node
                    H  = Hgp-Hi;                                            % Find the value of H
                    
                    Ba = [Nx(iN)*H     0;
                              0    Ny(iN)*H;
                          Ny(iN)*H Nx(iN)*H];
                    Benr(:,iB:(iB+1)) = Ba;
                    iB = iB+2;
                    
                    if iGP == 1
                        Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,2)-1) DISPLACEMENT(2*NN(iN,2))];
                        iLoc = iLoc+2;
                    end
                
                elseif NN(iN,4) ~= 0                                        % Crack tip node
                    a1gp = sqrt(r)*sin(theta/2);                            % Node 1 crack tip enrichment value
                    a2gp = sqrt(r)*cos(theta/2);                            % Node 2 crack tip enrichment value
                    a3gp = sqrt(r)*sin(theta)*sin(theta/2);                 % Node 3 crack tip enrichment value
                    a4gp = sqrt(r)*sin(theta)*cos(theta/2);                 % Node 4 crack tip enrichment value
                    
                    a1N = NN(iN,5); a2N = NN(iN,7); a3N = NN(iN,9); a4N = NN(iN,11);
                    a1  = a1gp-a1N; a2  = a2gp-a2N; a3  = a3gp-a3N; a4  = a4gp-a4N;
                    
                    c = 1/2/sqrt(r); ct = CCS(1,1); st = CCS(1,2);
                    
                    % Derivate of Phi_alpha with respect to X
                    Px = c*[-sin(theta/2)*ct              + cos(theta/2)*-st;...
                             cos(theta/2)*ct              + sin(theta/2)*-st;...
                            -sin(3*theta/2)*sin(theta)*ct + (sin(theta/2)+sin(3*theta/2)*cos(theta))*-st;...
                            -cos(3*theta/2)*sin(theta)*ct + (cos(theta/2)+cos(3*theta/2)*cos(theta))*-st];
                    
                    % Derivative of Phi_alpha with respect to Y
                    Py = c*[-sin(theta/2)*st              + cos(theta/2)*ct;...
                             cos(theta/2)*st              + sin(theta/2)*ct;...
                            -sin(3*theta/2)*sin(theta)*st + (sin(theta/2)+sin(3*theta/2)*cos(theta))*ct;...
                            -cos(3*theta/2)*sin(theta)*st + (cos(theta/2)+cos(3*theta/2)*cos(theta))*ct];
                    
                    B1 = [Nx(iN)*a1+N(iN)*Px(1)          0;
                                   0            Ny(iN)*a1+N(iN)*Py(1);
                          Ny(iN)*a1+N(iN)*Py(1) Nx(iN)*a1+N(iN)*Px(1)];
                    
                    B2 = [Nx(iN)*a2+N(iN)*Px(2)          0;
                                   0            Ny(iN)*a2+N(iN)*Py(2);
                          Ny(iN)*a2+N(iN)*Py(2) Nx(iN)*a2+N(iN)*Px(2)];
                    
                    B3 = [Nx(iN)*a3+N(iN)*Px(3)          0;
                                   0            Ny(iN)*a3+N(iN)*Py(3);
                          Ny(iN)*a3+N(iN)*Py(3) Nx(iN)*a3+N(iN)*Px(3)];
                    
                    B4 = [Nx(iN)*a4+N(iN)*Px(4)          0;
                                   0            Ny(iN)*a4+N(iN)*Py(4);
                          Ny(iN)*a4+N(iN)*Py(4) Nx(iN)*a4+N(iN)*Px(4)];
                    
                    Bb = [B1 B2 B3 B4];
                    Benr(:,iB:(iB+7)) = Bb;
                    iB = iB+8;
                    
                    if iGP == 1
                        index = 3;
                        for iAlpha = 1:4
                            Uenr(iLoc:(iLoc+1)) = [DISPLACEMENT(2*NN(iN,iAlpha+index)-1) DISPLACEMENT(2*NN(iN,iAlpha+index))];
                            index = index+1;
                            iLoc  = iLoc+2;
                        end
                    end
                end
            end

            C = Cm; G = Gm; k = km; 
            
            B  = [Bu Benr];
            Xe = [U Uenr]';
            lB = size(B,2);
            
            % Solve for the stress and strain at the current gauss point
            strain = B*Xe;                                                  % Strain at current gauss point
            stress = C*strain;                                              % Stress at current gauss point
            
            % Derivates of q with respect to X and Y
            q     = qNode(iElem,:);                                         % q values at nodes of current element
            gradq = q*[Nx(1) Ny(1);Nx(2) Ny(2);Nx(3) Ny(3);Nx(4) Ny(4)];    % The derivative of q with respect to X and Y
            
            % Derivatives of nodal displacements with respect to X and Y
            Ux = B(1,1:2:lB)*Xe(1:2:lB);                                    % Derivative of U with respect to X
            Uy = B(2,2:2:lB)*Xe(1:2:lB);                                    % Derivative of U with respect to Y
            Vx = B(1,1:2:lB)*Xe(2:2:lB);                                    % Derivative of V with respect to X
            Vy = B(2,2:2:lB)*Xe(2:2:lB);                                    % Derivative of V with respect to Y
            
            % Convert quanties into crack tip coordinate system
            GradQ     = CCS*gradq';                                         % Gradient of q in crack tip coordinate system
            GradDisp  = CCS*[Ux Uy;Vx Vy]*CCS';                             % Gradient of displacement in crack tip coordinate system
            CalStress = CCS*[stress(1) stress(3);stress(3) stress(2)]*CCS'; % Stresses in crack tip coordinate system
            
            K1 = 1.0;                                                       % Chosen KI for pure Mode I asymptotic field
            K2 = 1.0;                                                       % Chosen KII for pure Mode II asymptotic field
                                                 
                % Predefine variables to be used repeatidly
                SQR  = sqrt(r);                                             % The square-root of r
                CT   = cos(theta);                                          % The cosine of theta
                ST   = sin(theta);                                          % The sine of theta
                CT2  = cos(theta/2);                                        % The cosine of one half theta
                ST2  = sin(theta/2);                                        % The sine of one half theta
                C3T2 = cos(3*theta/2);                                      % The cosine of three halves theta
                S3T2 = sin(3*theta/2);                                      % The sine of three halves theta
                
                drdx =  CT;                                                 % Derivative of r with respect to X
                drdy =  ST;                                                 % Derivative of r with respect to Y
                dtdx = -ST/r;                                               % Derivative of theta with respect to X
                dtdy =  CT/r;                                               % Derivative of theta with respect to Y
                
                cAuxStress = sqrt(1/(2*pi));                                % Constant for auxiliary stress calculation
                cAuxDisp   = sqrt(1/(2*pi))/(2*G);                          % Constant for auxiliary displacement calculation
                
                AuxStress   = zeros(2,2);                                   % Initialize auxiliary stress matrix
                AuxGradDisp = zeros(2,2);                                   % Initialize gradient of displacement matrix
                for mode = 1:2
                    if mode == 1                                            % K1 = 1.0 and K2 = 0.0
                        AuxStress(1,1) = K1*cAuxStress/SQR*CT2*(1-ST2*S3T2);% Auxiliary stress component for Mode I loading in x-direction
                        AuxStress(2,2) = K1*cAuxStress/SQR*CT2*(1+ST2*S3T2);% Auxiliary stress component for Mode I loading in y-direction
                        AuxStress(1,2) = K1*cAuxStress/SQR*ST2*CT2*C3T2;    % Auxiliary shear stress component for Mode I loading
                        AuxStress(2,1) = AuxStress(1,2);                    % Auxiliary shear stress component for Mode I loading
                        
                        du1dr = K1*cAuxDisp*0.5/SQR*CT2*(k-CT);             % Derivative of auxiliary x-displacement with respect to r
                        du1dt = K1*cAuxDisp*SQR*(-0.5*ST2*(k-CT)+CT2*ST);   % Derivative of auxiliary x-displacement with respect to theta
                        du2dr = K1*cAuxDisp*0.5/SQR*ST2*(k-CT);             % Derivative of auxiliary y-displacement with respect to r
                        du2dt = K1*cAuxDisp*SQR*(0.5*CT2*(k-CT)+ST2*ST);    % Derivative of auxiliary y-displacement with respect to theta
                        
                        AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx;           % Auxiliary displacement gradient of u1 with respect to x
                        AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy;           % Auxiliary displacement gradient of u1 with respect to y
                        AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx;           % Auxiliary displacement gradient of u2 with respect to x
                        AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy;           % Auxiliary displacement gradient of u2 with respect to y
                        
                        AuxStrain(1,1) = AuxGradDisp(1,1);
                        AuxStrain(1,2) = 1/2*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
                        AuxStrain(2,1) = AuxStrain(1,2);
                        AuxStrain(2,2) = AuxGradDisp(2,2);
                    elseif mode == 2                                         % K1 = 0.0 and K2 = 1.0
                        AuxStress(1,1) = -K2*cAuxStress/SQR*ST2*(2+CT2*C3T2);% Auxiliary stress component for Mode II loading in x-direction
                        AuxStress(2,2) =  K2*cAuxStress/SQR*ST2*CT2*C3T2;    % Auxiliary stress component for Mode II loading in y-direction
                        AuxStress(1,2) =  K2*cAuxStress/SQR*CT2*(1-ST2*S3T2);% Auxiliary shear stress component for Mode II loading
                        AuxStress(2,1) =  AuxStress(1,2);                    % Auxiliary shear stress component for Mode II loading
                        
                        du1dr =  K2*cAuxDisp*0.5/SQR*ST2*(k+2+CT);          % Derivative of auxiliary x-displacement with respect to r
                        du1dt =  K2*cAuxDisp*SQR*(0.5*CT2*(k+2+CT)-ST2*ST); % Derivative of auxiliary x-displacement with respect to theta
                        du2dr = -K2*cAuxDisp*0.5/SQR*CT2*(k-2+CT);          % Derivative of auxiliary y-displacement with respect to r
                        du2dt = -K2*cAuxDisp*SQR*(-0.5*ST2*(k-2+CT)-CT2*ST);% Derivative of auxiliary y-displacement with respect to theta
                        
                        AuxGradDisp(1,1) = du1dr*drdx+du1dt*dtdx;           % Auxiliary displacement gradient of u1 with respect to x
                        AuxGradDisp(1,2) = du1dr*drdy+du1dt*dtdy;           % Auxiliary displacement gradient of u1 with respect to y
                        AuxGradDisp(2,1) = du2dr*drdx+du2dt*dtdx;           % Auxiliary displacement gradient of u2 with respect to x
                        AuxGradDisp(2,2) = du2dr*drdy+du2dt*dtdy;           % Auxiliary displacement gradient of u2 with respect to y
                        
                        AuxStrain(1,1) = AuxGradDisp(1,1);
                        AuxStrain(1,2) = 1/2*(AuxGradDisp(1,2)+AuxGradDisp(2,1));
                        AuxStrain(2,1) = AuxStrain(1,2);
                        AuxStrain(2,2) = AuxGradDisp(2,2);
                    end
                    
                    I1 = (CalStress(1,1)*AuxGradDisp(1,1)+CalStress(2,1)*AuxGradDisp(2,1))*GradQ(1)+...
                        (CalStress(1,2)*AuxGradDisp(1,1)+CalStress(2,2)*AuxGradDisp(2,1))*GradQ(2);
                    
                    I2 = (AuxStress(1,1)*GradDisp(1,1)+AuxStress(2,1)*GradDisp(2,1))*GradQ(1)+...
                        (AuxStress(1,2)*GradDisp(1,1)+AuxStress(2,2)*GradDisp(2,1))*GradQ(2);
                    
                    if (isnan(I1) == 1) || (isnan(I2) == 1), break, end
                    
                    StrainEnergy = 0;
                    for j = 1:2
                        for k = 1:2
                            StrainEnergy = StrainEnergy+CalStress(j,k)*AuxStrain(j,k);
                        end
                    end
                    
                    I(mode,1) = I(mode,1)+(I1+I2-StrainEnergy*GradQ(1))*detJ*W;
                end
                iGP = iGP + 1;
         
               
                iGP = iGP + 1;
        end
    end
    
    
        % Find the effective modulus
        if plane == 1                                                       % Plane stress
            Eeff = Em;
        elseif plane == 2
            Eeff  = Em/(1-vm^2);                                            % Plane strain
        end
        
        % Solve for the mixed-mode stress intensity factors
        Kcalc   = I*Eeff/2;                                                 % Solve for Mode I and Mode II SIF
        KI(iJ)  = Kcalc(1);                                                 % Mode I  SIF
        KII(iJ) = Kcalc(2);                                                 % Mode II SIF
    
end