function params = defineModelParameters()
    % Define material and model properties as a structured variable for clarity

    % Solid properties
    params.nu = 0.3;                % Poisson's ratio
    params.Yu = 50e+9;             % Young's modulus (Pa)
    params.rho_s = 2500;           % Rock density (kg/m^3)
    params.Cs = 900;               % Specific heat capacity of rock (J/(kg·K))
    params.fi = 0.1;               % Porosity
    params.K_perm = 2e-15;         % Intrinsic permeability (m^2)
    params.landa_s = 3;            % Rock thermal conductivity (W/(m·K))
    params.Bs = 9e-8;              % Volumetric thermal expansion of solid (1/K)
    params.com_s = 2e-11;          % Compressibility of solid (1/Pa)

    % Fluid properties
    params.mu = 0.001;             % Fluid viscosity (Pa·s)
    params.Mo = params.K_perm / params.mu; % Mobility (m^2/(Pa·s))
    params.com_f = 5e-9;        % Water compressibility (1/Pa)
    params.Cf = 4200;              % Specific heat capacity of fluid (J/(kg·K))
    params.alpha= 0.75;            % Biot-Willis coefficient
    params.rho_f = 1000;           % Fluid density (kg/m^3)
    params.landa_f = 0.6;          % Fluid thermal conductivity (W/(m·K))
    params.Bf = 1e-3;              % Volumetric thermal expansion of fluid (1/K)

    % Fracture properties
    params.Sc_f = 4.5e-10;         % Storage coefficient of fracture (1/Pa)
    params.Kn = 2e+10;            % Initial normal stiffness (Pa/m)
    params.Ks = 2e+10;             % Shear stiffness (Pa/m)
    params.Tw = 1e-7;              % Transversal conductivity coefficient (m/s)
    params.Dn_max = -1e-3;         % Maximum discontinuity closure (m)
    params.hc = 0.6;               % Transversal thermal conductivity of fracture (W/(m·K))
    params.W0 = 1e-3;              % fracture initial width

    % Other properties
    params.G = params.Yu / (2 * (1 + params.nu)); % Shear modulus (Pa)
    params.Cm = (1 - params.fi) * params.Cs + params.fi * params.Cf; % Specific heat capacity of porous media
    params.rho_m = (1 - params.fi) * params.rho_s + params.fi * params.rho_f; % Bulk density of porous media
    params.landa_m = (1 - params.fi) * params.landa_s + params.fi * params.landa_f; % Thermal conductivity of porous media
    params.TDC_m = params.landa_m / (params.rho_m * params.Cm); % Thermal diffusivity coefficient of porous media
    params.TDC_f = params.landa_f / (params.rho_f * params.Cf); % Thermal diffusivity coefficient of fluid
    params.n = (params.rho_f * params.Cf) / (params.rho_m * params.Cm); % Ratio of specific heat capacities
    params.landa2 = (params.alpha - params.fi) * params.Bs + params.fi * params.Bf; % Combined thermal expansion coefficient
    params.K = 2*params.G*(1+params.nu)/(3*(1-2*params.nu)); 
    params.Sc_m = params.fi * params.com_f + (params.alpha - params.fi) * (1 - params.alpha) / params.K; % Storage coefficient of porous media (1/Pa)

    % Dynamic properties
    params.cs = sqrt(params.G / params.rho_s); % Shear wave speed
    params.cp = sqrt(params.Yu * (1 - params.nu) / (params.rho_s * (1 + params.nu) * (1 - 2 * params.nu))); % Compressional wave speed
    params.DF_shear = 0.5 * params.rho_s * params.cs; % Damping factor for shear
    params.DF_normal = 0.5 * params.rho_s * params.cp; % Damping factor for normal
    params.natural_frequencies = [0.5; 15]; % Target frequencies (Hz)
    params.omega = 2 * pi * params.natural_frequencies; % Angular frequencies (rad/s)
    params.zeta = [0.05; 0.05]; % Damping ratio

    % Fault frictional properties
    params.uf0 = 0.6;
    params.a = 0.01;
    params.b = 0.02;
    params.Dc = 1e-4;
    params.V0 = 1e-9;
    params.theta0 = params.Dc / params.V0;
    params.epsilon = 0.001;
    params.dt_min = 1e-6;
    params.cs=(params.G/params.rho_s)^0.5;   % Shear Wave Speed
    params.cp=sqrt(params.Yu*(1-params.nu)/((1+params.nu)*(1-2*params.nu)*params.rho_m));
    % params.damp=params.G/(2*params.cs);          % Damping Factor
    params.damp=0;
    % In-situ stresses

    params.Sxxr=-60e+6;
    params.Syyr=-20e+6;
    params.Sxyr=0;
    params.P0=0;
    params.g=[0;0];
end