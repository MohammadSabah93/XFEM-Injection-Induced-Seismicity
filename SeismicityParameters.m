function [AS, M0, Mw, CFS, SSD, ST] = SeismicityParameters(interface_elements, interface_nodes, Vs, slip, Ss, Sn, uf, params)

    % Parameters
    Fault_Width = 1;                  % [m]
    V0 = 0.1;                        % [m/s] reference slip velocity
    Nt = size(Vs,2);              % number of time steps
    numElems = size(interface_elements, 1);
    element_length = zeros(numElems, 1);  % [m]
   
    % Compute element lengths
    for e = 1:numElems
        n1 = interface_elements(e, 1);
        n2 = interface_elements(e, 2);
        coord1 = interface_nodes(n1, :);
        coord2 = interface_nodes(n2, :);
        element_length(e) = norm(coord2 - coord1);
    end

    % Slip increment per time step
    slip_inc = [zeros(numElems, 1), diff(slip, 1, 2)];  % size: numElems × Nt

    % Initialize outputs
    M0 = zeros(Nt, 1);
    Mw = NaN(Nt, 1);   % pre-fill with NaN
    AS = zeros(Nt, 1);
    total_area = sum(element_length) * Fault_Width;

    % Loop over time steps to compute M0 and Mw
    M0_total = 0;

    for t = 1:Nt

    total_slip_t = sum(slip(:, t) .* element_length);
    AS(t) = total_slip_t / total_area;

        if all(abs(Vs(:, t)) <= V0)
            % Reset moment accumulation
            M0_total = 0;
            M0(t) = 0;
            Mw(t) = NaN;
        else
            % Start or continue moment accumulation
            dslip = slip_inc(:, t);  % incremental slip at time t
            dM0 = sum(dslip .* element_length) * Fault_Width * params.G;
            M0_total = M0_total + dM0;
            M0(t) = M0_total;
            Mw(t) = (2/3) * (log10(max(M0(t), eps)) - 9.05);
        end
    end

    % Identify seismic events (same logic as before)
    active = any(Vs > V0, 1);
    d = diff([0, active, 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    numEvents = numel(starts);
    ST = [starts', ends'];  % start and end indices for each event

    % Event-based outputs
    SSD = zeros(numEvents, 1);

    for k = 1:numEvents
        
        % Average stress drop
        tau_start = abs(Ss(:, starts(k)));
        tau_end   = abs(Ss(:, ends(k)));
        SSD(k) = sum((tau_start - tau_end) .* element_length) / sum(element_length);

    end

    % Compute Coulomb Failure Stress (element × time)
    CFS = abs(Ss) - abs(uf .* Sn);

    
end
