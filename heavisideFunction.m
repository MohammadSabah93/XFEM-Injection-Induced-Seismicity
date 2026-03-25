function H = heavisideFunction(phi, type)
    % Inputs:
    %   phi  - Level set function (scalar or vector).
    %   type - Type of Heaviside function ('type1' or 'type2').
    %
    % Outputs:
    %   H    - Heaviside function values corresponding to phi.

    if strcmp(type, 'type1')
        % Heaviside function: H(x) = 0 if phi(x)<0, 1 if phi(x)>0
        H = zeros(size(phi,2),1); % Initialize H
        H(phi >= 0) = 1;       % Assign 1 where phi > 0

    elseif strcmp(type, 'type2')
        % Heaviside function: H(x) = -1 if phi(x)<0, +1 if phi(x)>0
        H = zeros(size(phi,2),1); % Initialize H
        H(phi < 0) = -1;      % Assign -1 where phi < 0
        H(phi > 0) = 1;       % Assign +1 where phi > 0

    else
        error('Invalid type. Use "type1" or "type2".');
    end
end
