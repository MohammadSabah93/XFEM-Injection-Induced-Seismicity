
function [S] = slipfunction(params, Ss_trial, slipfix, thetafix, Sn, dt, damp, counter)

Ks=params.Ks;
Kn=params.Kn;
V0=params.V0;
Dc=params.Dc(counter);
a=params.a(counter);
b=params.b(counter);
uf0=params.uf0;
S=slipfix;
iter = 0;
Norm = Inf;
tolerance = 1e-8;
maxIter = 50;
while Norm > tolerance && iter < maxIter

iter = iter + 1;

V = (S-slipfix)/dt;

if abs(V)<V0
  V=V0; 
end
Vabs = abs(V);

theta = (thetafix + dt) / (1 + (dt * Vabs) / Dc);   
uf = uf0 + a * log(Vabs / V0) + b * log((V0 * theta) / Dc);  
r = abs(Ss_trial - Ks * sign(Ss_trial) * (S - slipfix)) - uf * abs(Sn) - damp * V; 
duf = (1./dt) .* ( a / V - (b / theta) * ( Dc * dt * (theta + dt) ) / ( (Dc + dt * V)^2 ) ); 
J = Ks - duf * Sn + damp / dt;
S = S + r / J;

if (abs(norm(r))<=tolerance)
    break;
end
end

