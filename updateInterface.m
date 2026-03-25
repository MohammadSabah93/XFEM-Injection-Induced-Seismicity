function [Ss,slip,V, uf, duf, theta, stickflag]=updateInterface(params, Ds, Ss, Sn, V, theta, uf, slip, stickflag, nt, dt, iter)

counter = 1;
Ks=params.Ks;
V0=params.V0;
Dc=params.Dc;
a=params.a;
b=params.b;
uf0=params.uf0;
damp=params.damp;
nElements = size(Ds,1);

for iElem = 1:nElements

Ss_trial = Ss(counter,nt) + Ks * (Ds(counter,nt+1)-Ds(counter,nt)); % Stick predictor

if iter==1
if abs(V(counter,nt))<V0
V(counter,nt)=V0;
end
mu = uf0 + a(counter)*log(abs(V(counter,nt))/V0) + b(counter)*log(V0*theta(counter,nt)/Dc(counter));
F = abs(Ss_trial) - mu * abs(Sn(counter,nt)) - damp * abs(V(counter,nt));
else
if abs(V(counter,nt+1))<V0
V(counter,nt+1)=V0;
end
theta(counter,nt+1)=(theta(counter,nt)+dt)./(1+(dt*V(counter,nt+1)/Dc(counter)));
mu = uf0 + a(counter)*log(abs(V(counter,nt+1))/V0) + b(counter)*log(V0*theta(counter,nt+1)/Dc(counter));
F = abs(Ss_trial) - mu * abs(Sn(counter,nt+1)) - damp * abs(V(counter,nt+1));
end

if F<0
slip(counter,nt+1) = slip(counter,nt);
Ss(counter,nt+1) = Ss_trial;
stickflag(counter,nt+1) = 0;
else
% Slip corrector
slip(counter,nt+1) = slipfunction(params, Ss_trial, slip(counter,nt), theta(counter,nt), Sn(counter,nt+1), dt, damp, counter);
Ss(counter,nt+1) = Ss(counter,nt) + Ks * ((Ds(counter,nt+1)-Ds(counter,nt)) - sign(Ss_trial) * (slip(counter,nt+1)-slip(counter,nt)));
stickflag(counter,nt+1) = 1;
end

V(counter,nt+1)=(slip(counter,nt+1)-slip(counter,nt))/dt;
if abs(V(counter,nt+1))<V0
V(counter,nt+1)=V0;
else
V(counter,nt+1)=abs(V(counter,nt+1)); 
end

theta(counter,nt+1)=(theta(counter,nt)+dt)./(1+(dt*V(counter,nt+1)/Dc(counter)));
uf(counter,nt+1)= uf0 + a(counter)*log(V(counter,nt+1)/V0) + b(counter)*log(V0*theta(counter,nt+1)/Dc(counter));
duf(counter,1)= 1/dt * (a(counter)/V(counter,nt+1) - b(counter)/theta(counter,nt+1)*(Dc(counter)*dt*(theta(counter,nt)+dt))/((Dc(counter)+dt*V(counter,nt+1))^2));
counter = counter +1;
end
end

% argument = (V(counter,nt+1) ./ (2 * V0)) * exp((uf0 + b * log(V0 * theta(counter,nt+1) / Dc)) ./ a);
% uf(counter,nt+1) = a .* asinh(argument);
% duf(counter,1)=((a/(abs(V(counter,nt+1))*dt))-(alpha*b/((1+alpha*abs(V(counter,nt+1)))*dt)))*tanh(uf(counter,nt+1)/a);
% K=(1/(2*V0))*exp((uf0+b*log(V0*theta(counter,nt+1)/Dc))/a);
% duf (counter,1)= 1/dt * (a*K/sqrt(1+(K*V(counter,nt+1))^2)-(b*K*V(counter,nt+1)/(theta(counter,nt+1)*sqrt(1+(K*V(counter,nt+1))^2)))*(Dc*dt*(theta(counter,nt)+dt)/(Dc+dt*V(counter,nt+1))^2));