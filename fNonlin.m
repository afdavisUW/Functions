function [ Zdot ] = fNonlin( tau, Z )

global Func 
% U1 - displacement
% U2 - velocity
% U3 - Acceleration

% [U1,U2,U3] = particleDynamics_spec(tau);
[U1,U2,U3] = estimateVel(tau);

Zdot = Func.f(Z',U1,U2,U3,tau);

end
