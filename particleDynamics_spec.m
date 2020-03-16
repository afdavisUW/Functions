function [U1,U2,U3] = particleDynamics_spec(time)
%%This function calculates the vertical dynamics of a water particle given a
%global spectrum and a specific time value.

% U1 - disp
% U2 - vel
% U3 - acc

global SpectrumData

omega = SpectrumData.W;
Mag = SpectrumData.amp;
Phi = SpectrumData.phase;
sum1 = 0; sum2 = 0; sum3 =0;

% displacement loop
for k = 1:length(omega)
    sum1 = sum1 + Mag(k)*sin(omega(k)*time+Phi(k));
end
U1 = sum1;

% velocity loop
for k = 1:length(omega)
    sum2 = sum2 + omega(k)*Mag(k)*cos(omega(k)*time+Phi(k));
end
U2 = sum2;

% acceleration loop
for k = 1:length(omega)
    sum3 = sum3 - (omega(k))^2*Mag(k)*sin(omega(k)*time+Phi(k));
end
U3 = sum3;


end