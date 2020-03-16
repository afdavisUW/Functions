function [U1,U2,U3] = estimateVel(tau)
%estimateVel is meant to provide an estimate of the water particle velocity
%based on the estimated wave excitation force. 2/28/2018 
% U1 - water particle position
% U2 - water particle velocity
% U3 - water particle acceleration

global Param Est smoothCounter Time
timeIndex = smoothCounter;

if Param.velEstimateVal == 0
    [U1,U2,U3] = particleDynamics_spec(tau);
    
elseif timeIndex == 1
    U2 = 0;
    
elseif Param.velEstimateVal == 1 % using method 1
    Fex2 = Est.states(3,timeIndex);
    Fex1 = Est.states(3,timeIndex-1);
    U2 = ((Fex2/Param.FexLumped)-(Fex1/Param.FexLumped))/Time.dT;
    Est.vel(timeIndex) = U2;
    
elseif Param.velEstimateVal == 2 
    
%     Param.numFexComponents
n = Param.numFexComponents;
Fex2 = sum(Est.states(3:2:2*(n-1)+3,timeIndex));
Fex1 = sum(Est.states(3:2:2*(n-1)+3,timeIndex-1));
U2 = ((Fex2/Param.FexLumped)-(Fex1/Param.FexLumped))/Time.dT;
Est.vel(timeIndex) = U2;

end
     


U1 = [];
U3 = [];

end

