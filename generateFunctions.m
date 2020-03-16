function [ func ] = generateFunctions(bodyNumber,vars,dynamics,output,Paths)
% generate functions populates the working foulder with functions
%   This should be 1 of 2 places that the function needs to be specified.
%   The other location is the nonlinear function that is used in the ode
%   solver.
global Param

% bodyNumber is code in as a placeholder for when we need multiple bodies
% simulated.

%% Define Input and Output Functions
% syms symT x1 x2 x3 x4 x5 x6 x7 U1 U2 U3 time
syms U1 U2 U3 symTime
% vars = [x1,x2,x3,x4,x5,x6,x7];

% extracting symbolic variables
% numNonRad = length(vars)-Param.sizeRadSS;
% for k = 1:length(vars)-numNonRad
%     symVec(k) = vars(k+numNonRad);
% end
%     

% Input and Output function needs to be written here
% dynamics = [x2;
%     1/Param.totalInertia*(Param.exciteTerm*U1 - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*x1 - Param.dragFactor*x3*(x2-U2)*abs(x2-U2) - (Param.C_r*[x4;x5;x6;x7]+Param.D_r*x2));
%     0;
%     Param.A_r*[x4;x5;x6;x7]+Param.B_r*x2];

% dynamics = [vars(2);
%     1/Param.totalInertia*(excitationForce_spec(time) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - Param.dragFactor*vars(3)*(vars(2)-U2)*abs(vars(2)-U2) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
%     0;
%     Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];

% testing without relative velocity
% dynamics = [vars(2);
%     1/Param.totalInertia*(excitationForce_spec(time) - (Param.heaveMoorStiff+Param.hydrostaticCoeff)*vars(1) - Param.dragFactor*vars(3)*vars(2) - (Param.body(bodyNumber).C_r*symVec'+Param.body(bodyNumber).D_r*vars(2)));
%     0;
%     Param.body(bodyNumber).A_r*symVec'+Param.body(bodyNumber).B_r*vars(2)];
 

% output = [vars(1);vars(2)];

%% Generate .m files for: f, F, h, and H
Jacob_f = jacobian(dynamics,vars);
Jacob_h = jacobian(output,vars);

tempDir = cd;
cd(Paths.outputDir)
func.f = matlabFunction(dynamics,'File','f','Vars',{vars,U1,U2,U3,symTime});
func.h = matlabFunction(output,'File','h','Vars',{vars});

func.F = matlabFunction(Jacob_f,'File','Jacob_f','Vars',{vars,U1,U2,U3,symTime});
func.H = matlabFunction(Jacob_h,'File','Jacob_h','Vars',{vars});
cd(tempDir)

end
