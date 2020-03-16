function [] = generateBody_BosmaFPS(bodyName, Paths, desiredOrder )
%Function meant to generate all mesh, freqency, SS, and IRF data
%   8/29/17 This function is meant to be the file that houses all aspects
%   concerning the data of a body in the water. This function will excecute
%   aximesh.m nemoh.m generate a SS realization and an IRF representation.
% TO DO: choose better saving mechanisms that reduce the amount of
% unaccessed files and make the integration limits more robust. As well as
% use trapz on the integrations

% checking the operating system
% as of 9/14/17 only windows has been shown to work with NEMOH
if isunix
    display('NEMOH and aximesh will not work on Linux')
    return
elseif ismac
    display('It is unknown if NEMOH and aximesh will work on Mac, run at own risk')
end
%% Generate Mesh for NEMOH
if strcmp(bodyName,'FPS') == 1
    Rb = (.2706+0.25)/2; % intermediate big radius for FPS
    Rs = (.1894+.175)/2; % small radius for FPS
    z = [0.15,0.15,0,-0.07,-0.075,-0.075];
    r = [0,Rb,Rb,Rs+(Rb-Rs)/3,Rs,0];
    display('FPS: Enter the Following: 36, storage, 0.05, 500')
    [Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,length(z));
elseif strcmp(bodyName,'OWC') == 1
%     Rb  used asd well.
    R = (0.4059+0.375)/2;
    Rb = (.2706+0.25)/2; % intermediate small radius for OWC
    tt = 0.01; %10mm thick polycrbonate
%     z = [0.75,0.75,0.50,0.44,0.06,0,-1.0,-1.0,0.75];
    z = [0.25,0.22,-0.22,-0.25];
%     r = [Rb-tt,Rb,Rb,R,R,Rb,Rb,Rb-tt,Rb-tt];
    r = [Rb,R,R,Rb];
    display('OWC: Enter the Following: 36, storage, 0.2, 500')
    [Mass,Inertia,KH,XB,YB,ZB] = axiMesh(r,z,length(z));
else
    display('bodyName does not select aximesh points')
    return
end


%% Excecute NEMOH

num_omega = 25;
max_period =10;
min_omega = 2*pi/max_period; 
min_period = .15;
max_omega = 2*pi/min_period;
% rounding omega to the nearest radian/second
min_omega = round(min_omega,1); max_omega = round(max_omega,1); % Warning!!!
% watch these round functions, because the most common error I am getting
% is when the round function takes one of the input frequencies to be zero
omega = linspace(min_omega,max_omega,num_omega);

num_angle = 1;
min_angle = 0;
max_angle = 0;
angle_specified = linspace(min_angle,max_angle,num_angle);
depth = 0; % deep water waves
[A,B,Fe] = Nemoh(omega,angle_specified,depth);

save(strcat(bodyName,'_ABFe_freq.mat'),'A','B','Fe','omega')

plotBool = [0]; % 1 if you want plots 0 if you don't
if plotBool == 1
%   Plotting results
    figure()
    B_BEM(1,:) = B(3,3,:);
    plot(omega,B_BEM(1,:))
    title('Radiation Damping')
    figure()
    plot(omega,squeeze(A(3,3,:)))
    title('Added Mass')
    figure()
    plot(omega,real(Fe(:,3)))
    title('real part of excitation force')
    figure()
    plot(omega,imag(Fe(:,3)))
    title('imag part of excitation force')
    figure()
    plot(omega,abs(Fe(:,3)))
    title('magnitude of Fex')
    figure()
    plot(omega,angle(Fe(:,3)))
    title('phase of Fex')
end

clc
display('Nemoh and mesh print statements cleared')
%% State-Space Realization and Impulse Response Function
NumStates = desiredOrder; % Number of states for the reduced SS realization

% Creating SS realization of the Radiation Damping Coefficient
b = squeeze(B(3,3,:));
% frequency depended excitation force magnitude
fe = [fliplr(conj(Fe(:,3))'),Fe(:,3)']; % this needs to be a 2 sided spectrum

figure()
plot(omega,b,'bo')
hold on
dW = 0.01;
xx = 0:dW:max(omega); % desired Resolution
yy = spline(omega,b,xx); % creating cubic spline data to fill in gaps
plot(xx,yy)
hold off

% compute the impulse response function using trapazoidal integration
h = 0.002;
t = [0:h:5];

Sum = zeros(1,length(t));
for jj = 1:length(t)
    for j = 1:length(xx)-1
        
        Trap = 1/2*(yy(j)+yy(j)) * dW;
        Trap = Trap * cos(xx(j)*t(jj));
        Sum(jj) = Sum(jj)+Trap;
    end    
end
IRF = 2/pi*Sum;
% figure
% plot(t,IRF)

% using imp2ss to determine the impulse response function
[G] = imp2ss(IRF*h,h); % propper scaling based on Nathan Tom's
% 
% global A_r B_r C_r D_r
A_r = G.a;
B_r = G.b;
C_r = G.c;
D_r = G.d;
% save('ABCD_r_SSreal.mat','A_r','B_r','C_r','D_r')

% analyzing the stability of A_r
size(A_r)
lambda = eig(A_r);
MAX = max(real(lambda));
AVG = mean(real(lambda));

% looking at the comparison between the original frequency response
H = freqresp(G,xx);
for j = 1:length(xx)
    Mag(j) = abs(real(H(:,:,j)));
end

figure()
plot(xx,yy,'b',xx,Mag,'r')
legend('WAMIT','realization')
title('F.R.F.')

% looking at the impulse response of the system
[Y,time] = impulse(G,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('I.R.F. Radiation')


%% Reducing the order

MRTYPE = 1; % stands for model reduction type
[AM,BM,CM,DM,TOTBND,HSV] = balmr(A_r,B_r,C_r,D_r,MRTYPE,NumStates);

% analyzing the stability of A_r
size(AM);
lambda = eig(AM);
MAX = max(real(lambda));
AVG = mean(real(lambda));

clear H Mag
reduced_model = ss(AM,BM,CM,DM);
% looking at the comparison between the original frequency response
H = freqresp(reduced_model,xx);
for j = 1:length(xx)
    Mag(j) = abs(real(H(1,1,j)));
end

figure()
plot(xx,yy,'b',xx,Mag,'r.')
legend('Nemoh','realization')
legend('Nemoh','realization')
title('F.R.F. Reduced System')

% looking at the impulse response of the system
[Y,time] = impulse(reduced_model,t);
figure()
plot(time,Y,'r.',t,IRF,'b')
title('IRF of the reduced system')

clear A_r B_r C_r D_r
% global A_r B_r C_r D_r
A_r = AM;
B_r = BM;
C_r = CM;
D_r = DM;
% save('ABCD_r_SSreal.mat','A_r','B_r','C_r','D_r')
save(strcat(bodyName,'_ABCD_ssReal.mat'),'A_r','B_r','C_r','D_r')


%% Excitation Force

% xx = -31:0.1:31; % desired Resolution (frequency)
dw = 0.001;
xx = [-max(omega):dw:max(omega)];
omega_2sided = [-fliplr(omega),omega];
yyR = spline(omega_2sided,real(fe),xx); % creating cubic spline data to fill in gaps
yyI = spline(omega_2sided,imag(fe),xx); %  

% computing the impulse response function of the excitation force
% dw = xx(2)-xx(1);
Time = -15:0.001:15;

% computing the values for integration
for k = 1:length(Time)
intArg1 = yyR.*cos(xx.*Time(k)) - yyI.*sin(xx.*Time(k)); 
intArg2 = yyR.*sin(xx.*Time(k)) + yyI.*cos(xx.*Time(k));

f_t(k) = dw/(2*pi)*(trapz(intArg1) + 1i*trapz(intArg2));

end

F_ex.IRF = real(f_t); % f_t is nearly purely real so applying the real funciton is appropriate
F_ex.time = Time;
F_ex.W = omega;

figure()
plot(Time,real(f_t))
title('impulse response function (Fex)')

save(strcat(bodyName,'FexIRF.mat'),'F_ex')

val = exist('storage','dir');
if val == 0
    display('NEMOH file "storage" not found! Did you name the NEMOH folder correctly?')
end


close all

movefile('storage',Paths.outputDir);
movefile('ID.dat',Paths.outputDir);
movefile('mesh.cal',Paths.outputDir);
movefile(strcat(bodyName,'_ABFe_freq.mat'),Paths.outputDir);
movefile(strcat(bodyName,'_ABCD_ssReal.mat'),Paths.outputDir);
movefile(strcat(bodyName,'FexIRF.mat'),Paths.outputDir);


end

