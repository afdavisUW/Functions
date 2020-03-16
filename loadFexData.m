function [FPS,OWC,deadTime] = loadFexData(oRingIndex,subNo,Tstart,Tend)
%This function is inteded to clean up some of the AutoRegressive scripts by
%allowing a single funtion to recal the important information so that it
%can be used as needed. This may disrupt the capabilities of parallel
%computing, however loading the files in the script would do the same.

Paths.dataDir = fullfile(cd, '..', filesep, 'Data', filesep, 'FexData');
addpath(Paths.dataDir)

waveRuns = [1,180,0.075,1.1;
            2,181,0.175,1.42;
            3,182,0.05,1.9;
            4,183,0.175,1.9;
            5,184,0.3,1.9;
            6,185,0.050,2.37;
            7,186,0.175,2.37;
            8,199,0.3,2.37;
            9,189,0.05,2.69;
            10,190,0.175,2.69;
            11,191,0.3,2.69;
            12,192,0.05,3.16;
            13,193,0.175,3.16;
            14,194,0.175,3.95;
            15,195,0.175,4.74;
            16,196,0.175,5.53;
            17,197,0.375,2.06;
            18,198,0.375,2.37];
runNo = num2str(waveRuns(oRingIndex,2));
wavePeriod = waveRuns(oRingIndex,4);
waveHs = waveRuns(oRingIndex,3);
distance_FPStoOWC = 7.5; % [m]
g = 9.81; % [m/s^2]
Celerity = g*wavePeriod/(2*pi); % [m/s]
deadTime = distance_FPStoOWC/Celerity; % time for wave to travel between FPS and OWC

% Loading estimated excitation data
Temp1 = open(string(strcat('Fex_FPS_Run',num2str(oRingIndex),'_Time',num2str(Tstart),'_',num2str(Tend),'.mat')));
FPS.Est = Temp1.FPS.Est(1:subNo:end);
FPS.Time = Temp1.FPS.Time(1:subNo:end);
FPS.conv = Temp1.FPS.Conv(1:subNo:end); 
FPS.start = Temp1.FPS.start;
FPS.end = Temp1.FPS.end; 
FPS.MSE = Temp1.FPS.MSE; 
FPS.components = Temp1.FPS.components(:,1:subNo:end); 
FPS.W = Temp1.FPS.W;

Temp2 = open(string(strcat('Fex_OWC_Run',num2str(oRingIndex),'_Time',num2str(Tstart),'_',num2str(Tend),'.mat')));
OWC.Est = Temp2.OWC.Est(1:subNo:end); 
OWC.Time = Temp2.OWC.Time(1:subNo:end); 
OWC.conv = Temp2.OWC.Conv(1:subNo:end); 
OWC.start = Temp2.OWC.start; 
OWC.end = Temp2.OWC.end; 
OWC.MSE = Temp2.OWC.MSE; 
OWC.components = Temp2.OWC.components(:,1:subNo:end);
OWC.W = Temp2.OWC.W;
% note the FPS is the input, u, and OWC is the output, y.


end

