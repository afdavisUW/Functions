function [measurementOutput] = check_measurements(measString,dataName,Paths,timeSpan)
%check_measurements: 1/16/18 is meant to be the decision making function
%that either loads experimental data or loads simulated data. The original
%implementation will be to load the 10th order state space realization data
%that was used as the reference in the first journal paper.

global Func Data Param Time

if strcmp(measString,'loadSim')
    load(dataName)
elseif strcmp(measString,'sim')
    simTime = [timeSpan(1):1/Param.samplingFreq.motion:timeSpan(2)];
    [~,simTrueStates] = ode45(Func.fODE, simTime, Param.IC.true);
    simTrueMeas = simTrueStates(:,1:2);
        
%     adding the state noise
%     simStates = simTrueMeas;
%     statesNoise = sqrt(Param.Q)*randn(length(simTime),Param.numOutputs);
%     simStates(:,1) = simTrueMeas(:,1)+Param.G(1)*statesNoise(:,1);
%     simStates(:,2) = simTrueMeas(:,2)+Param.G(2)*statesNoise(:,2);
    measNoise = [Param.Ra(1)*randn(length(Time.t),1),Param.Ra(2)*randn(length(Time.t),1)];
    stateNoise = [Param.Ga(1)*randn(length(Time.t),1),Param.Ga(2)*randn(length(Time.t),1)];
    Measurements = simTrueMeas + measNoise+stateNoise;

%     simStates = simTrueMeas;
    
%     adding measurement noise
%     measNoise = sqrt(Param.R(1,1))*randn(length(simTime),Param.numOutputs);
%     Measurements = simStates + measNoise;
    
elseif strcmp(measString,'loadExper')
    if Param.bodyNumber == 1
        [minimum,startIndex] = min(abs(Data.motion.time-Time.start));
        Measurements = [Data.motion.FPS.heave(startIndex:Time.subSampleNo:startIndex+Time.subSampleNo*length(Time.t)-1),...
            Data.motion.FPS.heaveVel(startIndex:Time.subSampleNo:startIndex+Time.subSampleNo*length(Time.t)-1)];
        Param.residual = Data.motion.FPS.residual(startIndex:Time.subSampleNo:startIndex+Time.subSampleNo*length(Time.t)-1);
    elseif Param.bodyNumber == 2
        [minimum,startIndex] = min(abs(Data.motion.time-Time.start));
        Measurements = [Data.motion.OWC.heave(startIndex:Time.subSampleNo:startIndex+Time.subSampleNo*length(Time.t)-1),...
            Data.motion.OWC.heaveVel(startIndex:Time.subSampleNo:startIndex+Time.subSampleNo*length(Time.t)-1)];
        Param.residual = Data.motion.OWC.residual(startIndex:Time.subSampleNo:startIndex+Time.subSampleNo*length(Time.t)-1);
    end        
        
%     artificial noise addition
    measNoise = [Param.Ra(1)*randn(length(Time.t),1),Param.Ra(2)*randn(length(Time.t),1)];
    stateNoise = [Param.Ga(1)*randn(length(Time.t),1),Param.Ga(2)*randn(length(Time.t),1)];
    Measurements = Measurements + measNoise + stateNoise;

else
    display('no data')
    return
end

measurementOutput = Measurements;
end

