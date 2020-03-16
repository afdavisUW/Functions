function [NMSE] = perform_ARvarK(forecastParams,oRingIndex)
% coded 3/28/19 based on perform_AR_dev
% original code is based on my ARX implementation 
% This function is meant to be used to produce a figure that
% non-dimensionalizes

% forecastParams contains the following parameters:
% methodNum, runstart, Duration, Fstart, subNo, Tstart, Tend, oRingIndex, na,

% perform_AR is a real time implementation of the AR function using
% the OWC and the FPS with various different approaches: distrbance,
% estimated force, analytical force determinatiion. to determine the
% performance. 

% AR implementation
predictMethodOptions = ["FexTrad","ConvTrad","FexDist"];
predictMethod = predictMethodOptions(forecastParams.methodNum);
TpOptions = [0,0,1.9,1.9,1.9,2.37,2.37,2.37,2.69,2.69,2.69,3.16,3.16,3.95,4.74,5.53,2.06,2.37]; % Tp in s
Tp = TpOptions(oRingIndex);
forecastTime = Tp/2;
Ts = forecastTime/forecastParams.K;
Fs = 1/Ts;
subNo = round(200/Fs);

% parameters
% Horizon = 1.5 % [s] forecast horizon
simTimeStart = forecastParams.runStart; % time into data set that the AR functions start
simDuration = forecastParams.Duration;
forecastStart = forecastParams.Fstart;
% subNo = forecastParams.subNo;
% Fs = 200/subNo;
% Ts = 1/Fs;
% distance = 7.5; % meters separating FPS and OWC
g = 9.81;
sampleStart = round(simTimeStart*Fs)+1; % MATLAB is a 1 index language
sampleEnd = round(simDuration*Fs)+1;
sampleFstart = round(forecastStart*Fs)+1;
a = sampleStart; b = sampleEnd; c = sampleFstart;
Tstart = forecastParams.Tstart; 
Tend = forecastParams.Tend;
% K = round(Horizon*Fs)

%% Setup Model
% oRingIndex = forecastParams.oRingIndex;
[FPS,OWC,deadTime] = loadFexData(oRingIndex,subNo,Tstart,Tend);
% DistOrder = length(OWC.W); 
% DistOrder = 3

% Fex -> traditional AR
if predictMethod == 'FexTrad' || predictMethod == 'ConvTrad'
    na = forecastParams.na; % number of poles
    modelOrders = [na];
    modelObj = recursiveAR(modelOrders);
    modelObj.ForgettingFactor = 0.995;
    nk = round(deadTime/Ts);
    
    %     selecting the number of forecasting steps
    K = forecastParams.K;
    
elseif predictMethod == 'FexDist'
    
    return
    
end

%% Implement recursive ARX

if predictMethod == 'FexDist'
    return
    
elseif predictMethod == 'ConvTrad'
    for kk = a:b
        [~,~] = modelObj(OWC.conv(kk)); % stepping through the use of a recursive regression function
%         K = nk;
        
        if kk >= c
            system = idpoly(modelObj);
            pastDataInLoop = iddata(OWC.conv(a:kk)',FPS.conv(a:kk)',Ts);
            futureInputsInLoop = iddata([],FPS.conv(kk+1:kk+K)',Ts);
            yforecastInLoop = forecast(system,pastDataInLoop,K,futureInputsInLoop);
            NMSEinLoop(kk) = goodnessOfFit(yforecastInLoop.OutputData,OWC.conv(kk+1:kk+K)','NMSE');
            yforecastK(kk) = yforecastInLoop.OutputData(end);
            percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
            disp(['Percent Complete: ',num2str(percentComplete),'%'])
        end
    
    end

    
elseif predictMethod == 'FexTrad'
    for kk = a:b
        [~,~] = modelObj(OWC.Est(kk)); % stepping through the use of a recursive regression function
%         K = nk;

        %         forcasting for each step
        if kk >= c
            system = idpoly(modelObj);
            pastDataInLoop = iddata(OWC.Est(a:kk)',FPS.Est(a:kk)',Ts);
            futureInputsInLoop = iddata([],FPS.Est(kk+1:kk+K)',Ts);
            yforecastInLoop = forecast(system,pastDataInLoop,K,futureInputsInLoop);
            NMSEinLoop(kk) = goodnessOfFit(yforecastInLoop.OutputData,OWC.Est(kk+1:kk+K)','NMSE');
            yforecastK(kk) = yforecastInLoop.OutputData(end);
            percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
            disp(['Percent Complete: ',num2str(percentComplete),'%'])
        end
    end
    
end

%% Visualizaing Results, Compute NMSE

% visualizations for Fex and Conv
if predictMethod == 'FexTrad' || predictMethod == 'ConvTrad'
    if predictMethod == 'FexTrad'
        Measured = OWC.Est(b+1:b+K)';
    elseif predictMethod == 'ConvTrad'
        Measured = OWC.conv(b+1:b+K)';
    end
    
    % comparing K step ahead predictions
    forecastIstart = round(forecastStart/Ts)+1; % round introduced 4/21/19 to allow non roots of 200 to be useed as sbuNo
    forecastIstop = round(simDuration/Ts)+1;
    hold on
    plot(OWC.Time(forecastIstart+K:forecastIstop+K),yforecastK(forecastIstart:forecastIstop),'b')
    plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.conv(forecastIstart+K:forecastIstop+K),'r')
    plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.Est(forecastIstart+K:forecastIstop+K),'m')
    if predictMethod == 'FexTrad'
        NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.Est(forecastIstart+K:forecastIstop+K)','NMSE');
    elseif predictMethod == 'ConvTrad'
        NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.conv(forecastIstart+K:forecastIstop+K)','NMSE');
    end
    legend('Fpredict','Fconv','Fex')
    title(['NMSE = ',num2str(NMSE),', Prediction: K = ',num2str(K),'(',num2str(K*Ts),'s)'])
    
elseif predictMethod == 'FexDist'
        
%         % comparing K step ahead predictions
%         forecastIstart = round(forecastStart/Ts)+1; % round introduced 4/21/19 to allow non roots of 200 to be useed as sbuNo
%         forecastIstop = round(simDuration/Ts)+1;
%         figure()
%         hold on
%         plot(OWC.Time(forecastIstart+K:forecastIstop+K),yforecastK(forecastIstart:forecastIstop),'b')
%         plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.Est(forecastIstart+K:forecastIstop+K),'m')
%         NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.Est(forecastIstart+K:forecastIstop+K)','NMSE');
%         legend('Fpredict','Fex')
%         title(['NMSE = ',num2str(NMSE),', Prediction: K = ',num2str(K),'(',num2str(K*Ts),'s)'])
    
end


end


