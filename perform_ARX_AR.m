function [NMSE] = perform_ARX_AR(forecastParams,oRingIndex)
% perfom_ARX coded 3/28/19 based on perform_ARX_dev
% forecastParams contains the following parameters:
% runstart, Duration, Fstart, subNo, Tstart, Tend, oRingIndex, methodNum,
% na, nb

% perform_ARMAX is a real time implementation of the ARMAX function using
% the OWC and the FPS with various different approaches: distrbance,
% estimated force, analytical force determinatiion. to determine the
% performance. 

% format compact
% addpath(fullfile(cd, '..', filesep, 'Functions'))

% ARX implementation
predictMethodOptions = ["FexTrad","ConvTrad","FexDist"];
predictMethod = predictMethodOptions(forecastParams.methodNum);

% parameters
% Horizon = 1.5 % [s] forecast horizon
% simTimeStart = 0; % time into data set that the AR functions start
simTimeStart = forecastParams.runStart;
% simDuration = 75;
simDuration = forecastParams.Duration;
% forecastStart = 65;
forecastStart = forecastParams.Fstart;
% subNo = 50;
subNo = forecastParams.subNo;
Fs = 200/subNo;
Ts = 1/Fs;
distance = 7.5; % meters separating FPS and OWC
g = 9.81;
sampleStart = round(simTimeStart*Fs)+1; % MATLAB is a 1 index language
sampleEnd = round(simDuration*Fs)+1;
sampleFstart = round(forecastStart*Fs)+1;
a = sampleStart; b = sampleEnd; c = sampleFstart;
% Tstart = 200; % for loading saved data
Tstart = forecastParams.Tstart; 
% Tend = 400; % for loading saved data
Tend = forecastParams.Tend;
% K = round(Horizon*Fs)

%% Setup Model
% oRingIndex = 6;
% oRingIndex = ;
[FPS,OWC,deadTime] = loadFexData(oRingIndex,subNo,Tstart,Tend);
forgetFactor = 0.99;
% 0.995 also common

% Fex -> traditional ARX
if predictMethod == 'FexTrad' || predictMethod == 'ConvTrad'

%     setup ARX model
    na = forecastParams.na;
    nk = round(deadTime*Fs); % number of dead samples between FPS and OWC
    nb = forecastParams.nb;
    modelOrdersARX = [na nb nk];
    modelObjARX = recursiveARX(modelOrdersARX);
    modelObjARX.ForgettingFactor = forgetFactor;
    
%     setup AR model
    modelOrdersAR = [na];
    modelObjAR = recursiveAR(modelOrdersAR);
    modelObjAR.ForgettingFactor = forgetFactor;
    
    if forecastParams.forecast == 0
        K = nk;
    else
        K = round(forecastParams.forecast/Ts);
    end

elseif predictMethod == 'FexDist'
    return
 
end

%% Implement recursive ARX

if predictMethod == 'FexDist'
    return

elseif predictMethod == 'ConvTrad'
    for kk = a:b
        [~,~,~] = modelObjARX(OWC.conv(kk),FPS.conv(kk)); % stepping through the use of a recursive regression function
        [~,~] = modelObjAR(OWC.conv(kk));

        if kk >= c
            if K <= nk % ARX prediction only
                
                system = idpoly(modelObjARX);
                pastDataInLoop = iddata(OWC.conv(a:kk)',FPS.conv(a:kk)',Ts);
                futureInputsInLoop = iddata([],FPS.conv(kk+1:kk+K)',Ts);
                yforecastInLoop = forecast(system,pastDataInLoop,K,futureInputsInLoop);
                NMSEinLoop(kk) = goodnessOfFit(yforecastInLoop.OutputData,OWC.conv(kk+1:kk+K)','NMSE');
                yforecastK(kk) = yforecastInLoop.OutputData(end);
                percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
                disp(['Percent Complete: ',num2str(percentComplete),'%'])
            elseif K > nk
                
                system1 = idpoly(modelObjARX);
                pastDataInLoop = iddata(OWC.conv(a:kk)',FPS.conv(a:kk)',Ts);
                futureInputsInLoop = iddata([],FPS.conv(kk+1:kk+nk)',Ts);
                yforecastARX = forecast(system1,pastDataInLoop,nk,futureInputsInLoop);
                
                %             update AR with ARX predictions
%                 for j = 1:length(yforecastARX.OutputData)
%                     [~,~] = modelObjAR(yforecastARX.OutputData(j));
%                 end
                
                system2 = idpoly(modelObjAR);
                pastData = iddata([OWC.conv(a:kk),yforecastARX.OutputData']',[],Ts);
                yforecastAR = forecast(system2,pastData,K-nk);
                yforecastK(kk) = yforecastAR.OutputData(end);
                
                percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
                disp(['Percent Complete: ',num2str(percentComplete),'%'])
                
            end
            
        end

    end
elseif predictMethod == 'FexTrad'
    for kk = a:b
        [~,~,~] = modelObjARX(OWC.Est(kk),FPS.Est(kk)); % stepping through the use of a recursive regression function
        [~,~] = modelObjAR(OWC.Est(kk));

        %         forcasting for each step
        if kk >= c
            if K <= nk % ARX prediction only
                system = idpoly(modelObjARX);
                pastDataInLoop = iddata(OWC.Est(a:kk)',FPS.Est(a:kk)',Ts);
                futureInputsInLoop = iddata([],FPS.Est(kk+1:kk+K)',Ts);
                yforecastInLoop = forecast(system,pastDataInLoop,K,futureInputsInLoop);
                NMSEinLoop(kk) = goodnessOfFit(yforecastInLoop.OutputData,OWC.Est(kk+1:kk+K)','NMSE');
                yforecastK(kk) = yforecastInLoop.OutputData(end);
                percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
                disp(['Percent Complete: ',num2str(percentComplete),'%'])
            elseif K > nk
                
                system1 = idpoly(modelObjARX);
                pastDataInLoop = iddata(OWC.Est(a:kk)',FPS.Est(a:kk)',Ts);
                futureInputsInLoop = iddata([],FPS.Est(kk+1:kk+nk)',Ts);
                yforecastARX = forecast(system1,pastDataInLoop,nk,futureInputsInLoop);
                
%                 %             update AR with ARX predictions
%                 for j = 1:length(yforecastARX.OutputData)
%                     [~,~] = modelObjAR(yforecastARX.OutputData(j));
%                 end
                
                system2 = idpoly(modelObjAR);
                pastData = iddata([OWC.Est(a:kk),yforecastARX.OutputData']',[],Ts);
                yforecastAR = forecast(system2,pastData,K-nk);
                yforecastK(kk) = yforecastAR.OutputData(end);
                
                percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
                disp(['Percent Complete: ',num2str(percentComplete),'%'])
                
            end
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
    NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.Est(forecastIstart+K:forecastIstop+K)','NMSE');
    if forecastParams.plotBool == 1
    figure()
    hold on
    plot(OWC.Time(forecastIstart+K:forecastIstop+K),yforecastK(forecastIstart:forecastIstop),'b')
    plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.conv(forecastIstart+K:forecastIstop+K),'r')
    plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.Est(forecastIstart+K:forecastIstop+K),'m')
    legend('Fpredict','Fconv','Fex')
    title(['NMSE = ',num2str(NMSE),', Prediction: K = ',num2str(K),'(',num2str(K*Ts),'s)'])
    end
elseif predictMethod == 'FexDist'
        
        % comparing K step ahead predictions
        forecastIstart = round(forecastStart/Ts)+1; % round introduced 4/21/19 to allow non roots of 200 to be useed as sbuNo
        forecastIstop = round(simDuration/Ts)+1;
        NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.Est(forecastIstart+K:forecastIstop+K)','NMSE');
        if forecastParams.plotBool == 1
        figure()
        hold on
        plot(OWC.Time(forecastIstart+K:forecastIstop+K),yforecastK(forecastIstart:forecastIstop),'b')
        plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.Est(forecastIstart+K:forecastIstop+K),'m')
        legend('Fpredict','Fex')
        title(['NMSE = ',num2str(NMSE),', Prediction: K = ',num2str(K),'(',num2str(K*Ts),'s)'])
        end

    
    
end

