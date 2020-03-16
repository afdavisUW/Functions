function [NMSE] = perform_ARMA(forecastParams,oRingIndex)
% coded 3/28/19 based on perform_AR_dev
% original code is based on my ARX implementation 

% forecastParams contains the following parameters:
% methodNum, runstart, Duration, Fstart, subNo, Tstart, Tend, oRingIndex, na,

% perform_AR is a real time implementation of the AR function using
% the OWC and the FPS with various different approaches: distrbance,
% estimated force, analytical force determinatiion. to determine the
% performance. 

% AR implementation
predictMethodOptions = ["FexTrad","ConvTrad","FexDist"];
predictMethod = predictMethodOptions(forecastParams.methodNum);

% parameters
% Horizon = 1.5 % [s] forecast horizon
simTimeStart = forecastParams.runStart; % time into data set that the AR functions start
simDuration = forecastParams.Duration;
forecastStart = forecastParams.Fstart;
subNo = forecastParams.subNo;
Fs = 200/subNo;
Ts = 1/Fs;
% distance = 7.5; % meters separating FPS and OWC
g = 9.81;
sampleStart = round(simTimeStart*Fs)+1; % MATLAB is a 1 index language
sampleEnd = round(simDuration*Fs)+1;
a = sampleStart; b = sampleEnd;
Tstart = forecastParams.Tstart; 
Tend = forecastParams.Tend;
% K = round(Horizon*Fs)

%% Setup Model
% oRingIndex = forecastParams.oRingIndex;
[FPS,OWC,deadTime] = loadFexData(oRingIndex,subNo,Tstart,Tend);
DistOrder = length(OWC.W); 
% DistOrder = 3

% Fex -> traditional AR
if predictMethod == 'FexTrad' || predictMethod == 'ConvTrad'
    na = forecastParams.na; % number of poles
    nc = forecastParams.nc;
    modelOrders = [na nc];
    modelObj = recursiveARMA(modelOrders);
    modelObj.ForgettingFactor = 0.99;
    nk = round(deadTime/Ts);

elseif predictMethod == 'FexDist'
        
%     na = 15; % number of poles
% %     tempVec = [35,40,70]
%     %     setting up an ARX with information from each frequency
%     for kk = 1:DistOrder
%         deadTime(kk) = distance/(g/FPS.W(kk)); % distance divided by celerity
% %         Kstep(kk) = round(deadTime(kk)*Fs);
% %         nb(1,kk) = 55; % number of zeros plus one
%         nk(1,kk) = round(deadTime(kk)*Fs); % number of dead samples between FPS and OWC
%         nb(1,kk) = 2*nk(1,kk); % number of zeros plus one
% 
%     end
%     modelOrders = [na nb nk]
%     modelObj = recursiveARX(modelOrders);
%     modelObj.ForgettingFactor = 0.995
  
end

%% Implement recursive ARX

if predictMethod == 'FexDist'
    
    for kk = a:b
%         [A,B,EstOutput] = modelObj(OWC.Est(kk),FPS.components(1:DistOrder,kk)); % stepping through the use of a recursive regression function
%         K = nk(1);
% 
%         if (kk-a) >= forecastStart/Ts
%             system = idpoly(modelObj);
%             pastDataInLoop = iddata(OWC.Est(a:kk)',FPS.components(1:DistOrder,a:kk)',Ts);
%             futureInputsInLoop = iddata([],FPS.components(1:DistOrder,kk+1:kk+nk(1))',Ts);
% %             note, using nk(1) as the number of steps as it is the limiting factor
%             yforecastInLoop = forecast(system,pastDataInLoop,K,futureInputsInLoop);
% %             NMSEinLoop(1,kk) = goodnessOfFit(yforecastInLoop.OutputData,OWC.Est(kk+1:kk+K)','NMSE');
%             yforecastK(1,kk) = yforecastInLoop.OutputData(end);
%             percentComplete = round(100*(kk-forecastStart/Ts)/(b-forecastStart/Ts));
%             disp(['Percent Complete: ',num2str(percentComplete),'%'])
%         end
%         
    end


elseif predictMethod == 'ConvTrad'
    for kk = a:b
        [~,~] = modelObj(OWC.conv(kk)); % stepping through the use of a recursive regression function
        K = nk;
        
        if (kk-a) >= forecastStart/Ts
            system = idpoly(modelObj);
            pastDataInLoop = iddata(OWC.conv(a:kk)',FPS.conv(a:kk)',Ts);
            futureInputsInLoop = iddata([],FPS.conv(kk+1:kk+nk)',Ts);
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
        K = nk;

        %         forcasting for each step
        if (kk-a) >= forecastStart/Ts
            system = idpoly(modelObj);
            pastDataInLoop = iddata(OWC.Est(a:kk)',FPS.Est(a:kk)',Ts);
            futureInputsInLoop = iddata([],FPS.Est(kk+1:kk+nk)',Ts);
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
    NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.Est(forecastIstart+K:forecastIstop+K)','NMSE');
    legend('Fpredict','Fconv','Fex')
    title(['NMSE = ',num2str(NMSE),', Prediction: K = ',num2str(K),'(',num2str(K*Ts),'s)'])
    
elseif predictMethod == 'FexDist'
        
        % comparing K step ahead predictions
        forecastIstart = round(forecastStart/Ts)+1; % round introduced 4/21/19 to allow non roots of 200 to be useed as sbuNo
        forecastIstop = round(simDuration/Ts)+1;
        figure()
        hold on
        plot(OWC.Time(forecastIstart+K:forecastIstop+K),yforecastK(forecastIstart:forecastIstop),'b')
        plot(OWC.Time(forecastIstart+K:forecastIstop+K),OWC.Est(forecastIstart+K:forecastIstop+K),'m')
        NMSE = goodnessOfFit(yforecastK(forecastIstart:forecastIstop)',OWC.Est(forecastIstart+K:forecastIstop+K)','NMSE');
        legend('Fpredict','Fex')
        title(['NMSE = ',num2str(NMSE),', Prediction: K = ',num2str(K),'(',num2str(K*Ts),'s)'])
    
end


end


