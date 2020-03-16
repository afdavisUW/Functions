function [F_ex ] = excitationConv(time,F_exIRF,IRFtime,Disp,DispTime)
%8/18/17 This function is meant to use the convolution approach to
%determnine the exact wave excitaiton force that is a non-causal funtion

% checking time domain boundaries
boundary = 5; % this specifies the minimum time in order to guarentee that 
% there is data to actually compute the convoltuion

if time <= boundary
    display('Specified time too early!!!')
    display('This results in NaN in timeDisp')
    retrun
end

% global IRFDisp IRFTime WaveDisp WaveTime
%     F_exIRF = IRFDisp;
%     IRFtime = IRFTime;
%     Disp = WaveDisp;
%     DispTime = WaveTime;
    
startIRF = -boundary;
endIRF = boundary;
startDisp = time-(-boundary);
endDisp = time-boundary;
dt = 0.005;
timeIRF = [startIRF:dt:endIRF];
timeDisp = [startDisp:-dt:endDisp];

termIRF = interp1(IRFtime,F_exIRF,timeIRF);
termDisp = interp1(DispTime,Disp,timeDisp);


Prod = termIRF.*termDisp;        


F_ex = dt*trapz(Prod);


end

