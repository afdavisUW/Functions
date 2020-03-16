function [DistVals] = chooseDistFreq(bodyNumber,pwrThreshold)
% choosingDistFreq is a script that will use the bretschnider spectrum to
% choose the frequencies, omega,  and the corresponding wave amplitudes 
% that will be used when estimating the excitation force as the disturbance to the EKF

global Param SpectrumData
numComponents = Param.numFexComponents; % this will be the number of spectral divisions
[peakVal,peakIndex] = max(SpectrumData.pwr);
[val,firstThresIndex] = min(abs(SpectrumData.pwr(1:peakIndex) - pwrThreshold*peakVal));
[val,secondThreshIndex] = min(abs(SpectrumData.pwr(peakIndex:end) - pwrThreshold*peakVal));

% generating a new spectrum with only 'high' energy frequencies
num = 10000;
minOmega = SpectrumData.W(firstThresIndex);
maxOmega = SpectrumData.W(secondThreshIndex);
frequency = linspace(minOmega/(2*pi),maxOmega/(2*pi),num);
dW = (frequency(2)-frequency(1))*2*pi;

for j = 1:length(frequency) % loop to calculate specta
    A = frequency(j)/(Param.wavePeriod^(-1));
    S(j) = 5*(Param.waveHs)^2/(16*(Param.wavePeriod)^(-1))*(1/A^5)*exp(-5/4*A^(-4));
end

PowerSpec = S;
% Amplitude = sqrt(2*S.*dW); not applicable selecting few frequencies
Omega = frequency.*2*pi;

% using SpectrumData to determine total energy
totalArea = trapz(PowerSpec); % if evenly spaced, only the total area is necessary
componentArea = totalArea/numComponents;

% vec1 = SpecrtrumData.amp-maxAmp/numComponents;
% PeriodStart = periodRange(1);
% PeriodEnd = periodRange(2);

if numComponents == 1
    clear frequency Amplitude
    DistVals.frequency = Param.wavePeriod^(-1);
    DistVals.amplitude = Param.waveHs;   
    DistVals.force = DistVals.amplitude*Param.FexLumped;
else
%     
    divisionIndex = 1; % starting at the first division of the power
    divisionStart = 1; % staring the with the first index of spectrum Data
    for k = 1:length(PowerSpec)
        if trapz(PowerSpec(1:k)) >= divisionIndex*componentArea
            avgDivPwr = mean(PowerSpec(divisionStart:k));
            [val,index] = min(abs(PowerSpec(divisionStart:k)-avgDivPwr));
            spectrumIndices(divisionIndex) = index+divisionStart-1;
            divisionDf(divisionIndex) = (frequency(k)-frequency(divisionStart));
            
%             increments
            divisionStart = k+1;
            divisionIndex = divisionIndex + 1;
        end
        
    end
    
%     if statement in case there is numerical error and a last frequency is
%     not chosen
    if length(spectrumIndices) < numComponents
        avgDivPwr = mean(PowerSpec(divisionStart:k));
       [val,index] = min(abs(PowerSpec(divisionStart:k)-avgDivPwr));
        spectrumIndices(divisionIndex) = index+divisionStart-1;
       divisionDf(divisionIndex) = (frequency(k)-frequency(divisionStart));

    end
    


for k = 1:numComponents % cycle through the choices of frequency
    DistVals.power(k) = PowerSpec(spectrumIndices(k));
    DistVals.amplitude(k) = sqrt(2*DistVals.power(k)*divisionDf(k));
    DistVals.frequency(k) = frequency(spectrumIndices(k));
    FexGain(k) = interp1(SpectrumData.frequency,SpectrumData.body(bodyNumber).FeMag,frequency(spectrumIndices(k)));
    DistVals.force(k) = DistVals.amplitude(k)*FexGain(k);
end

end % terminates 1 frequency condition
DistVals.W = 2*pi.*DistVals.frequency;
end

