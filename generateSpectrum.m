function [SpectrumOutput] = generateSpectrum(waveBool,bodyData)
%generateSpectrum: 11/1/18 is a function to generate the Bretscneider
%spectrum for a body. this functoin has number of sine waves and min and
%max omega hardcoded and outputs data that can be used to form global data
%for wave generation 
global Param

num = 1000; % number of sine waves summed to generate irregular waves
minOmega = 0.15;
maxOmega = 30;
frequency = linspace(minOmega/(2*pi),maxOmega/(2*pi),num);
dW = (frequency(2)-frequency(1))*2*pi;
dF = frequency(2)-frequency(1);

for j = 1:length(frequency) % loop to calculate specta
    A = frequency(j)/(Param.wavePeriod^(-1));
    S(j) = 5*(Param.waveHs)^2/(16*(Param.wavePeriod)^(-1))*(1/A^5)*exp(-5/4*A^(-4));
  if waveBool == 0
      S(j) = 0; % option preserved to generate no wave
  end
  
  if waveBool == 2
      S(j) = 0; % this will be used to make a regular wave.
  end

end

SpectrumOutput.dW = dW;
SpectrumOutput.frequency = frequency;
% The work in JP1 shows that the use of dF does not produce realistic wave
% heights, dW does, as of 4/10/18 dW will be used. 4/19/18 This is undone
% SpectrumOutput.amp = sqrt(2*S.*dW);
SpectrumOutput.amp = sqrt(2*S.*dF);
SpectrumOutput.pwr = S;
SpectrumOutput.phase = rand([1,num])*2*pi;
SpectrumOutput.W = frequency.*2*pi;

if waveBool == 2
    vec = frequency-Param.peakFrequency;
    [val,index] = min(abs(vec));
    SpectrumOutput.amp(index) = 1;
end

% formatting frequency excitation to useful form
for k = 1:Param.numBodies
SpectrumOutput.body(k).FeMag = spline(bodyData.body(k).W,bodyData.body(k).FeMag,SpectrumOutput.W);
SpectrumOutput.body(k).FePhase = spline(bodyData.body(k).W,bodyData.body(k).FePhase,SpectrumOutput.W);



end

