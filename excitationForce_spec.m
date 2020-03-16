function [ force ] = excitationForce_spec(time)
% Excitation force that is a sum of sine waves multiplied by magnitude

global SpectrumData
FX = SpectrumData.body(1).FeMag;
FxPhase = SpectrumData.body(1).FePhase;
W = SpectrumData.W;
Amp = SpectrumData.amp;
Phi = SpectrumData.phase;
dW = SpectrumData.dW;


for k = 1:length(W)

%   SpectrumData.magnitude has already been modified to be the
%   amplitude with sqrt(2*S*dW)
    product(k) = Amp(k)*FX(k)*exp(1i*(W(k)*time+FxPhase(k)+Phi(k))); % AS FAR AS I KNOW, THIS IS CORRECT!
% I am believe that + FxPhase is appropriate for how NEMOH calculates the
% excitation force, meaning phase initially declines.
    
%     product(k) = sqrt(2*Mag(k)*dW)*FX(k)*exp(1i*(W(k)*time+Phi(k)));
%       product(k) = Amp(k)*FX(k)*exp(1i*(W(k)*time+Phi(k))); % used to
%       duplicate JP1


end
force = imag(sum(product)); % imag extracts the sine component...
% force = real(trapz(product)); % duplicates the original submission for JP1

end

