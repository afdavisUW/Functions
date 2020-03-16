function [bodyData] = check_generateBody(bodyString,name,bodyFunc,Paths,desiredOrder)
%check_generateBody: 1/11/18 is meant to perform the logic checks on
%wheather or not to run the generateBody function or to simply load the
%data. This is meant to clean the main script file from needing to have
%many logic statements that need to be copied each time the simulation is
%changed.  Note that this function will adhere to my naming convention.
%individual bodies are named body 1,2,... but then the data is given the
%end of the filename which is hard coded below.

for k = 1:length(name.body)
    if strcmp(bodyString,'run') == 1
        cd(Paths.functionsDir)
        bodyFunc(char(name.body(k)),Paths,desiredOrder)
        cd(Paths.exampleDir)
    elseif strcmp(bodyString,'check') == 1
        cd(Paths.outputDir)
        val = exist(strcat(char(name.body(k)),'_ABFe_freq.mat'),'file');
        cd(Paths.exampleDir)
        if val == 0
            cd(Paths.functionsDir)
            bodyFunc(char(name.body(k)),Paths,desiredOrder)
            cd(Paths.exampleDir)
        end
    elseif strcmp(bodyString,'skip') == 0
        display('!! Invalid generateBody string !!')
    end
    % loading body data
    cd(Paths.outputDir)
    load(strcat(char(name.body(k)),'_ABCD_ssReal.mat'))
    bodyData.body(k).A = A_r; 
    bodyData.body(k).B = B_r; 
    bodyData.body(k).C = C_r;
    bodyData.body(k).D = D_r;    
    clear A_r B_r C_r D_r
    load(strcat(char(name.body(k)),'FexIRF.mat'))
    bodyData.body(k).FIRF = F_ex.IRF; % this is already the real component
    bodyData.body(k).irfTime = F_ex.time;    
    clear F_ex
    load(strcat(char(name.body(k)),'_ABFe_freq.mat'))
    bodyData.body(k).FeMag = abs(Fe(:,3));
    bodyData.body(k).W = omega';
    bodyData.body(k).FePhase = angle(Fe(:,3));
    clear A B Fe omega
    cd(Paths.exampleDir)


end

end

