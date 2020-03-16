function [DataOutput] = check_readBosmaData(dataString,name,Paths)
% check_readBosmaData: 1/11/18 is meant to run through the logic
% statements of whether or not to read Bosma data or not. This is meant
% to de-clutter the main scripts that are simulating the FPS experiment.
% NOTE: the sampling rates are hard coded into this logic checing function.

global Param
% samplingFreq.motion = Param.samplingFrew.motion;
% samplingFreq.wave = 128;
% samplingFreq.loadCell = 128;
mooringBool = [0]; % 0-> taught & 1-> catenary

if strcmp(dataString,'read')
    readBosmaData(name, Param.samplingFreq, mooringBool, Paths)
    cd(Paths.dataDir)
    load(strcat('BosmaData',char(name.runNo),'.mat'))
    cd(Paths.exampleDir)
elseif strcmp(dataString,'check')
    cd(Paths.dataDir)
    val = exist(strcat('BosmaData',char(name.runNo),'.mat'),'file');
    cd(Paths.exampleDir)
    if val == 0
        readBosmaData(name, Param.samplingFreq, mooringBool, Paths)
        cd(Paths.dataDir)
        load(strcat('BosmaData',char(name.runNo),'.mat'))
        cd(Paths.exampleDir)
    else
        cd(Paths.dataDir)
        load(strcat('BosmaData',char(name.runNo),'.mat'))
        cd(Paths.exampleDir)
    end
elseif stcmp(dataString,'sim') == 0
%     if function breaks on above line, did you specify a valid string?
end

DataOutput = Data;
end
