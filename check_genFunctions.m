function [Functions] = check_genFunctions(genFunString,bodyNumber,vars,dynamicsEq,outputEq,Paths)
%check_genFunctions: 1/18/18 is meant to check to see if there are
%functions present and if so, avoid the time consuming function generation
%step

if strcmp(genFunString,'gen')
    Functions = generateFunctions(bodyNumber,vars,dynamicsEq,outputEq,Paths);  % lowercase functions -> nonlinear
    
elseif strcmp(genFunString,'check')
    val1 = exist('f.m','file');
    val2 = exist('h.m','file');
    val3 = exist('Jacob_f.m','file');
    val4 = exist('Jacob_h.m','file');
    val = val1*val2*val3*val4; % if any val# is zero then val will be zero
    
    if val == 0 
        Functions = generateFunctions(bodyNumber,vars,dynamicsEq,outputEq,Paths);  % lowercase functions -> nonlinear
    else
        Functions.f = @f;
        Functions.h = @h;
        Functions.F = @Jacob_f;
        Functions.H = @Jacob_h;
    end
        
else
    display('invalid string!!!')
    return
end

    Functions.fODE = @fNonlin;

end

