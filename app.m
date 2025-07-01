%% reset space and clear artifacts
close all; clc; clear;

%% program option
runProg = "N-Params";
% runProg = "oneshot";

%% program settings

if runProg == "oneshot"    
    % oneshot Mode
    opFreq = 3e9;  
    qwav = 380;%mil, quarter wavelength pi/2
    simtime = 1e-09;
    config = "line";
    mode = "single";
    sourceport = 1; %sourceport is array and sourcePort is scalar part of sourceport
    sourcePort = sourceport;%since scalar for line test
    fftprobe = [1,2]; %all ports listed is all ports plotted for each source looped
    oneshot;
    

elseif runProg == "N-Params"
    % N-parame settings
    opFreq = 3e9;  
    qwav = 380;%mil, quarter wavelength pi/2
    simtime = 1e-09;
    mode = "single"; 
    config = "line";
    NetParams;
    

elseif runProg == "thinwire"
    disp("thinwire not designed")
elseif runProg == "conformal"
    disp("conformal not designed")
end