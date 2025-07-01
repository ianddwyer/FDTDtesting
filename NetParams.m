try
    
    %% RUNS PROGRAM 
    %SETUP MODEL IN config file
    %must be executed at root of app
    
    %fftprobe is number of lines per sim
    %sourceport runs a sim for each sourceport
    %1:4 on both sourceport and fft probe returns 16 traces with 4 sims
    %clearvars -except config clrmess fftprobe sourceport mode
    
    appDir = cd; %hold location where app dir is saved for later use
    
    %init object structures for outputs
    holdOutput = struct([]);
    holdOutput_line = struct([]);
    fVecs = struct([]);
    tVecs = struct([]);
    linespecs = struct([]); %holds line dimensions
    
    
    %obtain line parameters for reference
    config = "line";
    fftprobe = 1; %all ports listed is all ports plotted for each source looped
    sourceport = 1; %sourceport is array and sourcePort is scalar part of sourceport
    sourcePort = sourceport;%since scalar for line test
    mode = "single"; 
    wrkspc_name = sprintf(string(config)+"_wkspc.mat");
    
    cd program;
    Main3DFDTD_Accelerated; %runs main with new configuration for sourcePort 
    cd ..;
    
    save(wrkspc_name);

    %setup actual circuit to test
    config = "180hybrid";
    fftprobe = [1,2,3,4]; %all ports listed is all ports plotted for each source looped
    sourceport = [1,2,3,4]; %sourceport is array and sourcePort is scalar part of sourceport
    mode = "Sparam"; 
    wrkspc_name = sprintf(string(config)+"_wkspc.mat");
    
    
    
    %tests all ports claimed 
    %sourceport is array and sourcePort is scalar part of sourceport
    for sourcePort=sourceport %sets sourceport to run through fftprobes set
 
        cd program;
        Main3DFDTD_Accelerated; %runs main with new configuration for sourcePort
        cd ..;
    
    end
    save(wrkspc_name);
    
    ExtractParams;
    
    
    
    % clr_wkspc;
    
    %% Notes to self
    % disp(" ");
    % disp("NOTES TO SELF:")
    % disp("-Get S, V+, V-, Zin, Zo, epsEff, ")
    % disp("-PML configurations needed, effective values used...")
    % disp("-Automate port switching for test in app function")
    % disp("-Associate port objects with source and probes in readinputdata")
    % disp("-Use calculated values in ratio from MW circ final project")
    % disp("-Get Gedney's help with:")
catch error_msg
    cd(appDir) %returns to app if error
    save error_wkspc.mat
    rethrow(error_msg)
end