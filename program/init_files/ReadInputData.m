%load physical data
cd ../..; %back to config at UI
if config == "180hybrid"
        config_180hybrid;
    elseif config == "line"
        config_line;
    elseif config == "manual"
        config_manual;
    elseif confdtig == "manual2"
        config_manual2;
    elseif config == "3GHz_180hybrid"
        config_3GHz_180hybrid;
    elseif config == "3GHz_180hybrid_correct"
        config_3GHz_180hybrid_correct;
    elseif config == "good"
        config_3GHz_180hybrid_good;
    elseif config == "BroadBand"
        config_3GHz_180hybrid_BroadBand;
end
cd program/init_files; %return to internal config

%  Reads the PML parameters from text file
fName = 'PMLparams.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bounding box data:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = mesh(1).dx;
dy = mesh(1).dy;
dz = mesh(1).dz;

N_x = mesh(1).xmax;
N_y = mesh(1).ymax;
N_z = mesh(1).zmax;

% global grid dimensions:
D_x = dx*double(N_x-1);
D_y = dy*double(N_y-1);
D_z = dz*double(N_z-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step limit:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numOfPEC = length(PEC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PEC Blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CFLN = mesh(1).CFLN;
totSimTime = mesh(1).simtime;

% initialize bounds arrays:
if(numOfPEC > 0)
    iPEC1 = zeros(numOfPEC,1);
    jPEC1 = zeros(numOfPEC,1);
    kPEC1 = zeros(numOfPEC,1);

    iPEC2 = zeros(numOfPEC,1);
    jPEC2 = zeros(numOfPEC,1);
    kPEC2 = zeros(numOfPEC,1);
end
% Read in the bounds of each block:
for ii = 1:numOfPEC
    iPEC1(ii) = floor(PEC(ii).i(1));
    jPEC1(ii) = floor(PEC(ii).j(1));
    kPEC1(ii) = floor(PEC(ii).k(1));

    iPEC2(ii) = floor(PEC(ii).i(2));
    jPEC2(ii) = floor(PEC(ii).j(2));
    kPEC2(ii) = floor(PEC(ii).k(2));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dielectric Blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numOfDielBlocks = length(substrate);  % Number of Dielectric 
                                                   % Blocks
% initialize bounds arrays:
if(numOfDielBlocks > 0)
    iEps1 = zeros(numOfDielBlocks,1);
    jEps1 = zeros(numOfDielBlocks,1);
    kEps1 = zeros(numOfDielBlocks,1);
    iEps2 = zeros(numOfDielBlocks,1);
    jEps2 = zeros(numOfDielBlocks,1);
    kEps2 = zeros(numOfDielBlocks,1);
    epsrB = zeros(numOfDielBlocks,1);
    sigmaB= zeros(numOfDielBlocks,1);
end
% Read in the bounds of each block and relative permittivity & 
% conductivity of each block:
for i = 1:numOfDielBlocks
    iEps1(i) = substrate(i).i(1);
    jEps1(i) = substrate(i).j(1);
    kEps1(i) = substrate(i).k(1);
    iEps2(i) = substrate(i).i(2);
    jEps2(i) = substrate(i).j(2);
    kEps2(i) = substrate(i).k(2);
    epsrB(i) = substrate(1).epsr;
    sigmaB(i)= substrate(1).sig;
end

% Open the excel file worksheet:
fid = fopen(fName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PML Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  X-Minus  (left side of mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerXm = fgetl(fid);
npmlXMinus = uint16(str2double(fgetl(fid)));  % PML Depth in cells
mScaleXMinus    = str2double(fgetl(fid));     % mscaling of sigma, kappa
pSigmaOptXMinus = str2double(fgetl(fid));     % % of sigma-opt
kappaMaxXMinus  = str2double(fgetl(fid));     % kappa-max
aMaxXMinus      = str2double(fgetl(fid));     % a-max
aScalXMinus     = str2double(fgetl(fid));     % a scale
epsrEffXMinus   = str2double(fgetl(fid));     % epsilon effective
%%%%%%%  X-Plus  (right side of mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerXp = fgetl(fid);
npmlXPlus = uint16(str2double(fgetl(fid)));  % PML Depth in cells
mScaleXPlus    = str2double(fgetl(fid));     % mscaling of sigma, kappa
pSigmaOptXPlus = str2double(fgetl(fid));     % % of sigma-opt
kappaMaxXPlus  = str2double(fgetl(fid));     % kappa-max
aMaxXPlus      = str2double(fgetl(fid));     % a-max
aScalXPlus     = str2double(fgetl(fid));     % a scale
epsrEffXPlus   = str2double(fgetl(fid));     % epsilon effective

%%%%%%%  Y-Minus  (left-y side of mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerYm = fgetl(fid);
npmlYMinus = uint16(str2double(fgetl(fid)));  % PML Depth in cells
mScaleYMinus    = str2double(fgetl(fid));     % mscaling of sigma, kappa
pSigmaOptYMinus = str2double(fgetl(fid));     % % of sigma-opt
kappaMaxYMinus  = str2double(fgetl(fid));     % kappa-max
aMaxYMinus      = str2double(fgetl(fid));     % alhpa-max
aScalYMinus     = str2double(fgetl(fid));     % a scale
epsrEffYMinus   = str2double(fgetl(fid));     % epsilon effective
%%%%%%%  Y-Plus  (right-y side of mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerYp = fgetl(fid);
npmlYPlus = uint16(str2double(fgetl(fid)));  % PML Depth in cells
mScaleYPlus    = str2double(fgetl(fid));     % mscaling of sigma, kappa
pSigmaOptYPlus = str2double(fgetl(fid));     % % of sigma-opt
kappaMaxYPlus  = str2double(fgetl(fid));     % kappa-max
aMaxYPlus      = str2double(fgetl(fid));     % alhpa-max
aScalYPlus     = str2double(fgetl(fid));     % a scale
epsrEffYPlus   = str2double(fgetl(fid));     % epsilon effective

%%%%%%%  Z-Minus  (ceiling of mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerZm = fgetl(fid);
npmlZMinus = uint16(str2double(fgetl(fid)));  % PML Depth in cells
mScaleZMinus    = str2double(fgetl(fid));     % mscaling of sigma, kappa
pSigmaOptZMinus = str2double(fgetl(fid));     % % of sigma-opt
kappaMaxZMinus  = str2double(fgetl(fid));     % kappa-max
aMaxZMinus      = str2double(fgetl(fid));     % alhpa-max
aScalZMinus     = str2double(fgetl(fid));     % a scale
epsrEffZMinus   = str2double(fgetl(fid));     % epsilon effective
%%%%%%%  Z-Plus  (floor of mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
headerZp = fgetl(fid);
npmlZPlus = uint16(str2double(fgetl(fid)));  % PML Depth in cells
mScaleZPlus    = str2double(fgetl(fid));     % mscaling of sigma, kappa
pSigmaOptZPlus = str2double(fgetl(fid));     % % of sigma-opt
kappaMaxZPlus  = str2double(fgetl(fid));     % kappa-max
aMaxZPlus      = str2double(fgetl(fid));     % alhpa-max
aScalZPlus     = str2double(fgetl(fid));     % a scale
epsrEffZPlus   = str2double(fgetl(fid));     % epsilon effective


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source Excitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numSources = length(source);  % Number of Current
                                              % sources

% initialize bounds arrays:
if(numSources > 0)
    sourceType = zeros(numSources,1);
    sourceSigType = zeros(numSources,1);
    srct_w  = zeros(numSources,1);
    srct_o  = zeros(numSources,1);
    srcfreq = zeros(numSources,1);
    isrc1 = zeros(numSources,1);
    jsrc1 = zeros(numSources,1);
    ksrc1 = zeros(numSources,1);
    isrc2 = zeros(numSources,1);
    jsrc2 = zeros(numSources,1);
    ksrc2 = zeros(numSources,1);
    polInc  = zeros(numSources,1);
    samp   = zeros(numSources,1);
    thetInc = zeros(numSources,1);
    phiInc  = zeros(numSources,1);
    resInc  = zeros(numSources,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the source type, the bounds, and amplitude.  Also, read
% in specific source information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  typeOfSource = string indicating source type (see below)
%  polInc = the source polarization (x, y, or z for directed
%           sources, or vertical or horizontal for plane wave source (not
%           implemented
%  samp   = the source amplitude (multiplied times the signature)
%  (i1,j1,k1), (i2,j2,k2) = nodes of bounds of the source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source Signature  (for each source, s.t., each source can have a different
%                    time signature)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sourceSig = 'Gaussian' for gaussian pulse
%           = 'DiffGaussian' for differentiated Gaussian pulse
%           = 'ModGaussian' for modulated Gaussian pulse
%       t_w = pulse width in seconds
%      nt_o = delay multiplier (t_o = nt_o * t_w)
%      freq = modulating frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numSources

    typeOfSource = source(i).type;
    sourcePolarization = source(i).pol;
    if(strcmp(typeOfSource,'JCurrent'))
        sourceType(i) = 1;
    elseif(strcmp(typeOfSource,'MCurrent'))
        sourceType(i) = 2;
    elseif(strcmp(typeOfSource,'PlaneWave'))
        sourceType(i) = 3;
        thetInc(i) = source(i).thetInc;
        phiInc(i) = source(i).phiInc;
    elseif(strcmp(typeOfSource,'DeltaGap')) 
        sourceType(i) = 4;
    elseif(strcmp(typeOfSource,'HardVoltage')) 
        sourceType(i) = 5;
    elseif(strcmp(typeOfSource,'SoftVoltage'))
        sourceType(i) = 6;
    elseif(strcmp(typeOfSource,'Thevenin')) 
        sourceType(i) = 7;
        resInc(i) = source(i).resInc;  % Thevenin resistance
    else
        error('Illegal source type');
    end
    
    if(strcmp(sourcePolarization,'x'))   % Axis of current or voltage excitation
        polInc(i) = 1;
    elseif(strcmp(sourcePolarization,'y'))
        polInc(i) = 2;
    elseif(strcmp(sourcePolarization,'z'))
        polInc(i) = 3;
    elseif(strcmp(sourcePolarization,'Vertical'))  % If PW excitation
        polInc(i) = 1;
    elseif(strcmp(sourcePolarization,'Horizontal'))
        polInc(i) = 2;
    else
        error('Illegal source polarization type');
    end
    
    % source amplitude (weights the source signature)
    samp(i)  = source(i).amp;
    
    % (i,j,k) primary grid bounds of the source
    isrc1(i) = source(i).i(1);
    jsrc1(i) = source(i).j(1);
    ksrc1(i) = source(i).k(1);
    isrc2(i) = source(i).i(2);
    jsrc2(i) = source(i).j(2);
    ksrc2(i) = source(i).k(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  read in the source signature specifics for this source
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sourceSig = source(i).sig;
    srct_w(i)  = source(i).tw;
    nt_o = source(i).nto;
    srct_o(i) = nt_o * srct_w(i);
    
    if(strcmp(sourceSig,'Gaussian'))
        sourceSigType(i) = 1;
    elseif(strcmp(sourceSig,'DiffGaussian'))
        sourceSigType(i) = 2;
    elseif(strcmp(sourceSig,'ModGaussian'))
        sourceSigType(i) = 3;
        srcfreq(i) = source(i).freq;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Currently, user has the ability to output probed values of
%  fields (Ex,Ey,Ez,Hx,Hy,Hz), voltage, or current as a function of time.
%  Others can be added, such as NFFF.  The output can be plotted at 
%  the end of the run.  All output probes are dumped to a file vs time.
%  Additional directives could easily be added to FFT probes as well.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numOutputQty = length(probe);  % Number of Output
                                                 % fields

% initialize bounds arrays:
if(numOutputQty > 0)
    iOutput1 = zeros(numOutputQty,1);
    jOutput1 = zeros(numOutputQty,1);
    kOutput1 = zeros(numOutputQty,1);
    iOutput2 = zeros(numOutputQty,1);
    jOutput2 = zeros(numOutputQty,1);
    kOutput2 = zeros(numOutputQty,1);
    OutputType = zeros(numOutputQty,1);
    OutputAxis = zeros(numOutputQty,1);
    outputPlot = zeros(numOutputQty,1);
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the bounds of each output quantity type
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numOutputQty
    typeOfOutput = probe(i).type;
    if(strcmp(typeOfOutput,'Ex'))
        OutputType(i) = 1;
    elseif(strcmp(typeOfOutput,'Ey'))
        OutputType(i) = 2;
    elseif(strcmp(typeOfOutput,'Ez'))
        OutputType(i) = 3;
    elseif(strcmp(typeOfOutput,'Hx')) 
        OutputType(i) = 4;
    elseif(strcmp(typeOfOutput,'Hy')) 
        OutputType(i) = 5;
    elseif(strcmp(typeOfOutput,'Hz'))
        OutputType(i) = 6;
    elseif(strcmp(typeOfOutput,'Voltage')) 
        OutputType(i) = 7;
        % Define the axis for integration to evaluate voltage
        voltageAxis = probe(i).axis;
        if(strcmp(voltageAxis,'x'))
            OutputAxis(i) = 1;
        elseif(strcmp(voltageAxis,'y'))
            OutputAxis(i) = 2;
        elseif(strcmp(voltageAxis,'z'))
            OutputAxis(i) = 3;
        else
            error('Illegal axis type');
        end
    elseif(strcmp(typeOfOutput,'Current')) 
        OutputType(i) = 8;
        % Define the axis of current flow being measured
        % (this dictates the normal axis of the surface being integrated)
        currentAxis = probe(i).axis;
        if(strcmp(currentAxis,'x'))
            OutputAxis(i) = 1;
        elseif(strcmp(currentAxis,'y'))
            OutputAxis(i) = 2;
        elseif(strcmp(currentAxis,'z'))
            OutputAxis(i) = 3;
        else
            error('Illegal axis type');
        end
    
    else
        error('Illegal field probe type');
    end

    %  State if this should be plotted after simulation complete
    %    0 = don't plot
    %    1 = plot
    plotNoPlot = probe(i).plot;
    if(strcmp(plotNoPlot,'plot'))
        outputPlot(i) = 1;
    elseif(strcmp(plotNoPlot,'noplot'))
        outputPlot(i) = 0;
    else
        error('Illegal plot/noplot statement');
    end
    
    % i,j,k primary grid bounds of the output quantity
    %  Field bound descriptions:
    %   E fields: E-fields at edge centers, so will need bounds
    %             of the edge(s)  (one ahead of the field, one behind)
    %             Example:  ex(5,5,5) will have bounds (5,5,5), (6,5,5)
    %
    %   H-fields: H-fieldsare at face centers, so will need bounds of
    %             the face(s)  
    %             Example: Hx could have bounds (5,5,5),(5,6,6))
    %
    %   Voltage:  Voltages are integrals over E-field edges
    %             and expect bounds along a cardinal axis.  The bounds
    %             are nodes of the primary grid.
    %
    %   Current:  Current bounds are defined by the primary grid bounds
    %             of the conductor which is carrying the current.  The
    %             integral is done over secondary grid faces 1/2 cell AHEAD
    %             of the primary grid bounds along the axis of current
    %             flow.
    %             
    iOutput1(i) = probe(i).i(1);
    jOutput1(i) = probe(i).j(1);
    kOutput1(i) = probe(i).k(1);
    iOutput2(i) = probe(i).i(2);
    jOutput2(i) = probe(i).j(2);
    kOutput2(i) = probe(i).k(2);
    
end

clear header*;

% close file
fclose(fid);
%end