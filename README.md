# FDTDtesting
MATLAB FDTD vs ADS MoM


%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%

- 3DFDTD_finalProj is the root directory
- All user functions are managed in 3DFDTD_finalProj
- The program is configured with any config file in root
- The plots can be checked with 'showPlots'
- The program is run with 'app' from command window
- Workspace and probe settings are set in app
- If error on run from 'config_manual', return to root
- Use object structs to check physical parameters
- 'Rotate 3D' in figure toolbar unlocks 3D movement
- ExtractParams is not robust and must be modified based on config


%%%%%%%%%%%%%%%%% SHOW RESULTS %%%%%%%%%%%%%%%%%%%

1. Input configurations to config files
2. Set quarter wave value, simtime and sourceport in app
3. Run app 
4. After running simulation ExtractParams may be run 
   at any time to reload objects for use and show circuit
    -'obj_180hybrid' may be used in place of nearly all of 
    the previous workspace's variables in static calcs

%%%%%%%%%%%%%%%%% RUN PROGRAM %%%%%%%%%%%%%%%%%%%%

interface use:     
    -   app;
    -   ExtractParams;
    -   config_<name>;
output object:
    - 'obj_180hybrid' the first element is the line structures 
    and the second is the 180hybrid's structures
configFile(string var): 
    -   the saved config_<name>.m file used for object design. 
        Current inputs are 'manual', 'line' and '180hybrid'
    -   Config files can be run independently to check object
    -   If new config, add conditional to run from the 'app'
        'program/init_files/ReadInputData.m'
    - For network parameters, the line must match the test object
    - Use 180hybrid as template for other configurations
probe(uint8 vector): 
    -   uint8 array of probe numbers for existing probes
sourceport(uint8 vector):
    -   port(s) excited for testing
clrmess: 
    -   bool value. true to keep all workspace false 
        to keep only objects. Not suggested to use in dev
load <name>_wkspc.mat: 
    -   used to load workspace saved in program where 
        'name' is same as . Useful for revisiting results

%%%%%%%%%%%%%%  COMMONLY USED  %%%%%%%%%%%%%%%%%%%

- 'app'
- 'ExtractParams'
- 'load 180hybrid_wkspc.mat'
- 'load line_wkspc.mat'
- 'obj_180hybrid(2).trace(#)'
- 'obj_180hybrid(2).port(#)'
- 'obj_180hybrid(2).source(#)'
- 'obj_180hybrid(2).probe(#)'

