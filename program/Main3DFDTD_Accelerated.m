% 3D FDTD Code for MATLAB
% Author:   Stephen D. Gedney, University of Colorado Denver
%                              Electrical Engineering Department
% Date:     11/9/2015
% Distribution:   Authorized for distribution only to UCD students
% Disclaimer:     This code is not claimed to be error free.  Please
%                 report any bugs to stephen.gedney@ucdenver.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = tic;
% Pre-processing:
%--------------------------------------------



cd init_files;
% Input mesh dimensions, materials, conformal boundaries, sources,
% time-signature, simulation control, output control, etc.
ReadInputData;
% Initializations:
%--------------------------------------------
% initialize the time-simulation (time-steps, simulation time)
InitilizeTimeSimulation;
% Average the material parameters:
AverageMaterialParams;
% initialize PML parameters:
InitializePML;
% initialize main update coefficients:
InitializeUpdateCoefficients;
% intialize the fields and auxiliary variables throughout the grid:
InitializeFields;
% TimeStepping:
%--------------------------------------------
% Advance the fields through nMax time iterations:
cd ..;



cd loop;
for n = 1:nMax
    clc
    fprintf('\n%s mode \ntesting source: %i \n%3.0f%%  \n',config,sourcePort,n/nMax*100);
    
    % Explicit update of the electric fields:
    %---------------------------------------
    % update electric fields throughout the entire domain:
    ex(1:nx-1,2:ny-1,2:nz-1) = ex(1:nx-1,2:ny-1,2:nz-1) +...
        (hz(1:nx-1,2:ny-1,2:nz-1)-hz(1:nx-1,1:ny-2,2:nz-1)-...
         hy(1:nx-1,2:ny-1,2:nz-1)+hy(1:nx-1,2:ny-1,1:nz-2)).*...
         coefex(1:nx-1,2:ny-1,2:nz-1);

    ey(2:nx-1,1:ny-1,2:nz-1) = ey(2:nx-1,1:ny-1,2:nz-1) +...
        (hx(2:nx-1,1:ny-1,2:nz-1)-hx(2:nx-1,1:ny-1,1:nz-2)-...
         hz(2:nx-1,1:ny-1,2:nz-1)+hz(1:nx-2,1:ny-1,2:nz-1)).*...
         coefey(2:nx-1,1:ny-1,2:nz-1);
       
    ez(2:nx-1,2:ny-1,1:nz-1) = ez(2:nx-1,2:ny-1,1:nz-1) +...
        (hy(2:nx-1,2:ny-1,1:nz-1)-hy(1:nx-2,2:ny-1,1:nz-1)-...
         hx(2:nx-1,2:ny-1,1:nz-1)+hx(2:nx-1,1:ny-2,1:nz-1)).*...
         coefez(2:nx-1,2:ny-1,1:nz-1);
    
    % Auxiliary electric field updates (material, local, PML):
    EPMLFieldUpdates;
    
    % Electric field source injection:
    ETypeSourceInjection;
    
    % Explicit update of the magnetic fields:
    %---------------------------------------
    % update magnetic fields throughout the entire domain:
    hx(1:nx,1:ny-1,1:nz-1) =  hx(1:nx,1:ny-1,1:nz-1)-...
        (ez(1:nx,2:ny,1:nz-1)-ez(1:nx,1:ny-1,1:nz-1)-...
         ey(1:nx,1:ny-1,2:nz)+ey(1:nx,1:ny-1,1:nz-1)).*...
         coefhx(1:nx,1:ny-1,1:nz-1);
        
    hy(1:nx-1,1:ny,1:nz-1)  = hy(1:nx-1,1:ny,1:nz-1)-...
        (ex(1:nx-1,1:ny,2:nz)-ex(1:nx-1,1:ny,1:nz-1)-...
         ez(2:nx,1:ny,1:nz-1)+ez(1:nx-1,1:ny,1:nz-1)).*....
         coefhy(1:nx-1,1:ny,1:nz-1);
       
    hz(1:nx-1,1:ny-1,1:nz)  = hz(1:nx-1,1:ny-1,1:nz)-...
        (ey(2:nx,1:ny-1,1:nz)-ey(1:nx-1,1:ny-1,1:nz)-...
         ex(1:nx-1,2:ny,1:nz)+ex(1:nx-1,1:ny-1,1:nz)).*...
         coefhz(1:nx-1,1:ny-1,1:nz);

    % Auxiliary magnetic field updates (material, local, PML):
    HPMLFieldUpdates;
    
    % Magnetic field source injection:
    %HTypeSourceInjection;
    
    %processing
    OutputData;    
end
cd ..;



cd post_processing\;
PostProcess;   
tfin = toc(tstart);
fprintf('Simulation complete after %6.2f seconds \n', tfin');

cd ..; %back to main

