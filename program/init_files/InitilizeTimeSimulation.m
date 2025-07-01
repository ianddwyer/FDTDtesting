% This function Initializes the grid dimensions and the time step
%
% Free-space Material property definitions:
eps0 = 8.854187817620389850544e-12; % free-space permittivity (F/m)
mu0 = 12.56637061435917295384e-7;   % free-space permeability (H/m)

c0 = 1/sqrt(mu0*eps0);  % speed of light (m/s)
eta0 = sqrt(mu0/eps0);  % free space wave impedance (ohms)

% Compute the local grid cell sizes:
% dx = double(D_x)/(double(N_x)-1.0);
% dy = double(D_y)/(double(N_y)-1.0);
% dz = double(D_z)/(double(N_z)-1.0);

nx = N_x;
ny = N_y;
nz = N_z;

% Compute the time step:
dt = CFLN/(c0*sqrt(1.0/dx^2 + 1.0/dy^2 + 1.0/dz^2));

% Number of time steps:
nMax = floor(totSimTime/dt)

disp("");
