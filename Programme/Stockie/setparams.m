% SETPARAMS: Set up various physical parameters for the atmospheric
%    dispersion problem.

% Contaminant parameters (zinc):
grav  = 9.8;          % gravitational acceleration (m/s^2)
mu    = 1.8e-5;       % dynamic viscosity of air (kg/m.s)
rho   = 7140;         % density of zinc (kg/m^3)
R     = 0.45e-6;      % diameter of zinc particles (m).  See Gatz (1975) 
Wdep  = 0.0062;       % Zn deposition velocity (m/s), in the range [5e-4,1e-2]  
Wset  = 2*rho*grav*R^2 / (9*mu); % settling velocity (m/s): Stokes law
Mzn   =  65.4e-3;     % molar mass of zinc (kg/mol)

% Other parameters:
dia   = 0.162;        % receptor diameter (m)
A     = pi*(dia/2)^2; % receptor area (m^2) 

% Stack emission source data:
source.n = 4;                         % # of sources
source.x = [288, 308, 900, 1093];     % x-location (m)
source.y = [ 77, 207, 293,  186];     % y-location (m)
source.z = [ 15,  35,  15,   15];     % height (m)
source.label=[' S1'; ' S2'; ' S3'; ' S4'];
tpy2kgps = 1.0 / 31536;               % conversion factor (tonne/yr to kg/s)
source.Q = [35, 80, 5, 5] * tpy2kgps; % emission rate (kg/s)

% Set locations of receptors where deposition measurements are made:
recept.n = 9;                                                 % # of receptors
recept.x = [  60,  76, 267, 331, 514, 904, 1288, 1254, 972 ]; % x location (m)
recept.y = [ 130,  70, -21, 308, 182,  75,  116,  383, 507 ]; % y location (m)
recept.z = [   0,  10,  10,   1,  15,   2,    3,   12,  12 ]; % height (m)
recept.label=[ ' R1 '; ' R2 '; ' R3 '; ' R4 '; ' R5 '; ' R6 '; ...
               ' R7 '; ' R8 '; ' R9 ' ];