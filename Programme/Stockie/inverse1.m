% INVERSE1: Solve the inverse atmospheric dispersion 
%    problem for a constant wind with no constraints.
%    (example in Section 4.1).
%    
%    Namely, given the amount of Zn deposited in each receptor
%    (in kg), estimate the source emission rates (in kg/s). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all parameter values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
setparams;            % read parameters

slist = [1];          % list of sources
Ns    = length(slist);% number of sources
rlist = [1:9];        % list of receptors
Nr    = length(rlist);% number of receptors

Uwind = 5.0;          % constant wind speed (m/s) 
dtime = 30*24*60*60;  % total time (sec) = one month
source.Q0 = 1 + 0*source.Q;  % dummy source values (kg/s)

% Receptor measurements (kg)
recept.Dzn = 1e-6 * [ 8.4, 68, 33, 4.2, 11, 8.2, 2.9, 2.2, 0.93 ]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve the inverse problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( 1, 'Constructing matrix ...\n' );
Glsq = zeros( Nr, Ns );

for is = 1 : Ns,
  % Calculate Ermak solution and accumulate in G matrix 
  % (units of G entries are s/m^3).
  warning( 'OFF', 'MATLAB:divideByZero' );
  ss = slist(is);
  Glsq(:,is) = Glsq(:,is) + ...
      ermak( recept.x(rlist) - source.x(ss), recept.y(rlist) - source.y(ss), ...
             recept.z(rlist), source.z(ss), source.Q0(ss), Uwind, Wset, Wdep )' / source.Q0(ss);
  warning( 'ON', 'MATLAB:divideByZero' );
end

% Set up the linear least squares solve, assuming Nr >= Ns:
%   * Glsq is an (Nr x Ns)-matrix,
%   * rhs is an (Nr x 1)-vector,
%   * solution Q is an (Ns x 1)-vector.
rhs = recept.Dzn(rlist) / (A * dtime * Wdep); 

% Solve the system Glsq * Q = rhs using linear least squares.  
[Q, res] = lsqlin( Glsq, rhs );
Q = Q';    % convert to a column vector

% Display Q, rescaled in units of T/yr.
fprintf( '\nZn emission rates (T/yr):\n' );
for is = 1 : Ns,
  fprintf( '  %2s: %9.4f\n', source.label(slist(is),:), Q(is)/tpy2kgps );
end
fprintf( '       ----------\nTotal: %9.4f T/yr\n', sum(Q)/tpy2kgps );
fprintf( 'LSQLIN residual = %e\n\n', res );