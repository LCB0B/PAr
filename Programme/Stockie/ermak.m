function C = ermak( x, y, z, H, Q, U, Wset, Wdep )
% ERMAK: Compute contaminant concentration (kg/m^3) using the
%   Gaussian plume model, modified for a deposition and settling
%   velocity.  This code handles a single source (located at the
%   origin) and multiple receptors.
%
% Input parameters:
%
%     x - receptor locations: distance along the wind direction, with
%         the source at x=0 (m) 
%     y - receptor locations: cross-wind direction (m)
%     z - receptor locations: vertical height (m)
%     H - source height (m)
%     Q - contaminant emission rate (kg/s)
%     U - wind velocity (m/s)
%     Wset - gravitational settling velocity (m/s)
%     Wdep - deposition velocity (m/s)
%
% Output:
%
%     C - contaminant concentration (kg/m^3)
%
% References: Ermak (1977), Winges (1990/1992).

% First, define the cut-off velocity, below which concentration = 0.
Umin      = 0.0;

% Determine the sigma coefficients based on stability class C --
% slightly unstable (3-5 m/s). 
ay = 0.34;  by = 0.82;  az = 0.275; bz = 0.82;  
sigmay = ay*abs(x).^by .* (x > 0);
sigmaz = az*abs(x).^bz .* (x > 0);

% Calculate the eddy diffusivity (m^2/s).
Kz = 0.5*az*bz*U*abs(x).^(bz-1) .* (x > 0);  % K = 0.5*U*d(sigma^2)/dx

% Calculate the contaminant concentration (kg/m^3) using Ermak's formula.
if U < Umin,
  C = 0 * z;
else
  Wo = Wdep - 0.5*Wset;
  C  = Q ./ (2*pi*U*sigmay.*sigmaz) .* exp( -0.5*y.^2./sigmay.^2 ) .* ...
       exp( -0.5*Wset*(z-H)./Kz - Wset^2*sigmaz.^2/8./Kz.^2 ) .* ...
       ( exp( -0.5*(z-H).^2./sigmaz.^2 ) + ...
	 exp( -0.5*(z+H).^2./sigmaz.^2 ) - sqrt(2*pi)*Wo*sigmaz./Kz .* ...
	 exp( Wo*(z+H)./Kz + 0.5*Wo^2*sigmaz.^2./Kz.^2 ) .* ...
	 erfc( Wo*sigmaz/sqrt(2)./Kz + (z+H)./sqrt(2)./sigmaz ) );
  ii = find(isnan(C) | isinf(C));
  C(ii) = 0;   % Set all NaN and inf values to zero.
end