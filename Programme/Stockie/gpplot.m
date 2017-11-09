%
% gpplot.m: Plot the Gaussian plume solution with constant eddy diffusivity.
%
%
Q = 1;     % source emission rate (kg/s)
U = 1;     % wind speed (m/s)
H = 0;     % source height (m)
K = 1;     % eddy diffusion coefficient (m^2/s)

rmax = 2;  % size of domain in r (m)
xmax = rmax*U/K;
ymax = 5;  % size of domain in y (m)
zmax = 10; % size of domain in z (m)
N = 400;   % number of plotting points

% Ground-level maximum concentration:
x0 = U*H^2/(4*K)           % x coordinate 
c0 = 2*Q/(pi*U*H^2*exp(1)) % concentration value

% Plot the centerline concentration (y=0):
[rr,zz] = meshgrid( [1:N]*rmax/N, [1:N]*zmax/N );
xx = rr*U/K;
cc = Q./(4*pi*U*rr) .* ...
     ( exp(-(zz-H).^2./(4*rr)) + exp(-(zz+H).^2./(4*rr)) ); 
clist = [1e-5, 1e-4, 1e-3, 1e-2, 0.025, 0.05, 0.1];
[cs,ch] = contourf( xx, zz, cc, clist );
colormap jet 
clabel(cs)
xlabel('x'), ylabel('z')
hold on 
plot(x0+eps, 0+eps, 'ko')
plot(0+eps,  H+eps, 'rs')
arrow( [0.05*xmax, 0.8*zmax], [0.2*xmax,0.8*zmax], 'Width', 2 )  
text( 0.05*xmax, 0.88*zmax, 'wind' )
hold off
set(gca,'XLim', xmax*[-0.01,1])
shg
print -djpeg 'gplume1.jpg'
print -depsc 'gplume1.eps'

pause

% Plot the ground-level concentration (z=0);
[rr,yy] = meshgrid( [1:N]*rmax/N, [-N:N]*ymax/N );
xx = rr*U/K;
cc = 2 * Q./(4*pi*U*rr) .* exp(-yy.^2./(4*rr)) .* ...
     exp(-H^2./(4*rr));
[cs,ch] = contourf( xx, yy, cc, clist );
clabel(cs)
xlabel('x'), ylabel('y')
hold on
plot(x0+eps, 0+eps, 'ko')
plot(0+eps,  0+eps, 'rs')
arrow( [0.02*xmax, 0.7*ymax], [0.15*xmax,0.7*ymax], 'Width', 2 )  
text( 0.02*xmax, 0.85*ymax, 'wind' )
hold off
set(gca,'XLim', [-0.01,1]*xmax)
colorbar
shg
print -djpeg 'gplume2.jpg'
print -depsc 'gplume2.eps'