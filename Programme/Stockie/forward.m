% FORWARD: Solve the forward atmospheric dispersion problem using the
%    **standard Gaussian plume solution**.  That is, given the source
%    emission rates for Zn (in kg/s), calculate and plot the
%    ground-level Zn concentration (in mg/m^3). 

clear all
setparams;   % read parameters from a file
Uwind = 1;   % wind speed (m/s) 

% Set plotting parameters.
nx = 100;
ny = nx;
xlim = [   0, 2000];
ylim = [-100,  400];
x0 = xlim(1) + [0:nx]/nx * (xlim(2)-xlim(1)); % distance along wind direction (m)
y0 = ylim(1) + [0:ny]/ny * (ylim(2)-ylim(1)); % cross-wind distance (m)
[xmesh, ymesh] = meshgrid( x0, y0 );          % mesh points for contour plot
smallfont = 14;

glc = 0;
warning( 'OFF', 'MATLAB:divideByZero' );
for i = 1 : source.n, 
  % Sum up ground-level Zn concentrations from each source at all mesh
  % points, shifting the (x,y) coordinates so the source location is
  % at the origin.
  glc = glc + gplume( xmesh-source.x(i), ymesh-source.y(i), 0.0, ...
                      source.z(i), source.Q(i), Uwind );
end
warning( 'ON', 'MATLAB:divideByZero' );

% Plot contours of ground-level Zn concentration.
figure(1)
clist = [ 0.001, 0.01, 0.02, 0.05, 0.1 ]; 
glc2 = glc*1e6;  % convert concentration to mg/m^3
[c2, h2] = contourf( xmesh, ymesh, glc2, clist );
% axis equal     % (for plots in paper)
clabel(c2, h2, 'FontSize', smallfont-2 )
colormap(1-winter)  % These colors make the labels more readable
colorbar
set(gca, 'XLim', xlim ), set(gca, 'YLim', ylim )
xlabel('x (m)'), ylabel('y (m)')
title(['Zn concentration (mg/m^3), max = ', sprintf('%5.2f', max(glc2(:)))])
grid on

% Draw and label the source locations.
hold on
plot( source.x, source.y, 'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r' )
text( source.x, source.y, source.label, 'FontSize', smallfont, 'FontWeight','bold' );
hold off
shg
print -djpeg 'glc.jpg'