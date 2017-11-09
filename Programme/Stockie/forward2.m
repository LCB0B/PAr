% FORWARD2: Solve the forward atmospheric dispersion problem using
%    **Ermak's solution** with deposition and settling.  That is, given 
%    the source emission rates for zinc (in kg/s), calculate and plot the
%    ground-level Zn concentration (in mg/m^3) and the total deposition 
%    (in mg) at a given set of receptor jar locations.

clear all
setparams;   % read parameters from a file
Uwind = 5;   % wind speed (m/s) 
dt    = 30 * (24 * 3600);  % on month (in seconds)

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
dep = 0;
warning( 'OFF', 'MATLAB:divideByZero' );
for i = 1 : source.n, 

  % Sum up ground-level Zn concentrations from each source at all mesh
  % points, shifting the (x,y) coordinates so the source location is
  % at the origin.
  glc = glc + ermak( xmesh-source.x(i), ymesh-source.y(i), 0.0, ...
                     source.z(i), source.Q(i), Uwind, Wdep, Wset );

  % Do the same for the concentration at each receptor location and
  % scale the result by (A * dt * Wdep) to obtain a total deposition
  % (in mg) over the time interval 'dt'.
  dep = dep + (A * dt * Wdep) * ...
        ermak( recept.x-source.x(i), recept.y-source.y(i), recept.z, ...
               source.z(i), source.Q(i), Uwind, Wdep, Wset );

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

% Draw and label the source and receptor locations.
hold on
plot( source.x, source.y, 'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r' )
text( source.x, source.y, source.label, 'FontSize', smallfont, 'FontWeight','bold' );
greencolor = [0,0.816,0];
plot( recept.x, recept.y, 'g^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', greencolor )
text( recept.x, recept.y, recept.label, 'FontSize', smallfont, 'FontWeight', 'bold' )
hold off
shg
print -djpeg 'glc.jpg'

% Generate a second plot showing mass of Zn deposited at each receptor.
dep2 = dep * 1e6;   % deposition in mg
fprintf( 1, '\nDeposition in each receptor (mg):\n' );
disp(dep)
fprintf( 1, 'Total deposited in all receptors (mg):\n' );
disp(sum(dep))

figure(2)
bar(dep)
set(gca, 'XTick', [1:recept.n]);
set(gca, 'XTickLabel', recept.label)
xlabel('Receptor'), ylabel('Amount deposited (mg)')
grid on
set(gca,'XGrid', 'off')
shg
print -djpeg 'depbars.jpg'