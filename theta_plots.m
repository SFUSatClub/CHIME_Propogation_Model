syms re ra l theta lamda
syms dre dra dl dtheta
rl = sqrt(re^2+(re+ra)^2-2*re*(re+ra)*cos(l))
phi = theta - asin((re/rl)*sin(l))
P = 2*exp(-0.5*(phi/0.445)^2) %+20*log(lamda/(4*pi*rl))

%the fractional error of power transmitted from satellite (not including propogation to ground)
frac = sqrt((diff(P, l)*dl)^2+(diff(P, ra)*dra)^2+(diff(P, theta)*dtheta)^2)

placed=subs(frac, [lamda re ra dra l], [700E6 6400E3 800E3 10000 0]) %worst case numbers

%surffuncdra=solve(placed==0.01, dtheta) %looking for 1% fractional error solution, only works when isolating dtheta

surffuncdra=solve(placed==0.01, dtheta)
syms x y
%surfacexy = subs(surffuncdra, [dl dra], [x y]) %converted to xy so that fsurf can plot it
%fsurf(surfacexy, [0 0.17 1000 10000], 'MeshDensity',10)

surfacexy = subs(surffuncdra, [dl theta], [x y]) %converted to xy so that fsurf can plot it
fsurf(surfacexy, [0 0.01 0 1.5], 'MeshDensity', 40)

% This is just changing the plot appearance please ignore
xlabel('dl (rad)','Interpreter','latex')
ylabel('theta (rad)','Interpreter','latex')
zlabel('dtheta (rad)','Interpreter','latex')
zlim([0 0.1])
camlight(110,70)
h.EdgeColor = 'none';
h.AmbientStrength = 0.7;
a = gca;
a.TickLabelInterpreter = 'latex';
a.Box = 'on';
a.BoxStyle = 'full';
title_latex = ['$' latex(surfacexy) '$'];
title('L = 0, dra = 10km, dP = 0.01' ,'Interpreter','latex')