syms re ra l theta lamda
syms dre dra dl dtheta
rl = sqrt(re^2+(re+ra)^2-2*re*(re+ra)*cos(l))
phi = theta - asin((re/rl)*sin(l))
P = 2*exp(-0.5*(phi/0.445)^2) %+20*log(lamda/(4*pi*rl))

%the fractional error of power transmitted from satellite (not including propogation to ground)
frac = sqrt((diff(P, l)*dl)^2+(diff(P, ra)*dra)^2+(diff(P, theta)*dtheta)^2)/P


%subs(frac, [re ra l theta dre dra dl dtheta], [6400 800 0.48 0.1 1 1 0.01 0.01])
placed=subs(frac, [lamda re ra l theta], [700E6 6400E3 800E3 0.4759 1]) %worst case numbers

double(subs(placed, [dl dra dtheta], [0.01 1000 0.01]))%just checking the value looks good


surffuncdra=solve(placed==0.0047, dtheta) %looking for 1% fractional error solution, only works when isolating dtheta

double(subs(surffuncdra, [dra dl], [1000 0.01]))%just checking the value looks good

syms x y
surfacexy = subs(surffuncdra(2), [dra dl], [x y]) %converted to xy so that fsurf can plot it
fsurf(surfacexy, [0 10E3 0 0.1])
xlabel('dra')
ylabel('dl')
zlabel('dtheta')