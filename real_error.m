syms re ra l theta lamda
syms dre dra dl dtheta
rl = sqrt(re^2+(re+ra)^2-2*re*(re+ra)*cos(l))
phi = theta - asin((re/rl)*sin(l))
P = 2*exp(-0.5*(phi/0.445)^2) %+20*log(lamda/(4*pi*rl))

%the fractional error of power transmitted from satellite (not including propogation to ground)
frac = sqrt((diff(P, l)*dl)^2+(diff(P, ra)*dra)^2+(diff(P, theta)*dtheta)^2)/P


%subs(frac, [re ra l theta dre dra dl dtheta], [6400 800 0.48 0.1 1 1 0.01 0.01])
placed=subs(frac, [lamda re ra l theta], [700E6 6400E3 800E3 0.4759 1]) %worst case numbers
%double(subs(rl, [re ra l], [6400E3 800E3 0.4759]))
double(subs(placed, [dl dra dtheta], [0.01 1000 0.01]))
surffuncdra=solve(placed==0.0048, dra) %looking for 1% fractional error solution
double(subs(surffuncdra(1), [dl dtheta], [0.01 0.2]))
fsurf(surffuncdra(2), [0 0.1 0 0.1])