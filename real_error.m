syms re ra l theta
syms dre dra dl dtheta
rl = sqrt(re^2+(re+ra)^2-2*re*(re+ra)*cos(l))
phi = theta + asin((re/rl)*sin(l))
P = 2*exp(-0.5*(phi/0.4)^2)

%the fractional error of power transmitted from satellite (not including propogation to ground)
frac = sqrt((diff(P, l)*dl)^2+(diff(P, ra)*dra)^2+(diff(P, theta)*dtheta)^2)/P

digits(4)
vpa(subs(frac, [re ra l theta dre dra dl dtheta], [6400 800 0.1 0.1 0.1 0.1 0.1 0.1]))

%the regions where we look for the solution space
slice = 4
re_num = 6400000
ra_num = 800000
theta_num = 0
l_num = 0:(pi-0)/slice:pi
dl = 0:(pi-0)/slice:pi
dtheta = 0:(pi-0)/slice:pi

f1 = subs(frac, [re ra theta dra], [re_num ra_num theta_num 0])
f2 = subs(f1, l, l_num)
f2(2)
fsurf(f2(2)) %why is this not  plotting anything?

for i = 1:slice
   %double(subs(f2, [dl dtheta], []))
end