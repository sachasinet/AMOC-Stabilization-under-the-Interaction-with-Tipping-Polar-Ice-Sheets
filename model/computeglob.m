function dglob=computeglob(dTe)
% To compute the global warming from the equatorial warming
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

% syms dTN dT1 dT dTes
% ampN = 2;
% amp1 = 1.3;
% surftot=4*pi;
% surfN = 2*pi*(2-sqrt(3))/2;
% surf1 = 2*pi*(-1+sqrt(3))/2;
% glob = dTN*surfN/surftot+dT1*surf1/surftot +dTes*(surftot-surf1-surfN)/surftot;
% glob = subs(glob,[dTN,dT1],[ampN*dTes,amp1*dTes]);


dglob = (3*dTe)/4 + (52616922412711097*dTe)/(45035996273704960*pi);

end