function [out] = VolN (state,pars)
% Computation of the adimensional volume of Greenland (without theta function)
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

[NN,drN] = setvalues(pars);


rN   =  -1:drN:1;
rN = rN(length(rN)/2+1:end-1);
rNend = 1; 

out = 0;
out = out + 1/2*(state(1)+state(1))*(rN(1)^2-0^2);
for j = 1:NN-1
    out = out + 1/2*(state(j)+state(j+1))*(rN(j+1)^2-rN(j)^2);
end
out = out + 1/2*(state(NN)+0)*(rNend^2-rN(NN)^2);
out = pi*out;

end