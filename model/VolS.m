function Vol=VolS(state,parsS)
% Computation of the adimensional volume of the MIS
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl

[NS,ds]=setvalues(parsS);
heights = state(1:NS);
xg = state(2*NS+1);

Vol = 0;
Vol = Vol + heights(1)*ds/2 + heights(end)*ds/2;
for i=2:(length(heights)-1)
    Vol = Vol + heights(i)*ds;
end
Vol=Vol*xg;
Vol=2*Vol; % For symmetry
end