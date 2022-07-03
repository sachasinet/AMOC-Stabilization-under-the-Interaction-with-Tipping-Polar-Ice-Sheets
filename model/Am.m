function [A,der] = Am(T_ocean)
%Defines the coupling from T3 to the WAIS
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
global T0 T_oceaneq cAm A0S

T_ocean = T_ocean*T0;

A = A0S + cAm.*A0S./T_oceaneq.*(T_ocean-T_oceaneq);
der = cAm.*A0S./T_oceaneq;
end