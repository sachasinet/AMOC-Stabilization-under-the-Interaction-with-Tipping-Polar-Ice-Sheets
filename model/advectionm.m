function [ad] = advectionm(sol,pars)
% This is the modified (smooth) advection term in the rooth model
% Author: Sacha Sinet, 2021-2022, contact -> s.a.m.sinet@uu.nl
ad = zeros(5,1);

[VR,~,~,~,~,~,Stot]=setvalues(pars);
T1=sol(1);
T2=sol(2);
T3=sol(3);
S1=sol(4);
S3=sol(5);

q = (T3-T1)+(S1-S3);

S2 = (Stot-S1-S3)/VR;
k=200;

ad(1) = q/(1-exp(-k*q))*(T2-T1) + q/(1-exp(k*q))*(T1-T3);
ad(2) = q/(1-exp(-k*q))*(T3-T2)/VR + q/(1-exp(k*q))*(T2-T1)/VR;
ad(3) = q/(1-exp(-k*q))*(T1-T3) + q/(1-exp(k*q))*(T3-T2);
ad(4) = q/(1-exp(-k*q))*(S2-S1) + q/(1-exp(k*q))*(S1-S3);
ad(5) = q/(1-exp(-k*q))*(S1-S3) + q/(1-exp(k*q))*(S3-S2);


end