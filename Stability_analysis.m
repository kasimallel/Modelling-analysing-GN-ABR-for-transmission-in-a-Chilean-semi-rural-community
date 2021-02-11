% STABILITY SIR WITH BIRTHS AND DEATHS
% Analyse the stability of a system 
clear all;
%Solve a system of equations with MATLAB   https://uk.mathworks.com/help/symbolic/solve.html#inputarg_vars
syms S E U Z I1 I2 D1 D2 birth death recovresist recovsuscep treatmentrat1  clearance probbact exposedprop riskcoeff betaC  badtreatinfect  enhancmentinf drugresimort drugsensmort clearinffect passawaytime expo_suscep antibioticexp exptime                                         % define S, I and R as symbolic https://uk.mathworks.com/help/symbolic/syms.html
eqns = [birth*(S+E+U+Z+I1+I2)-death*S+recovresist*I2+I1*recovsuscep+(treatmentrat1+clearance)*Z-(1-probbact)*exposedprop*S*expo_suscep-(1-probbact)*(1-exposedprop)*S*expo_suscep-probbact*S*expo_suscep==0, -death*E+((1-probbact)*exposedprop*S*expo_suscep)-riskcoeff*(betaC*E*S*antibioticexp*exptime)==0, -death*U+((1-probbact)*(1-exposedprop)*S*expo_suscep)-(betaC*U*Z*antibioticexp*exptime)==0, -death*Z-(treatmentrat1+clearance)*Z+riskcoeff*(betaC*E*Z*antibioticexp*exptime)+betaC*U*Z*antibioticexp*exptime-probbact*(1-clearance)*Z+(1-recovresist)*clearinffect*I2==0, -death*I1+probbact*S*expo_suscep-recovsuscep*I1-(1-recovsuscep)*drugsensmort*I1*passawaytime+(1-recovresist)*enhancmentinf*I2-(1-recovsuscep)*badtreatinfect*I1==0, -death*I2+probbact*(1-clearance)*Z-(1-recovresist)*clearinffect*I2-recovresist*I2-(1-recovresist)*drugresimort*I2*passawaytime-(1-recovresist)*enhancmentinf*I2+(1-recovsuscep)*badtreatinfect*I1==0, (1-recovsuscep)*drugsensmort*I1*passawaytime==0, (1-recovresist)*drugresimort*I2*passawaytime==0]; % define the equations
Sol=solve(eqns, [S E U Z I1 I2 D1 D2]);                              % solve the equations according to variables S, I and R i.e. calculate the steady state solutions

%Look at your solutions writing in the command window Sol.I 

%%
J = jacobian([birth*(S+E+U+Z+I1+I2)-death*S+recovresist*I2+I1*recovsuscep+(treatmentrat1+clearance)*Z-(1-probbact)*exposedprop*S*expo_suscep-(1-probbact)*(1-exposedprop)*S*expo_suscep-probbact*S*expo_suscep, -death*E+((1-probbact)*exposedprop*S*expo_suscep)-riskcoeff*(betaC*E*Z*antibioticexp*exptime), -death*U+((1-probbact)*(1-exposedprop)*S*expo_suscep)-(betaC*U*Z*antibioticexp*exptime),-death*Z-(treatmentrat1+clearance)*Z+riskcoeff*(betaC*E*Z*antibioticexp*exptime)+betaC*U*Z*antibioticexp*exptime-probbact*(1-clearance)*Z+(1-recovresist)*clearinffect*I2,-death*I1+probbact*S*expo_suscep-recovsuscep*I1-(1-recovsuscep)*drugsensmort*I1*passawaytime+(1-recovresist)*enhancmentinf*I2-(1-recovsuscep)*badtreatinfect*I1,-death*I2+probbact*(1-clearance)*Z-(1-recovresist)*clearinffect*I2-recovresist*I2-(1-recovresist)*drugresimort*I2*passawaytime-(1-recovresist)*enhancmentinf*I2+(1-recovsuscep)*badtreatinfect*I1,(1-recovsuscep)*drugsensmort*I1*passawaytime,(1-recovresist)*drugresimort*I2*passawaytime], [S, E, U, Z, I1, I2, D1, D2]);  % calculate the Jacobian matrix  https://uk.mathworks.com/help/symbolic/jacobian.html
JF1 = subs(J,[S E U Z I1 I2 D1 D2],[Sol.S(1) Sol.E(1) Sol.U(1) Sol.Z(1) Sol.I1(1) Sol.I2(1) Sol.D1(1) Sol.D2(1)]);  % evaluate the Jacobian matrix in the dfe https://uk.mathworks.com/help/symbolic/subs.html
e = eig(JF1);                                        % calculate eigenvalues of matrix J in the dfe https://uk.mathworks.com/help/symbolic/jacobian.html





% You can see that the three eigenvalues are the same we found manually in class.S = vpasolve(z1^8 + (8019929151628087*z1^7)/4503599627370496 + (8817334758532073*z1^6)/9007199254740992 + (3891456001138161*z1^5)/18014398509481984 + (4822875502299711*z1^4)/288230376151711744 + (3808955651380341*z1^3)/37778931862957161709568 - (5104528027707415*z1^2)/1237940039285380274899124224, z1)

% The third is the largest eigenvalue which represents R0.
% The dfe is stable when this number is lower than 1. In this case, the
% disease will die out.
%badtreatinfect - birth + clearance + clearinffect + 6*death + enhancmentinf + expo_suscep + probbact + recovsuscep + recovresist + treatmentrat1 - clearance*probbact - badtreatinfect*recovsuscep + drugresimort*passawaytime + drugsensmort*passawaytime - clearinffect*recovresist - enhancmentinf*recovresist - drugsensmort*passawaytime*recovsuscep - drugresimort*passawaytime*recovresist
%You can study the endemic equilibrium by repeating these same steps for
%the first solution of the system i.e. (Sol.S(1),Sol.I(1),Sol.R(1))