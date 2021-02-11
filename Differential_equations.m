function dydt = Diss_equations(t,y,birth, death,recovresist, recovsuscep,treatmentrat1, clearance,probbact,exposedprop,riskcoeff,betaC, badtreatinfect, enhancmentinf, drugresimort, drugsensmort,clearinffect,passawaytime,expo_suscep,antibioticexp,exptime)  
 
    dydt1 =  birth*(y(1)+y(2)+y(3)+y(4)+y(5)+y(6))-death*y(1)+recovresist*y(6)+y(5)*recovsuscep+(treatmentrat1+clearance)*y(4)-(1-probbact)*exposedprop*y(1)*expo_suscep-(1-probbact)*(1-exposedprop)*y(1)*expo_suscep-probbact*y(1)*expo_suscep  ;
    dydt2 = -death*y(2)+((1-probbact)*exposedprop*y(1)*expo_suscep)-riskcoeff*(betaC*y(2)*y(4)*antibioticexp*exptime) ;
    dydt3 = -death*y(3)+((1-probbact)*(1-exposedprop)*y(1)*expo_suscep)-(betaC*y(3)*y(4)*antibioticexp*exptime) ;
    dydt4 = -death*y(4)-(treatmentrat1+clearance)*y(4)+riskcoeff*(betaC*y(2)*y(4)*antibioticexp*exptime)+betaC*y(3)*y(4)*antibioticexp*exptime-probbact*(1-clearance)*y(4)+(1-recovresist)*clearinffect*y(6); 
    dydt5 = -death*y(5)+probbact*y(1)*expo_suscep-recovsuscep*y(5)-(1-recovsuscep)*drugsensmort*y(5)*passawaytime+(1-recovresist)*enhancmentinf*y(6)-(1-recovsuscep)*badtreatinfect*y(5);
    dydt6 = -death*y(6)+probbact*(1-clearance)*y(4)-(1-recovresist)*clearinffect*y(6)-recovresist*y(6)-(1-recovresist)*drugresimort*y(6)*passawaytime-(1-recovresist)*enhancmentinf*y(6)+(1-recovsuscep)*badtreatinfect*y(5);
    dydt7 = (1-recovsuscep)*drugsensmort*y(5)*passawaytime;
    dydt8 = (1-recovresist)*drugresimort*y(6)*passawaytime;
 
    dydt = [dydt1; dydt2; dydt3; dydt4; dydt5; dydt6; dydt7; dydt8];
    
end