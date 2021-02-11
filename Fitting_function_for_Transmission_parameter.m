%SIR MODEL
% use the function ODE45 to calculate differential equations http://uk.mathworks.com/help/matlab/ref/ode45.html
%clear all

%FIRST WE NEED TO RUN THIS ONE TO GET THE TRUE VALUE OF BETA, THEN MOVE TO EXERCISE 1_MAIN:
%USE CURVE FITTING FOR THIS ONE by importing the function described below
%and using months and cases as x and y respectively.
% months=1:1:10;
% cases=[30 35 38 34 32 33 28 25 20 22]

function results = Fit_function(betaC,x)

%to solve differential equations we need to define initial conditions and time interval
y0=[159 0 0 158 0 0 0 0];  %1000 susceptibles, 30 exposed, 20 infected, 0 recovered. Try modifying initial values, what happen to the plots?
tspan=[1 9];
years=1:1:9;
cases=[6.32 12.69 59.639 54.18 65.01 60.54 136.94 108.65 47.52];
birth= 0.000032;
death= 0.000019;
probbact= 0.011;
drugsensmort= 0.2;
drugresimort= 0.6;
treatmentrat1= 0.2;
clearance= 0.1;
recovresist= 0.1;
recovsuscep=0.2;
enhancmentinf= 0.05;
badtreatinfect= 0.0135;
clearinffect= 0.05;
riskcoeff= 1.27;
exposedprop= 0.518;
%define the parameters
%recovery rate, infectious period = 1/gamma
antibioticexp= 0.012;
expo_suscep=1; 
passawaytime= 0.1; %1/10
exptime=1;


%cases=[30 35 38 34 32 33 28 25 20 22]; %remember to run this line in the command window if not in the main file otherwise matlab does not recognise the input cases

[t,y] = ode45(@(t,y) Diss_equations(t,y,birth, death,recovresist, recovsuscep,treatmentrat1, clearance,probbact,exposedprop,riskcoeff,betaC, badtreatinfect, enhancmentinf, drugresimort, drugsensmort,clearinffect,passawaytime,expo_suscep,antibioticexp,exptime),tspan,y0);

results = zeros(size(x, 1), size(x, 2));
%Fill the output matrix
    for r=1:size(results,1)
        for c=1:size(results,2)
           results(r,c) = interp1q(t,y(:,2),x(r, c));
        end
    end
end