function [t1, y]= testMatlab_misc(tv, initialvalues, rateconstants)

[t1, y]= testMatlab_misc_ode_event(@testMatlab_misc_odes, [1.0 5.0 9.0], tv, initialvalues, odeset('MaxStep', 0.01), rateconstants);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = testMatlab_misc_odes(t, y, rateconstants)

% state vector to node mapping
X = y(1);

% constants
f1 = rateconstants(1);




% differential equations for independent species
dydt(1)= + f1 ;
dydt = dydt(:);

