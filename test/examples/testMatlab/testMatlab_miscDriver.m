% initial values (free nodes only)
X = 1.66057788110262e-08;
ivalues = [X];

% rate constants
f1= 1;
rates= [f1];

% time interval
t0= 0;
tf= 15.0;

% call solver routine 
[t,y]= testMatlab_misc([t0:0.2:tf], ivalues, rates);

% map free node state vector names
X = y(:,1); 





% issue done message for calling/wrapper scripts
disp('Facile driver script done');

