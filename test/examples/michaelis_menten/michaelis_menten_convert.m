% converts simulation output
% assumes simulation time series data is in simdata variable

t= simdata(:,1);
E= simdata(:,2);
S= simdata(:,3);
C= simdata(:,4);
P= simdata(:,5);

% concentrations
if size(simdata,2) > 5
    Econc= simdata(:,6);
    Sconc= simdata(:,7);
    Cconc= simdata(:,8);
    Pconc= simdata(:,9);
end;

