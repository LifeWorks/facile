% converts simulation output
% assumes simulation time series data is in simdata variable

t= simdata(:,1);
GEF= simdata(:,2);
GTPase_0= simdata(:,3);
GEF_GTPase_0= simdata(:,4);
GTPase_1= simdata(:,5);
GTPase_1b= simdata(:,6);

% concentrations
if size(simdata,2) > 6
    GEFconc= simdata(:,7);
    GTPase_0conc= simdata(:,8);
    GEF_GTPase_0conc= simdata(:,9);
    GTPase_1conc= simdata(:,10);
    GTPase_1bconc= simdata(:,11);
end;

