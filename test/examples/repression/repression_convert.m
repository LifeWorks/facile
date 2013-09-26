% converts simulation output
% assumes simulation time series data is in simdata variable

t= simdata(:,1);
D= simdata(:,2);
C= simdata(:,3);
T= simdata(:,4);
mR= simdata(:,5);
mC1= simdata(:,6);
mC= simdata(:,7);
mC2= simdata(:,8);
mT= simdata(:,9);
A= simdata(:,10);
R= simdata(:,11);
Z= simdata(:,12);

% concentrations
if size(simdata,2) > 12
    Dconc= simdata(:,13);
    Cconc= simdata(:,14);
    Tconc= simdata(:,15);
    mRconc= simdata(:,16);
    mC1conc= simdata(:,17);
    mCconc= simdata(:,18);
    mC2conc= simdata(:,19);
    mTconc= simdata(:,20);
    Aconc= simdata(:,21);
    Rconc= simdata(:,22);
    Zconc= simdata(:,23);
end;

