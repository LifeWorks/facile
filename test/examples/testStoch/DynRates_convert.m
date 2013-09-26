% converts simulation output
% assumes simulation time series data is in simdata variable

t= simdata(:,1);
X= simdata(:,2);
Y= simdata(:,3);
XY= simdata(:,4);
A1= simdata(:,5);
B= simdata(:,6);
A2= simdata(:,7);
A3= simdata(:,8);
A4= simdata(:,9);
D1= simdata(:,10);
E1= simdata(:,11);
D2= simdata(:,12);
E2= simdata(:,13);
D3= simdata(:,14);
E3= simdata(:,15);
D4= simdata(:,16);
E4= simdata(:,17);

% concentrations
if size(simdata,2) > 17
    Xconc= simdata(:,18);
    Yconc= simdata(:,19);
    XYconc= simdata(:,20);
    A1conc= simdata(:,21);
    Bconc= simdata(:,22);
    A2conc= simdata(:,23);
    A3conc= simdata(:,24);
    A4conc= simdata(:,25);
    D1conc= simdata(:,26);
    E1conc= simdata(:,27);
    D2conc= simdata(:,28);
    E2conc= simdata(:,29);
    D3conc= simdata(:,30);
    E3conc= simdata(:,31);
    D4conc= simdata(:,32);
    E4conc= simdata(:,33);
end;

