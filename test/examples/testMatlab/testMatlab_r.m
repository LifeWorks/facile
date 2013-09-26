function n= testMatlab_r(a)

if strcmp(a, 'sourceS0f')
	n= 1;
elseif strcmp(a, 'fa0')
	n= 2;
elseif strcmp(a, 'fb0')
	n= 3;
elseif strcmp(a, 'b1')
	n= 4;
elseif strcmp(a, 'f1')
	n= 5;
elseif strcmp(a, 'f2')
	n= 6;
elseif strcmp(a, 'b2')
	n= 7;
elseif strcmp(a, 'sinkCf')
	n= 8;
elseif strcmp(a, 'sourceXb')
	n= 9;
elseif strcmp(a, 'sourceXf')
	n= 10;
elseif strcmp(a, 'sinkYb')
	n= 11;
elseif strcmp(a, 'sinkYf')
	n= 12;
elseif strcmp(a, 'sinkZb')
	n= 13;
elseif strcmp(a, 'sinkZf')
	n= 14;
elseif strcmp(a, 'gg')
	n= 15;
elseif strcmp(a, 'hh')
	n= 16;
elseif strcmp(a, 'ii')
	n= 17;
elseif strcmp(a, 'm1')
	n= 18;
elseif strcmp(a, 'm2')
	n= 19;
elseif strcmp(a, 'v2')
	n= 20;
else
	disp('ERROR!');
	n= -1;
end;