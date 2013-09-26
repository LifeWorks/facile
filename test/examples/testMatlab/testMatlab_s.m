function n= testMatlab_s(a)

if strcmp(a, 'S0')
	n= 1;
elseif strcmp(a, 'S1')
	n= 2;
elseif strcmp(a, 'A')
	n= 3;
elseif strcmp(a, 'B')
	n= 4;
elseif strcmp(a, 'C')
	n= 5;
elseif strcmp(a, 'D')
	n= 6;
elseif strcmp(a, 'X')
	n= 7;
elseif strcmp(a, 'Y')
	n= 8;
elseif strcmp(a, 'W')
	n= 9;
elseif strcmp(a, 'Z')
	n= 10;
elseif strcmp(a, 'E')
	n= 11;
elseif strcmp(a, 'F')
	n= 12;
elseif strcmp(a, 'G')
	n= 13;
elseif strcmp(a, 'M1')
	n= 14;
elseif strcmp(a, 'M2')
	n= 15;
else
	disp('ERROR!');
	n= -1;
end;