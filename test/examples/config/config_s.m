function n= config_s(a)

if strcmp(a, 'A')
	n= 1;
elseif strcmp(a, 'B')
	n= 2;
elseif strcmp(a, 'X')
	n= 3;
elseif strcmp(a, 'Y')
	n= 4;
else
	disp('ERROR!');
	n= -1;
end;