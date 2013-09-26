function n= config_r(a)

if strcmp(a, 'f1')
	n= 1;
elseif strcmp(a, 'f2')
	n= 2;
else
	disp('ERROR!');
	n= -1;
end;