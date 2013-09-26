function n= steady_state_s(a)

if strcmp(a, 'X')
	n= 1;
elseif strcmp(a, 'Y')
	n= 2;
elseif strcmp(a, 'Z')
	n= 3;
else
	disp('ERROR!');
	n= -1;
end;