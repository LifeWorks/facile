function n= poisson_s(a)

if strcmp(a, 'D')
	n= 1;
elseif strcmp(a, 'A')
	n= 2;
else
	disp('ERROR!');
	n= -1;
end;