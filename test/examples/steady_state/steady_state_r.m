function n= steady_state_r(a)

if strcmp(a, 'k_sinkx')
	n= 1;
elseif strcmp(a, 'k_sinky')
	n= 2;
elseif strcmp(a, 'k_sinkz')
	n= 3;
else
	disp('ERROR!');
	n= -1;
end;