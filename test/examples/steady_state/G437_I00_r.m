function n= G437_I00_r(a)

if strcmp(a, 'fb00')
	n= 1;
elseif strcmp(a, 'bb00')
	n= 2;
elseif strcmp(a, 'kp00')
	n= 3;
elseif strcmp(a, 'fb01')
	n= 4;
elseif strcmp(a, 'bb01')
	n= 5;
elseif strcmp(a, 'kp01')
	n= 6;
elseif strcmp(a, 'fb02')
	n= 7;
elseif strcmp(a, 'bb02')
	n= 8;
elseif strcmp(a, 'fb03')
	n= 9;
elseif strcmp(a, 'bb03')
	n= 10;
elseif strcmp(a, 'fb04')
	n= 11;
elseif strcmp(a, 'bb04')
	n= 12;
elseif strcmp(a, 'bb05')
	n= 13;
elseif strcmp(a, 'bb06')
	n= 14;
elseif strcmp(a, 'fb05')
	n= 15;
elseif strcmp(a, 'clamp_sink_LG0000_x')
	n= 16;
else
	disp('ERROR!');
	n= -1;
end;