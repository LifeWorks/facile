function n= G437_I43_s(a)

if strcmp(a, 'G1137_0x')
	n= 1;
elseif strcmp(a, 'TG0000_1')
	n= 2;
elseif strcmp(a, 'G1137_0x_TG0000_1_i00')
	n= 3;
elseif strcmp(a, 'TG0000_0')
	n= 4;
elseif strcmp(a, 'G3354_x')
	n= 5;
elseif strcmp(a, 'G3354_x_TG0000_0_i00')
	n= 6;
elseif strcmp(a, 'G0166_xxxx')
	n= 7;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_i00')
	n= 8;
elseif strcmp(a, 'G0166_xxxx_G3354_x_i00')
	n= 9;
elseif strcmp(a, 'G3786_x')
	n= 10;
elseif strcmp(a, 'G0166_xxxx_G3786_x_i00')
	n= 11;
elseif strcmp(a, 'LG0000_x')
	n= 12;
elseif strcmp(a, 'G0166_xxxx_LG0000_x_i00')
	n= 13;
elseif strcmp(a, 'G0166_xxxx_LG0000_x_i01')
	n= 14;
elseif strcmp(a, 'G3786_x_LG0000_x_i00')
	n= 15;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_G3354_x_i00')
	n= 16;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_LG0000_x_i00')
	n= 17;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_LG0000_x_i01')
	n= 18;
elseif strcmp(a, 'G0166_xxxx_G3354_x_G3786_x_i00')
	n= 19;
elseif strcmp(a, 'G0166_xxxx_G3354_x_LG0000_x_i00')
	n= 20;
elseif strcmp(a, 'G0166_xxxx_G3354_x_LG0000_x_i01')
	n= 21;
elseif strcmp(a, 'G0166_xxxx_G3786_x_LG0000_x_i00')
	n= 22;
elseif strcmp(a, 'G0166_xxxx_G3786_x_LG0000_x_i01')
	n= 23;
elseif strcmp(a, 'G0166_xxxx_LG0000_x_LG0000_x_i00')
	n= 24;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00')
	n= 25;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01')
	n= 26;
elseif strcmp(a, 'G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00')
	n= 27;
elseif strcmp(a, 'G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00')
	n= 28;
elseif strcmp(a, 'G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01')
	n= 29;
elseif strcmp(a, 'G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00')
	n= 30;
elseif strcmp(a, 'G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00')
	n= 31;
else
	disp('ERROR!');
	n= -1;
end;