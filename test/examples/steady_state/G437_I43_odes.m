% File generated by Facile version 0.53
%
function dydt = G437_I43_odes(t, y)

global event_flags;
global event_times

% state vector to node mapping
G1137_0x = y(1);
TG0000_1 = y(2);
G1137_0x_TG0000_1_i00 = y(3);
TG0000_0 = y(4);
G3354_x = y(5);
G3354_x_TG0000_0_i00 = y(6);
G0166_xxxx = y(7);
G0166_xxxx_G1137_0x_i00 = y(8);
G0166_xxxx_G3354_x_i00 = y(9);
G3786_x = y(10);
G0166_xxxx_G3786_x_i00 = y(11);
LG0000_x = y(12);
G0166_xxxx_LG0000_x_i00 = y(13);
G0166_xxxx_LG0000_x_i01 = y(14);
G3786_x_LG0000_x_i00 = y(15);
G0166_xxxx_G1137_0x_G3354_x_i00 = y(16);
G0166_xxxx_G1137_0x_LG0000_x_i00 = y(17);
G0166_xxxx_G1137_0x_LG0000_x_i01 = y(18);
G0166_xxxx_G3354_x_G3786_x_i00 = y(19);
G0166_xxxx_G3354_x_LG0000_x_i00 = y(20);
G0166_xxxx_G3354_x_LG0000_x_i01 = y(21);
G0166_xxxx_G3786_x_LG0000_x_i00 = y(22);
G0166_xxxx_G3786_x_LG0000_x_i01 = y(23);
G0166_xxxx_LG0000_x_LG0000_x_i00 = y(24);
G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 = y(25);
G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 = y(26);
G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 = y(27);
G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 = y(28);
G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 = y(29);
G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 = y(30);
G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 = y(31);

% constants and constant expressions
global ode_rate_constants;
fb00 = ode_rate_constants(1);
bb00 = ode_rate_constants(2);
kp00 = ode_rate_constants(3);
fb01 = ode_rate_constants(4);
bb01 = ode_rate_constants(5);
kp01 = ode_rate_constants(6);
fb02 = ode_rate_constants(7);
bb02 = ode_rate_constants(8);
fb03 = ode_rate_constants(9);
bb03 = ode_rate_constants(10);
fb04 = ode_rate_constants(11);
bb04 = ode_rate_constants(12);
bb05 = ode_rate_constants(13);
bb06 = ode_rate_constants(14);
fb05 = ode_rate_constants(15);
clamp_sink_LG0000_x = ode_rate_constants(16);



% dynamic rate expressions
clamp_source_LG0000_x = (+(event_flags(1) && ~event_flags(4))*min((t-event_times(1))/1000, 1)*4*4.0+event_flags(4)*max(1-(t-event_times(4))/1000, 0)*4*4.0+(event_flags(2) && ~event_flags(5))*min((t-event_times(2))/1000, 1)*4*4.0+event_flags(5)*max(1-(t-event_times(5))/1000, 0)*4*4.0+(event_flags(3) && ~event_flags(6))*min((t-event_times(3))/1000, 1)*4*4.0+event_flags(6)*max(1-(t-event_times(6))/1000, 0)*4*4.0);

% differential equations for independent species
dydt(size(y,1),1) = 0;
dydt(1)= + bb00*G1137_0x_TG0000_1_i00 + kp00*G1137_0x_TG0000_1_i00 + bb02*G0166_xxxx_G1137_0x_i00 + bb02*G0166_xxxx_G1137_0x_G3354_x_i00 + bb02*G0166_xxxx_G1137_0x_LG0000_x_i00 + bb02*G0166_xxxx_G1137_0x_LG0000_x_i01 + bb02*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 + bb02*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 + bb02*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 - fb00*G1137_0x*TG0000_1 - fb02*G0166_xxxx*G1137_0x - fb02*G1137_0x*G0166_xxxx_G3354_x_i00 - fb02*G1137_0x*G0166_xxxx_LG0000_x_i00 - fb02*G1137_0x*G0166_xxxx_LG0000_x_i01 - fb02*G1137_0x*G0166_xxxx_G3354_x_LG0000_x_i00 - fb02*G1137_0x*G0166_xxxx_G3354_x_LG0000_x_i01 - fb02*G1137_0x*G0166_xxxx_LG0000_x_LG0000_x_i00 ;
dydt(2)= + bb00*G1137_0x_TG0000_1_i00 + kp01*G3354_x_TG0000_0_i00 - fb00*G1137_0x*TG0000_1 ;
dydt(3)= + fb00*G1137_0x*TG0000_1 - bb00*G1137_0x_TG0000_1_i00 - kp00*G1137_0x_TG0000_1_i00 ;
dydt(4)= + kp00*G1137_0x_TG0000_1_i00 + bb01*G3354_x_TG0000_0_i00 - fb01*G3354_x*TG0000_0 ;
dydt(5)= + bb01*G3354_x_TG0000_0_i00 + kp01*G3354_x_TG0000_0_i00 + bb03*G0166_xxxx_G3354_x_i00 + bb03*G0166_xxxx_G1137_0x_G3354_x_i00 + bb03*G0166_xxxx_G3354_x_G3786_x_i00 + bb03*G0166_xxxx_G3354_x_LG0000_x_i00 + bb03*G0166_xxxx_G3354_x_LG0000_x_i01 + bb03*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 + bb03*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 + bb03*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 + bb03*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 + bb03*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 - fb01*G3354_x*TG0000_0 - fb03*G0166_xxxx*G3354_x - fb03*G3354_x*G0166_xxxx_G1137_0x_i00 - fb03*G3354_x*G0166_xxxx_G3786_x_i00 - fb03*G3354_x*G0166_xxxx_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_LG0000_x_i01 - fb03*G3354_x*G0166_xxxx_G1137_0x_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_G1137_0x_LG0000_x_i01 - fb03*G3354_x*G0166_xxxx_G3786_x_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_G3786_x_LG0000_x_i01 - fb03*G3354_x*G0166_xxxx_LG0000_x_LG0000_x_i00 ;
dydt(6)= + fb01*G3354_x*TG0000_0 - bb01*G3354_x_TG0000_0_i00 - kp01*G3354_x_TG0000_0_i00 ;
dydt(7)= + bb02*G0166_xxxx_G1137_0x_i00 + bb03*G0166_xxxx_G3354_x_i00 + bb04*G0166_xxxx_G3786_x_i00 + bb05*G0166_xxxx_LG0000_x_i00 + bb06*G0166_xxxx_LG0000_x_i01 - fb02*G0166_xxxx*G1137_0x - fb03*G0166_xxxx*G3354_x - fb04*G0166_xxxx*G3786_x - fb04*G0166_xxxx*LG0000_x - fb01*G0166_xxxx*LG0000_x ;
dydt(8)= + fb02*G0166_xxxx*G1137_0x + bb03*G0166_xxxx_G1137_0x_G3354_x_i00 + bb05*G0166_xxxx_G1137_0x_LG0000_x_i00 + bb06*G0166_xxxx_G1137_0x_LG0000_x_i01 - bb02*G0166_xxxx_G1137_0x_i00 - fb03*G3354_x*G0166_xxxx_G1137_0x_i00 - fb04*LG0000_x*G0166_xxxx_G1137_0x_i00 - fb01*LG0000_x*G0166_xxxx_G1137_0x_i00 ;
dydt(9)= + fb03*G0166_xxxx*G3354_x + bb02*G0166_xxxx_G1137_0x_G3354_x_i00 + bb04*G0166_xxxx_G3354_x_G3786_x_i00 + bb05*G0166_xxxx_G3354_x_LG0000_x_i00 + bb06*G0166_xxxx_G3354_x_LG0000_x_i01 - bb03*G0166_xxxx_G3354_x_i00 - fb02*G1137_0x*G0166_xxxx_G3354_x_i00 - fb04*G3786_x*G0166_xxxx_G3354_x_i00 - fb04*LG0000_x*G0166_xxxx_G3354_x_i00 - fb01*LG0000_x*G0166_xxxx_G3354_x_i00 ;
dydt(10)= + bb04*G0166_xxxx_G3786_x_i00 + bb01*G3786_x_LG0000_x_i00 + bb04*G0166_xxxx_G3354_x_G3786_x_i00 + bb04*G0166_xxxx_G3786_x_LG0000_x_i00 + bb04*G0166_xxxx_G3786_x_LG0000_x_i01 + bb04*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 + bb04*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 + bb04*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 - fb04*G0166_xxxx*G3786_x - fb05*G3786_x*LG0000_x - fb04*G3786_x*G0166_xxxx_G3354_x_i00 - fb04*G3786_x*G0166_xxxx_LG0000_x_i00 - fb04*G3786_x*G0166_xxxx_LG0000_x_i01 - fb04*G3786_x*G0166_xxxx_G3354_x_LG0000_x_i00 - fb04*G3786_x*G0166_xxxx_G3354_x_LG0000_x_i01 - fb04*G3786_x*G0166_xxxx_LG0000_x_LG0000_x_i00 ;
dydt(11)= + fb04*G0166_xxxx*G3786_x + bb03*G0166_xxxx_G3354_x_G3786_x_i00 + bb05*G0166_xxxx_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_G3786_x_LG0000_x_i01 - bb04*G0166_xxxx_G3786_x_i00 - fb03*G3354_x*G0166_xxxx_G3786_x_i00 - fb04*LG0000_x*G0166_xxxx_G3786_x_i00 - fb01*LG0000_x*G0166_xxxx_G3786_x_i00 ;
dydt(12)= + bb05*G0166_xxxx_LG0000_x_i00 + bb06*G0166_xxxx_LG0000_x_i01 + bb01*G3786_x_LG0000_x_i00 + bb05*G0166_xxxx_G1137_0x_LG0000_x_i00 + bb06*G0166_xxxx_G1137_0x_LG0000_x_i01 + bb05*G0166_xxxx_G3354_x_LG0000_x_i00 + bb06*G0166_xxxx_G3354_x_LG0000_x_i01 + bb05*G0166_xxxx_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_G3786_x_LG0000_x_i01 + bb06*G0166_xxxx_LG0000_x_LG0000_x_i00 + bb05*G0166_xxxx_LG0000_x_LG0000_x_i00 + bb05*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 + bb06*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 + bb06*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 + bb05*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 + bb05*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 + bb06*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 + bb05*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 + bb06*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 + bb05*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 + clamp_source_LG0000_x - fb04*G0166_xxxx*LG0000_x - fb01*G0166_xxxx*LG0000_x - fb05*G3786_x*LG0000_x - fb04*LG0000_x*G0166_xxxx_G1137_0x_i00 - fb01*LG0000_x*G0166_xxxx_G1137_0x_i00 - fb04*LG0000_x*G0166_xxxx_G3354_x_i00 - fb01*LG0000_x*G0166_xxxx_G3354_x_i00 - fb04*LG0000_x*G0166_xxxx_G3786_x_i00 - fb01*LG0000_x*G0166_xxxx_G3786_x_i00 - fb01*LG0000_x*G0166_xxxx_LG0000_x_i00 - fb04*LG0000_x*G0166_xxxx_LG0000_x_i01 - fb04*LG0000_x*G0166_xxxx_G1137_0x_G3354_x_i00 - fb01*LG0000_x*G0166_xxxx_G1137_0x_G3354_x_i00 - fb01*LG0000_x*G0166_xxxx_G1137_0x_LG0000_x_i00 - fb04*LG0000_x*G0166_xxxx_G1137_0x_LG0000_x_i01 - fb04*LG0000_x*G0166_xxxx_G3354_x_G3786_x_i00 - fb01*LG0000_x*G0166_xxxx_G3354_x_G3786_x_i00 - fb01*LG0000_x*G0166_xxxx_G3354_x_LG0000_x_i00 - fb04*LG0000_x*G0166_xxxx_G3354_x_LG0000_x_i01 - fb01*LG0000_x*G0166_xxxx_G3786_x_LG0000_x_i00 - fb04*LG0000_x*G0166_xxxx_G3786_x_LG0000_x_i01 - clamp_sink_LG0000_x*LG0000_x ;
dydt(13)= + fb04*G0166_xxxx*LG0000_x + bb02*G0166_xxxx_G1137_0x_LG0000_x_i00 + bb03*G0166_xxxx_G3354_x_LG0000_x_i00 + bb04*G0166_xxxx_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_LG0000_x_LG0000_x_i00 - bb05*G0166_xxxx_LG0000_x_i00 - fb02*G1137_0x*G0166_xxxx_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_LG0000_x_i00 - fb04*G3786_x*G0166_xxxx_LG0000_x_i00 - fb01*LG0000_x*G0166_xxxx_LG0000_x_i00 ;
dydt(14)= + fb01*G0166_xxxx*LG0000_x + bb02*G0166_xxxx_G1137_0x_LG0000_x_i01 + bb03*G0166_xxxx_G3354_x_LG0000_x_i01 + bb04*G0166_xxxx_G3786_x_LG0000_x_i01 + bb05*G0166_xxxx_LG0000_x_LG0000_x_i00 - bb06*G0166_xxxx_LG0000_x_i01 - fb02*G1137_0x*G0166_xxxx_LG0000_x_i01 - fb03*G3354_x*G0166_xxxx_LG0000_x_i01 - fb04*G3786_x*G0166_xxxx_LG0000_x_i01 - fb04*LG0000_x*G0166_xxxx_LG0000_x_i01 ;
dydt(15)= + fb05*G3786_x*LG0000_x - bb01*G3786_x_LG0000_x_i00 ;
dydt(16)= + fb02*G1137_0x*G0166_xxxx_G3354_x_i00 + fb03*G3354_x*G0166_xxxx_G1137_0x_i00 + bb05*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 + bb06*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 - bb02*G0166_xxxx_G1137_0x_G3354_x_i00 - bb03*G0166_xxxx_G1137_0x_G3354_x_i00 - fb04*LG0000_x*G0166_xxxx_G1137_0x_G3354_x_i00 - fb01*LG0000_x*G0166_xxxx_G1137_0x_G3354_x_i00 ;
dydt(17)= + fb02*G1137_0x*G0166_xxxx_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G1137_0x_i00 + bb03*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 + bb06*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 - bb02*G0166_xxxx_G1137_0x_LG0000_x_i00 - bb05*G0166_xxxx_G1137_0x_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_G1137_0x_LG0000_x_i00 - fb01*LG0000_x*G0166_xxxx_G1137_0x_LG0000_x_i00 ;
dydt(18)= + fb02*G1137_0x*G0166_xxxx_LG0000_x_i01 + fb01*LG0000_x*G0166_xxxx_G1137_0x_i00 + bb03*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 + bb05*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 - bb02*G0166_xxxx_G1137_0x_LG0000_x_i01 - bb06*G0166_xxxx_G1137_0x_LG0000_x_i01 - fb03*G3354_x*G0166_xxxx_G1137_0x_LG0000_x_i01 - fb04*LG0000_x*G0166_xxxx_G1137_0x_LG0000_x_i01 ;
dydt(19)= + fb03*G3354_x*G0166_xxxx_G3786_x_i00 + fb04*G3786_x*G0166_xxxx_G3354_x_i00 + bb05*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 - bb03*G0166_xxxx_G3354_x_G3786_x_i00 - bb04*G0166_xxxx_G3354_x_G3786_x_i00 - fb04*LG0000_x*G0166_xxxx_G3354_x_G3786_x_i00 - fb01*LG0000_x*G0166_xxxx_G3354_x_G3786_x_i00 ;
dydt(20)= + fb03*G3354_x*G0166_xxxx_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G3354_x_i00 + bb02*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 + bb04*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 - bb03*G0166_xxxx_G3354_x_LG0000_x_i00 - bb05*G0166_xxxx_G3354_x_LG0000_x_i00 - fb02*G1137_0x*G0166_xxxx_G3354_x_LG0000_x_i00 - fb04*G3786_x*G0166_xxxx_G3354_x_LG0000_x_i00 - fb01*LG0000_x*G0166_xxxx_G3354_x_LG0000_x_i00 ;
dydt(21)= + fb03*G3354_x*G0166_xxxx_LG0000_x_i01 + fb01*LG0000_x*G0166_xxxx_G3354_x_i00 + bb02*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 + bb04*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 + bb05*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 - bb03*G0166_xxxx_G3354_x_LG0000_x_i01 - bb06*G0166_xxxx_G3354_x_LG0000_x_i01 - fb02*G1137_0x*G0166_xxxx_G3354_x_LG0000_x_i01 - fb04*G3786_x*G0166_xxxx_G3354_x_LG0000_x_i01 - fb04*LG0000_x*G0166_xxxx_G3354_x_LG0000_x_i01 ;
dydt(22)= + fb04*G3786_x*G0166_xxxx_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G3786_x_i00 + bb03*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 + bb06*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 - bb04*G0166_xxxx_G3786_x_LG0000_x_i00 - bb05*G0166_xxxx_G3786_x_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_G3786_x_LG0000_x_i00 - fb01*LG0000_x*G0166_xxxx_G3786_x_LG0000_x_i00 ;
dydt(23)= + fb04*G3786_x*G0166_xxxx_LG0000_x_i01 + fb01*LG0000_x*G0166_xxxx_G3786_x_i00 + bb03*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 + bb05*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 - bb04*G0166_xxxx_G3786_x_LG0000_x_i01 - bb06*G0166_xxxx_G3786_x_LG0000_x_i01 - fb03*G3354_x*G0166_xxxx_G3786_x_LG0000_x_i01 - fb04*LG0000_x*G0166_xxxx_G3786_x_LG0000_x_i01 ;
dydt(24)= + fb01*LG0000_x*G0166_xxxx_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_LG0000_x_i01 + bb02*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 + bb03*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 + bb04*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 - bb06*G0166_xxxx_LG0000_x_LG0000_x_i00 - bb05*G0166_xxxx_LG0000_x_LG0000_x_i00 - fb02*G1137_0x*G0166_xxxx_LG0000_x_LG0000_x_i00 - fb03*G3354_x*G0166_xxxx_LG0000_x_LG0000_x_i00 - fb04*G3786_x*G0166_xxxx_LG0000_x_LG0000_x_i00 ;
dydt(25)= + fb02*G1137_0x*G0166_xxxx_G3354_x_LG0000_x_i00 + fb03*G3354_x*G0166_xxxx_G1137_0x_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G1137_0x_G3354_x_i00 - bb02*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 - bb03*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 - bb05*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i00 ;
dydt(26)= + fb02*G1137_0x*G0166_xxxx_G3354_x_LG0000_x_i01 + fb03*G3354_x*G0166_xxxx_G1137_0x_LG0000_x_i01 + fb01*LG0000_x*G0166_xxxx_G1137_0x_G3354_x_i00 - bb02*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 - bb03*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 - bb06*G0166_xxxx_G1137_0x_G3354_x_LG0000_x_i01 ;
dydt(27)= + fb02*G1137_0x*G0166_xxxx_LG0000_x_LG0000_x_i00 + fb01*LG0000_x*G0166_xxxx_G1137_0x_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G1137_0x_LG0000_x_i01 - bb02*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 - bb06*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 - bb05*G0166_xxxx_G1137_0x_LG0000_x_LG0000_x_i00 ;
dydt(28)= + fb03*G3354_x*G0166_xxxx_G3786_x_LG0000_x_i00 + fb04*G3786_x*G0166_xxxx_G3354_x_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G3354_x_G3786_x_i00 - bb03*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 - bb04*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 - bb05*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i00 ;
dydt(29)= + fb03*G3354_x*G0166_xxxx_G3786_x_LG0000_x_i01 + fb04*G3786_x*G0166_xxxx_G3354_x_LG0000_x_i01 + fb01*LG0000_x*G0166_xxxx_G3354_x_G3786_x_i00 - bb03*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 - bb04*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 - bb06*G0166_xxxx_G3354_x_G3786_x_LG0000_x_i01 ;
dydt(30)= + fb03*G3354_x*G0166_xxxx_LG0000_x_LG0000_x_i00 + fb01*LG0000_x*G0166_xxxx_G3354_x_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G3354_x_LG0000_x_i01 - bb03*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 - bb06*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 - bb05*G0166_xxxx_G3354_x_LG0000_x_LG0000_x_i00 ;
dydt(31)= + fb04*G3786_x*G0166_xxxx_LG0000_x_LG0000_x_i00 + fb01*LG0000_x*G0166_xxxx_G3786_x_LG0000_x_i00 + fb04*LG0000_x*G0166_xxxx_G3786_x_LG0000_x_i01 - bb04*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 - bb06*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 - bb05*G0166_xxxx_G3786_x_LG0000_x_LG0000_x_i00 ;


