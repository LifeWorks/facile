function [t, y, l] = ode_event(odefun, events, tv, y0, opt, varargin)

% ODE_EVENT(@ODEFUN, EVENTS, Y0, OPT, VARARGIN) calls ode23s to solve the 
% ode file specified by odefun. The ode23s solver often skips over events 
% which happen too quickly due to its procedure for choosing step sizes to
% minimize computational power required. 
%
% The inputs:
%
% ODEFUN - an ode file, just as would be input to the odesolver itself.
% 
% EVENTS - a vector of times at which quick changes in the conditions of
% the ode file may be missed.  Do not supply start and end times as they
% are provided in TV below.
%       ex. if   events = [event1 event2 ... eventn];
%           and      tv = [t0 tf]
%           then the solver will be called over the intervals
%               [t0     event1]
%               [event1 event2]
%                  ...
%               [eventn tf]
%
% TV  - a vector of times at which the result is to be sampled,
%       usually just [t0 tf], but can be e.g. [t0:0.01:tf]
%
% Y0 - a vector of initial conditions for the ode file at time t_start.
%
% OPT - options to pass to the odesolver. If no options are desired, the
% input should be opt = [].
%
% VARARGIN - any other variables to pass to the ode file, such as parameter
% values.
%
% The ouputs:
%
% [t, y] - the solution to the ode file on the interval [t_start t_end].
%
% l - a vector containing the length of each segment solved by ode23s, ie.
% the number of points generated in [t, y] between events. This is useful
% for working with only certain parts of the solution. For example, you may
% give a model time to come to equilibrium before perturbing the system and
% throw that time away as a transient. 

% Initialize the ouputs:
t = [];
y = [];

if (events(1) <= tv(1))
  disp('ERROR: ode_event -- events must be within specified integration vector, compare t0 and first event');
  return
end
if (events(end) >= tv(end))
  disp('ERROR: ode_event -- events must be within specified integration vector, compare tf and last event');
  return
end

events=[tv(1) events tv(end)];

% Call the solver at each event time specified.
for i = 1:length(events)-1
    str=sprintf('ode_event: integrating from %f to %f', events(i), events(i+1));
    disp(str);
    TV = [events(i) tv(find(tv > events(i) & tv < events(i+1))) events(i+1)];
% uncomment this line to debug integration vector
%    str=[sprintf('ode_event: current integration vector is ') sprintf('%.2f ', TV)];
    disp(str);
    [T, Y] = ode45(odefun, TV, y0, opt, varargin{:});
    y0 = Y(end,:);
    y((length(t) + 1):(length(t) + length(T)), 1:length(y0)) = Y;
    t((length(t) + 1):(length(t) + length(T)), 1) = T;
    l(i) = length(T);
end

