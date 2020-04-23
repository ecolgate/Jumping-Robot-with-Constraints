%% main.m
%
% Description:
%   Application entry point.
%
% Inputs: none
%
% Outputs: none
%
% Notes:

function main

%% Initialize environment
clear all; % a little inefficient but MATLAB's memory management is buggy.
close all;
clc;

init_env();

%% Initialize parameters
params = init_params;
variables = init_variables(params);

%% Set up events using odeset
options = odeset('Events',@(t,x) robot_events_nested(t,x,variables,params));

%% Simulate the robot forward in time 
x_IC = params.x_IC';    % initial conditions
tnow = 0.0;             % starting time

% start with null matrices for holding results -- we'll be adding to these
% with each segment of the simulation
tsim = [];
xsim = [];
F_list = [];

while tnow < params.sim.tfinal

    tspan = [tnow params.sim.tfinal];
    [tseg, xseg, te, xe, ie] = ode45(@(t,x) ...
                                     robot_dynamics(t,x,variables,params),...
                                     tspan, x_IC, options);

    % augment tsim and xsim; renew ICs
    tsim = [tsim;tseg];
    xsim = [xsim;xseg];
    tnow = tsim(end);
    x_IC = xsim(end,:); 
    
    % compute the constraint forces (why don't I use the
    % robot_events_nested function which already computes the forces?  I
    % could, but I believe that function is used to find zero crossings and
    % therefore is called at times that actually violate the constraints. I
    % could also use an output function in ode45, but I think that amounts
    % to the same thing as doing it here.)
    [Fseg] = constraint_forces(tseg,xseg',variables,params);
    F_list = [F_list,Fseg];
    
    % if simulation terminated before tfinal, determine which constaints
    % are still active, then continue integration
    if tseg(end) < params.sim.tfinal  % termination was triggered by an event
        nc = variables.foot.constraints(1) + 2*variables.foot.constraints(2);  % number of active contacts
        switch nc
            case 1  % the left foot was on the ground prior to termination
                %[F_list] = constraint_forces(tsim,xsim,variables,params);
                variables.foot.constraints(1) = 0;  % now the left foot is off
            case 2 % the right foot only was on the ground prior to termination
                 %[F_list] = constraint_forces(tsim,xsim,variables,params);
                 variables.foot.constraints(2) = 0;  % now the right foot is off
            case 3 % both feet were on the ground prior to termination
                 %[F_list] = constraint_forces(tsim,xsim,variables,params);
                 if ie == 1
                    variables.foot.constraints(1) = 0;  % now the left foot is off
                 else
                    variables.foot.constraints(2) = 0;  % now the right foot is off
                 end
        end
    end
end

figure;
weight = (params.model.dyn.foot.m+params.model.dyn.spine.m+params.model.dyn.body.m)*params.model.dyn.g*ones(1,length(tsim));
plot(tsim,F_list(2,:)+F_list(4,:),'b-','LineWidth',2);
hold on
plot(tsim,weight,'r-','LineWidth',1);
ylabel('Ground Reaction vs Weight (N)')
xlabel('time (sec)')
hold off

%% Un-comment if you want to animate the "jump"
pause;

% Let's resample the simulator output so we can animate with evenly-spaced
% points in (time,state).
% 1) deal with possible duplicate times in tsim:
% (https://www.mathworks.com/matlabcentral/answers/321603-how-do-i-interpolate-1d-data-if-i-do-not-have-unique-values
tsim = cumsum(ones(size(tsim)))*eps + tsim;

% 2) resample the duplicate-free time vector:
t_anim = 0:params.viz.dt:tsim(end);

% 3) resample the state-vs-time array:
x_anim = interp1(tsim,xsim,t_anim);
x_anim = x_anim'; % transpose so that xsim is 10xN (N = number of timesteps)

% 4) resample the constraint forces-vs-time array:
F_anim = interp1(tsim,F_list',t_anim);
F_anim = F_anim';

animate_robot(x_anim(1:5,:),F_anim,params,'trace_foot_com',true,...
    'trace_body_com',true,'trace_spine_tip',true,'show_constraint_forces',true,'video',true);
fprintf('Done!\n');

end