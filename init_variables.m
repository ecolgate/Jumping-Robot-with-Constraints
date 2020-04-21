%% init_variables.m
%
% Description:
%   Whereas parameters tend to be defined once and then stay that way, variables change.
%   The most important variables, such as time and state, are handled
%   naturally by ode45.  Other variables, however, need to be passed as
%   function arguments.  Rather than put everything into params, I am
%   choosing to put things that vary (e.g., motor torques, constraint
%   forces, etc) into variables.
%   
% Inputs:
%   params

% Outputs:
%   variables: a struct with many elements

function variables = init_variables(params)
    % motor torques:
        % create a timeline for use here
        tfinal = params.sim.tfinal;
        dt = params.sim.dt;
        % push the spine cw against the spring
        variables.motor.spine.time = 0:dt:tfinal;
        variables.motor.spine.torque = 0*params.sim.motor.spine.peaktorque*ones(1,length(variables.motor.spine.time)); % push the spine cw against the spring
        variables.motor.spine.torque(1) = 0;  % make the initial torque zero
        % push the body upward as hard as possible then brake
        variables.motor.body.time1 = 0:dt:.32;
        variables.motor.body.torque1 = params.sim.motor.body.peaktorque*ones(1,length(variables.motor.body.time1)); 
        variables.motor.body.torque1(1) = params.model.dyn.body.m*params.model.dyn.g*params.model.geom.body.r;  % make the initial torque enough to hold up body
        variables.motor.body.time2 = .32+dt:dt:.35;
        variables.motor.body.torque2 = -8*params.sim.motor.body.peaktorque*ones(1,length(variables.motor.body.time2));
        variables.motor.body.time3 = .35+dt:dt:tfinal;
        variables.motor.body.torque3 = params.model.dyn.body.m*params.model.dyn.g*params.model.geom.body.r*ones(1,length(variables.motor.body.time3));
        variables.motor.body.time = horzcat(variables.motor.body.time1,variables.motor.body.time2,variables.motor.body.time3);
        variables.motor.body.torque = horzcat(variables.motor.body.torque1,variables.motor.body.torque2,variables.motor.body.torque3);
        
    % constraint forces:  1 if active; 0 if inactive
    variables.foot.constraints = [1,1];   % initially, both left and right constraints are active 
    

end