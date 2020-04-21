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
clear;
close all;
clc;

init_env();

%% Initialize parameters
params = init_params;
variables = init_variables(params);

%% Set up events using odeset
options = odeset('Events',@robot_events_nested);

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
    [tseg, xseg, te, xe, ie] = ode45(@(t,x) robot_dynamics(t,x,variables,params), tspan, x_IC, options);

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

xsim = xsim'; % transpose so that xsim is 10xN (N = number of timesteps)


figure;
weight = (params.model.dyn.foot.m+params.model.dyn.spine.m+params.model.dyn.body.m)*params.model.dyn.g*ones(1,length(tsim));
plot(tsim,F_list(2,:)+F_list(4,:),'b-','LineWidth',2);
hold on
plot(tsim,weight,'r-','LineWidth',1);
ylabel('Ground Reaction vs Weight (N)')
xlabel('time (sec)')
hold off

%% Un-comment if you want to animate the "jump"
% pause;
% 
% animate_robot(xsim(1:5,:),F_list,params,'trace_foot_com',true,...
%     'trace_body_com',true,'trace_spine_tip',true,'show_constraint_forces',true,'video',true);
% fprintf('Done!\n');


%% Event function for ODE45%
% Description:
%   Event function that is called when a constraint becomes inactive (or, in the future, active) 
%
% Inputs:
%   tsim: the array of time values selected by ode45
%   xsim: the 10xlength(tsim) array of states computed by ode45;
%   variables:  a struct with many elements
%   params:  a struct with many elements
%
% Outputs:
%   value
%   isterminal
%   direction
function [value,isterminal,direction] = robot_events_nested(t,x)

    % for convenience, define q and q_dot
    nq = numel(x)/2;    % assume that x = [q;q_dot];
    q = x(1:nq);
    q_dot = x(nq+1:2*nq);

    % solve for control inputs at this instant
    tau_s = interp1(variables.motor.spine.time,variables.motor.spine.torque,t);
    tau_m = interp1(variables.motor.body.time,variables.motor.body.torque,t);
    Q = [0;0;0;tau_s;tau_m];

    % find the parts that don't depend on constraint forces
    H = H_eom(x,params);
    M = mass_matrix(x,params);
    Minv = inv_mass_matrix(x,params);

    % build the constraints and forces 
    %event_variables = variables.foot.constraints
    nc = variables.foot.constraints(1) + 2*variables.foot.constraints(2);  % number of active contacts
    switch nc  
        case 0      % both feet are off the ground
            value = 1;
            isterminal = 0;
            direction = 0;
        case 1      % left foot is on the ground and right is off
            [A_all,Hessian] = constraint_derivatives(x,params);
            A = A_all([1,2],:);
            Adotqdot = [q_dot'*Hessian(:,:,1)*q_dot;
                        q_dot'*Hessian(:,:,2)*q_dot ];
            F = inv(A*Minv*A')*(A*Minv*(Q - H) + Adotqdot);
            value = F(2);
            isterminal = 1;
            direction - -1;
        case 2      % right foot is on the ground and left is off
            [A_all,Hessian] = constraint_derivatives(x,params);
            A = A_all([3,4],:);
            Adotqdot = [q_dot'*Hessian(:,:,3)*q_dot;
                        q_dot'*Hessian(:,:,4)*q_dot ];
            F = inv(A*Minv*A')*(A*Minv*(Q - H) + Adotqdot);
            value = F(2);
            isterminal = 1;
            direction = -1;
        case 3      % both feet are on the ground
            [A_all,Hessian] = constraint_derivatives(x,params);
            A = A_all([1,2,4],:);
            Adotqdot = [q_dot'*Hessian(:,:,1)*q_dot;
                        q_dot'*Hessian(:,:,2)*q_dot;
                        q_dot'*Hessian(:,:,4)*q_dot ];
            F = inv(A*Minv*A')*(A*Minv*(Q - H) + Adotqdot);
            value = F(2:3);
            isterminal = ones(2,1);
            direction = -ones(2,1);
    end
    

end
%% End of nested function


end