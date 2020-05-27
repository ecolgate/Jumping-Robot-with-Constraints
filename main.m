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

%% Set max time step using odeset
options = odeset('MaxStep',params.control.dt);

%% Set the initial equilibrium pose of the robot
x_eq = zeros(6,1);
% Set foot angle to 30 degrees and spine to -5 degrees
x_eq(1) = pi/6;
x_eq(2) = -pi/36;
x_eq(3) = equilibrium_motor_angle(x_eq,params);

% Code for testing equilibrium 
%dist_up_spine = params.model.geom.body.r*x_eq(3)
%[FK_com_f,FK_com_s,FK_com_b] = fk_com(x_eq,params);
%should_be_zero = params.model.dyn.foot.m*FK_com_f(1) + ...
%    params.model.dyn.spine.m*FK_com_s(1) + ...
%    params.model.dyn.body.m*FK_com_b(1)
%x_wheel = wheel_position(x_eq,params)  % should be zero

%% Show the robot in equilibrium
plot_robot(x_eq(1:3),params,'new_fig',true);

%% Set the initial equilibrium motor torques
G = conservative_forces(x_eq,params);
u_eq = [G(2);G(3)];              % initial command
u = u_eq;

%% Design the LQR Controller
% first, find the linearized equations: (more info in
% derive_equations_JR.mlx)
%   x_dot = [ 0,  I; 0, -M(q_eq)\G_jac(q_eq)]*x + ...
%       M(q_eq)\[0,0;eye(2)]*u_lin
M_eq = mass_matrix(x_eq,params);    % mass matrix at equilibrium
G_jac_eq = derivative_conservative_forces(x_eq,params); % gravitational and spring forces at equilibrium
% create A and B matrices of linear system
A = [zeros(3,3),eye(3);-M_eq\G_jac_eq,zeros(3,3)];
B = [zeros(3,2);M_eq\[0,0;eye(2)]];
% then, set up Q and R matrices using Bryson's Rule
% I had to play with the weights quite a bit to get something reasonable.
% I ended up putting very low penalty on the angles, none on the angular
% velocities, and max penalty on the actuation (R_lqr).  This effectively
% turns down the gain, which is important for digital control.  Too much
% gain, and things get unstable
Q_lqr = diag([.0000001,.0000001,.0000001,0,0,0]);
R_lqr = eye(2);
% then solve for optimal gains
[Gains,~,Poles] = lqr(A,B,Q_lqr,R_lqr);
% Poles  % uncomment this line if you want to see the closed loop poles
% pause

%% initialize the state
x_IC = x_eq;
% perturb slightly so that something happens
x_IC(3) = x_IC(3) + 2*pi;

%% initialize controller memory
memory.u_eq = u_eq;  % not really memory ... just passing the equilibrium 
memory.x_eq = x_eq;  % not really memory ... just passing the equilibrium
memory.y = x_IC(2:4);
memory.est1_theta_f = x_IC(1);  % perfect initial knowledge!

%% set up control timing
t_write = 0.0;          % starting time
dt = params.control.dt; % controller time step

%% start with null matrices for holding results 
tsim = [];
xsim = [];
usim = [];
tcontrol = [];

%% the main loop
while t_write < params.sim.tfinal
    
    % simulate the plant from this write to the next write
    tspan = [t_write, t_write+dt];
    [tseg,xseg] = analog_plant(tspan,u,x_IC);
    
    % find the state and sensor measurements at the time a read was made
    t_read = tseg(end) - params.control.delay;
    x_read = interp1(tseg,xseg,t_read);
    y = sensor(x_read);
    
    % compute the control 
    [u,memory] = digital_controller(y,Gains,memory);
    
    % update t_write and x_IC for next iteration
    t_write = tseg(end);
    x_IC = xseg(end,:); 

    % store variables for plotting   
    % variables based off of tsim
    tsim = [tsim;tseg];
    xsim = [xsim;xseg];
    % variables based off of tcontrol, not tsim
    tcontrol = [tcontrol;t_read];
    usim = [usim;u];
        
end


%%  Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot joint angles
figure;
plot(tsim,xsim(:,1:3),'LineWidth',2);  
ylabel('Joint Angles')
xlabel('time (sec)')

% Animate Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's resample the simulator output so we can animate with evenly-spaced 
% points in (time,state). 
% 1) deal with possible duplicatetimes in tsim: 
% (https://www.mathworks.com/matlabcentral/answers/321603-how-do-i-interpolate-1d-data-if-i-do-not-have-unique-values
tsim = cumsum(ones(size(tsim)))*10*eps + tsim;

% 2) resample the duplicate-free time vector: 
t_anim = 0:params.viz.dt:tsim(end);

% 3) resample the state-vs-time array: 
x_anim = interp1(tsim,xsim,t_anim); x_anim = x_anim';  % transpose so that xsim is 12xN (N = number of timesteps)

animate_robot(x_anim([1:3],:),params,'trace_foot_com',true,...
    'trace_body_com',true,'trace_spine_tip',true,'video',true);

fprintf('Done!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end of main.m, except for nested functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% BELOW HERE ARE THE NESTED FUNCTIONS, SENSOR, DIGITAL_CONTROLLER, ANALOG_PLANT, AND ROBOT_DYNAMICS
%% THEY HAVE ACCESS TO ALL VARIABLES IN MAIN

%% sensor.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Computes the sensor values that are sent to the digital controller 
%
% Inputs:
%   t_read: time at which sensors are read
%   x_read: the 12x1 state vector at the time of a read (note that not all
%   of the state is actually read ... we must determine the sensor readings
%   from x_read)
%   u: the control inputs (used here because we have to call
%   robot_dynamics)
%
% Outputs:
%   y: the sensor values

function [y] = sensor(x_read)
    
    % NOTE:  right now, sensors are "perfect" -- no noise or quantization.
    % That *should* be added!
    y = zeros(3,1);
    % assume encoders for spine angle and body motor angle
    y(1:2) = x_read(2:3);   % theta_s and theta_m
    % assume gyro for foot angular velocity 
    y(3) = x_read(4);
    
end
%% end of sensor.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%% digital_controller.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Computes the torque commands to the two motors 
%
% Inputs:
%   y: the 5x1 output vector at the time of a read 
%   memory: a struct that holds stored variables
%
% Outputs:
%   u: the two torque commands
%   memory: a struct that holds stored variables

function [u,memory] = digital_controller(y,Gains,memory)
    
    % estimate theta_s_dot and theta_m_dot by backwards difference
    % differentiation
    % NOTE: some low pass filtering should probably be added
    est_theta_s_dot = (y(1) - memory.y(1))/params.control.dt;
    est_theta_m_dot = (y(2) - memory.y(2))/params.control.dt;
    
    % estimate theta_f from IMU readings
    % two approaches:  1) integrate theta_f_dot;  2) use the measured
    % accelerations and the knowledge of track shape to estimate it.  This
    % approach depends on a model, but doesn't have the problem of drift
    %
    % Approach 1:  Integrate gyro
    est1_theta_f = memory.est1_theta_f + y(3)*params.control.dt;
    
    % package up the state estimate
    x = [est1_theta_f;y(1);y(2);y(3);est_theta_s_dot;est_theta_m_dot];
    
    % compute the controls
    error = x - memory.x_eq;
    u = memory.u_eq - Gains*error;
    
    % I've removed saturation to make it work.  If I put this back, things
    % go haywire after a while.  Probably by increasing the gain a bit, we
    % could bring this back.  For now, it might be better just to plot
    % actuator torques
    
    % Saturate u (i.e., observe actuator limitations!)
%     if abs(u(1)) > params.motor.spine.peaktorque
%         u(1) = params.motor.spine.peaktorque*sign(u(1));
%     end
%     if abs(u(2)) > params.motor.body.peaktorque
%         u(2) = params.motor.body.peaktorque*sign(u(2));
%     end
    
    % Update memory (these are values that the Tiva would store)
    memory.y = y;
    memory.est1_theta_f = est1_theta_f;
    
end
%% end of digital_controller.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% analog_plant.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Updates the state of the robot (and, in future, motor torques) over a
%   time span typically consisting of one ZOH period.
%
%   Also computes the constraint forces
%
% Inputs:
%   tspan: [start time, end time]
%   u: torque commands
%   x_IC: the 6x1 state vector at the start time
%
% Outputs:
%   tseg: the vector of times at which states were computed
%   xseg: the states at those times

function [tseg,xseg] = analog_plant(tspan,u,x_IC)

    % NOTE: consider adding a model of motor dynamics:  right now, we
    % assume that applied torque = commanded torque
    
    % integrate the robot dyamics
    [tseg, xseg] = ode45(@(t,x) robot_dynamics(t,x,u), tspan, x_IC);

end
%% end of analog_plant.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% robot_dynamics.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%
%   Computes the derivative of the state:
%       x_dot(1:3) = x(4:6)
%       x_dot(4:6) = inv(M)*(Q - H)
%
% Inputs:
%   t: time (scalar)
%   x: the 6x1 state vector
%   u: the 2x1 motor torque commands
%
% Outputs:
%   dx: derivative of state x with respect to time.

function [dx] = robot_dynamics(~,x,u)

    dx = zeros(numel(x),1);
    nq = 3;    
    % for convenience, define q_dot 
    q_dot = x(nq+1:2*nq);

    % set up generalized forces based on motor torques
    Q = [0;u];

    % find the parts that don't depend on constraint forces
    H = H_eom(x,params);
    M = mass_matrix(x,params); % NOTE: eliminated use of symbolic inverse mass matrix since it takes so long to compute

    dx(1:nq) = q_dot;  
    dx(nq+1:2*nq) = M\(Q - H);

end
%% end of robot_dynamics.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
%% End of main.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%