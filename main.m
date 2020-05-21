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

%% Simulate the robot forward in time 
x_IC = params.x_IC';    % initial conditions
u = [0;0];              % initial command

%% initialize controller memory
memory.u = u;
memory.y = zeros(5,1);
memory.est1_theta_f = 0;
memory.filt1_theta_f = 0;
memory.est2_theta_f = 0;
memory.filt2_theta_f = 0;

%% set up control timing
t_write = 0.0;          % starting time
dt = params.control.dt; % controller time step

%% start with null matrices for holding results 
tsim = [];
xsim = [];
Psim = [];
F_list = [];
usim = [];
tcontrol = [];
theta_f_estc = [];

%% the main loop
while t_write < params.sim.tfinal
    
    % simulate the plant from this write to the next write
    tspan = [t_write, t_write+dt];
    [tseg,xseg,Fseg,Pseg] = analog_plant(tspan,u,x_IC);
    
    % find the state at the time a read was made
    t_read = tseg(end) - params.control.delay;
    x_read = interp1(tseg,xseg,t_read);
    
    % compute the control 
    [u,memory] = digital_controller(t_read,x_read,memory);
    
    % update t_write and x_IC for next iteration
    t_write = tseg(end);
    x_IC = xseg(end,:); 

    % store variables for plotting   
    % variables based off of tsim
    tsim = [tsim;tseg];
    xsim = [xsim;xseg];
    F_list = [F_list;Fseg'];
    Psim = [Psim;Pseg];
    % variables based off of tcontrol, not tsim
    tcontrol = [tcontrol;t_read];
    usim = [usim;u];
    theta_f_estc = [theta_f_estc;memory.filt1_theta_f+memory.filt2_theta_f];  % CHECK THIS
        
end

%% compute total energy as a function of time
ntsim = length(tsim);
Energy = zeros(1,ntsim-1);
KineticEnergy = zeros(1,ntsim-1);
PotentialEnergy = zeros(1,ntsim-1);
NetEnergy = zeros(1,ntsim-1);
actuator_work = 0;
for i=1:ntsim-1
    actuator_work = actuator_work + Psim(i+1)*(tsim(i+1)-tsim(i));
    Energy(i) = total_energy(xsim(i,:)',params);
    [KineticEnergy(i),PotentialEnergy(i)] = energy(xsim(i,:)',params);
    NetEnergy(i) = Energy(i) - actuator_work;
end


%%  Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot energy
figure;
plot(tsim(1:ntsim-1),Energy,'k-','LineWidth',2);
hold on;
plot(tsim(1:ntsim-1),NetEnergy,'r-','LineWidth',2);
plot(tsim(1:ntsim-1),KineticEnergy,'g-','LineWidth',2);
plot(tsim(1:ntsim-1),PotentialEnergy,'b-','LineWidth',2);
hold off
ylabel('Energy (J)')
xlabel('time (sec)')
legend('Total Energy','Net Energy','Kinetic Energy','Potential Energy')

% plot foot angle versus estimate
figure;
plot(tsim,xsim(:,3),'k-','LineWidth',2);  
hold on
plot(tcontrol,theta_f_estc,'r-','LineWidth',1);
ylabel('Foot Angle versus Estimate')
xlabel('time (sec)')
legend('Foot Angle','Estimate')
hold off

%% Animate Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's resample the simulator output so we can animate with evenly-spaced 
% points in (time,state). 
% 1) deal with possible duplicatetimes in tsim: 
% (https://www.mathworks.com/matlabcentral/answers/321603-how-do-i-interpolate-1d-data-if-i-do-not-have-unique-values
tsim = cumsum(ones(size(tsim)))*10*eps + tsim;

% 2) resample the duplicate-free time vector: 
t_anim = 0:params.viz.dt:tsim(end);

% 3) resample the state-vs-time array: 
x_anim = interp1(tsim,xsim,t_anim); x_anim = x_anim';  % transpose so that xsim is 12xN (N = number of timesteps)

% 4) resample the constraint forces-vs-time array: 
F_anim = interp1(tsim,F_list,t_anim); F_anim = F_anim';

animate_robot(x_anim([1:5,11,12],:),F_anim,params,'trace_foot_com',true,...
    'trace_body_com',true,'trace_spine_tip',true,'show_constraint_forces',true,'video',true);

fprintf('Done!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% end of main.m, except for nested functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% BELOW HERE ARE THE NESTED FUNCTIONS, DIGITAL_CONTROLLER, ANALOG_PLANT, AND ROBOT_DYNAMICS
%% THEY HAVE ACCESS TO ALL VARIABLES IN MAIN

%% digital_controller.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Computes the torque commands to the two motors 
%
% Inputs:
%   t_read: time at which sensors are read
%   x_read: the 12x1 state vector at the time of a read (note that not all
%   of the state is actually read ... we must determine the sensor readings
%   from x_read)
%   memory: a struct that holds stored variables
%
% Outputs:
%   u: the two torque commands
%   memory: a struct that holds stored variables

function [u,memory] = digital_controller(t_read,x_read,memory)
    
    % figure out what sensors can see
    % NOTE:  right now, sensors are "perfect" -- no noise or quantization.
    % That *should* be added!
    y = zeros(5,1);
    % assume encoders for spine angle and body motor angle
    y(1:2) = x_read(4:5);   % theta_s and theta_m
    % assume gyro for foot angular velocity 
    y(3) = x_read(8);
    % assume that two axes of acceleration are measured in the frame of the
    % foot.  y(4) is the x-axis acceleration when the foot is flat
    % y(5) is the z-axis acceleration (including gravity) when the foot is
    % flat
    theta_f = x_read(3);
    rot = [cos(theta_f), sin(theta_f); -sin(theta_f), cos(theta_f)];
    [dx,~] = robot_dynamics(t_read,x_read',memory.u);
    accel_actual = [dx(6);dx(7)-params.model.dyn.g];
    y(4:5) = rot*accel_actual;       % measured accelerations
    
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
    %
    % Approach 2:  Use gyro + acceleration + model of track
    theta_f_ddot = (y(3) - memory.y(3))/params.control.dt;  % foot angular acceleration
    rtrack = params.model.geom.track.r;  % radius of track
    est2_theta_f = atan2(rtrack*theta_f_ddot - y(4),rtrack*y(3)^2 - y(5));
    %
    % Combined Estimate: high pass est1, low pass est2, and add
    filt1_theta_f = params.control.foot.alpha*(memory.filt1_theta_f + est1_theta_f - memory.est1_theta_f);
    filt2_theta_f = params.control.foot.alpha*memory.filt2_theta_f + (1-params.control.foot.alpha)*est2_theta_f;
    est_combined_theta_f = filt1_theta_f + filt2_theta_f;
    
    % PUMP!!!  Based on theta_f and theta_f_dot, move body up or down
    if est_combined_theta_f>=0 && y(3)>=0
        theta_m_des = params.control.theta_m_high;
    elseif est_combined_theta_f>=0 && y(3)<0
        theta_m_des = params.control.theta_m_low;
    elseif est_combined_theta_f<0 && y(3)<0
        theta_m_des = params.control.theta_m_high;
    else
        theta_m_des = params.control.theta_m_low;  
    end
        
    % Compute u with PD control
    % NOTE: the desired velocity here is always zero.  Basically, we are
    % moving the robot from one fixed pose to another.  LQR is also a good
    % approach to this.
    u(2,1) = params.control.body.kp*(theta_m_des - y(2)) + ...
        params.control.body.kd*(0 - est_theta_m_dot);
    u(1,1) = params.control.spine.kp*(params.control.theta_s_des - y(1)) + ...
        params.control.spine.kd*(0 - est_theta_s_dot);  % - u(2,1);
    
    % Saturate u (i.e., observe actuator limitations!)
    if abs(u(1)) > params.motor.spine.peaktorque
        u(1) = params.motor.spine.peaktorque*sign(u(1));
    end
    if abs(u(2)) > params.motor.body.peaktorque
        u(2) = params.motor.body.peaktorque*sign(u(2));
    end
    
    % Update memory (these are values that the Tiva would store)
    memory.u = u;
    memory.y = y;
    memory.est1_theta_f = est1_theta_f;
    memory.est2_theta_f = est2_theta_f;
    memory.filt1_theta_f = filt1_theta_f;
    memory.filt2_theta_f = filt2_theta_f;
    
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
%   x_IC: the 12x1 state vector at the start time
%
% Outputs:
%   tseg: the vector of times at which states were computed
%   xseg: the states at those times
%   Fseg: the constraint forces at those times

function [tseg,xseg,Fseg,Pseg] = analog_plant(tspan,u,x_IC)

    % NOTE: consider adding a model of motor dynamics:  right now, we
    % assume that applied torque = commanded torque
    
    % integrate the robot dyamics
    [tseg, xseg] = ode45(@(t,x) robot_dynamics(t,x,u), tspan, x_IC);

    % compute the actuator power that was applied/drawn during segment
    Pseg = u(1)*xseg(:,9) + u(2)*xseg(:,10);

    % compute the constraint forces that were active during the segment
    Fseg = zeros(4,length(tseg));
    for i=1:length(tseg)
        [~,Fseg(:,i)] = robot_dynamics(tseg(i),xseg(i,:)',u);
    end

end
%% end of analog_plant.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% robot_dynamics.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description:
%   Computes the constraint forces: 
%       Fnow = inv(A*Minv*A')*(A*Minv*(Q-H) + Adotqdot)
%
%   Also computes the derivative of the state:
%       x_dot(1:5) = (I - A'*inv(A*A')*A)*x(6:10)
%       x_dot(6:10) = inv(M)*(Q - H - A'F)
%       x_dot(11:12) = [ds_left; ds_right];
%
% Inputs:
%   t: time (scalar)
%   x: the 10x1 state vector
%   u: the 2x1 motor torque commands
%
% Outputs:
%   dx: derivative of state x with respect to time.
%   F:  4x1 vector of constraint forces

function [dx,F] = robot_dynamics(~,x,u)

    dx = zeros(numel(x),1);
    nq = numel(x)/2-1;    % assume that x = [q;q_dot;sl;sr];
    % for convenience, define q_dot and s
    q_dot = x(nq+1:2*nq);
    s = x(2*nq+1:2*nq+2);

    % set up generalized forces based on motor torques
    Q = [0;0;0;u];

    % find the parts that don't depend on constraint forces
    H = H_eom(x,params);
    M = mass_matrix(x,params); % NOTE: eliminated use of symbolic inverse mass matrix since it takes so long to compute
    Minv = inv(M); 
    [A,Hessian] = constraint_derivatives(x,params);

    % solve for constraint forces, if any; compute velocities and accelerations
    n_active_constraints = sum(params.sim.constraints.active);
    A = A(params.sim.constraints.active,:);   % eliminate inactive rows
    Hessian_active = Hessian(:,:,params.sim.constraints.active);  % eliminate inactive hessians
    Adotqdot = zeros(n_active_constraints,1);
    for ic=1:n_active_constraints
        Adotqdot(ic) = q_dot'*Hessian_active(:,:,ic)*q_dot;
    end

    % compute the constraint forces and accelerations
    F_active = (A*Minv*A')\(A*Minv*(Q - H) + Adotqdot); % these are the constraint forces
    dx(1:nq) = (eye(nq) - A'*((A*A')\A))*q_dot;  
    % NOTE:  important change to next line:  subtract off an error term to
    % ensure that velocities converge to values that observe constraints
    dx(nq+1:2*nq) = Minv*(Q - H - A'*F_active) - A'*((A*A')\A)*q_dot/params.control.dt;

    % 4x1 vector of constraint forces
    F = zeros(4,1);
    F(params.sim.constraints.active) = F_active;   

    % compute the pathlength derivatives for constraint tracking
    [r_l,r_r] = foot_coordinates(x,params);
    [v_l,v_r] = foot_velocities(x,params);
    [p_l,t_l,~] = track(s(1),params);
    [p_r,t_r,~] = track(s(2),params);
    dx(2*nq+1) = (v_l + params.sim.constraints.gain*(r_l - p_l))'*t_l;
    dx(2*nq+2) = (v_r + params.sim.constraints.gain*(r_r - p_r))'*t_r;

end
%% end of robot_dynamics.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
%% End of main.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%