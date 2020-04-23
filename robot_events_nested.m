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
function [value,isterminal,direction] = robot_events_nested(t,...
                                                            x,...
                                                            variables,...
                                                            params)

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