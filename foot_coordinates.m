%% foot_coordinates.m
%
% Description:
%   Wrapper function for autogen_fk_feet.m
%   Computes the locations of the left and right feet of the jumping robot.
%
% Inputs:
%   x: the state vector, x = [q; q_dot; s];
%   params: a struct with many elements, generated by calling init_params.m
%
% Outputs:
%   r_l, r_r = foot coordinates

function [r_l,r_r] = foot_coordinates(x,params)

x_f = x(1);
z_f = x(2);
theta_f = x(3);

h_f = params.model.geom.foot.hbot;
w_f = params.model.geom.foot.w;

[r_l,r_r] = autogen_fk_feet(h_f,theta_f,w_f,x_f,z_f);

end
