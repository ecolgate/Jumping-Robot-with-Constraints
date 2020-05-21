%% constraints.m
%
% Description:
%   Wrapper function for autogen_constraints.m
%   Computes the four constraint equations (left foot, right foot, top of
%   spine, bottom of spine)
%
% Inputs:
%   x: the state vector, x = [q; q_dot];
%   params: a struct with many elements, generated by calling init_params.m
%
% Outputs:
%   C_unilateral = values of four constraint equations

function [C_unilateral] = constraints(x,params)

x_f = x(1);
z_f = x(2);
theta_f = x(3);
theta_m = x(5);

h_b = params.model.geom.body.h;
h_f = params.model.geom.foot.hbot;
l_spine = params.model.geom.spine.h;
r = params.model.geom.body.r;
w_f = params.model.geom.foot.w;

[p_l,~,n_l] = track(x(11),params);
[p_r,~,n_r] = track(x(12),params);

p_lx = p_l(1);
p_lz = p_l(2);
p_rx = p_r(1);
p_rz = p_r(2);
n_lx = n_l(1);
n_lz = n_l(2);
n_rx = n_r(1);
n_rz = n_r(2);

C_unilateral = autogen_constraints(h_b,h_f,l_spine,n_lx,n_lz,n_rx,n_rz,p_lx,p_lz,p_rx,p_rz,r,theta_f,theta_m,w_f,x_f,z_f);
end

