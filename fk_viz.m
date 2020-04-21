%% fk_viz.m
%
% Description:
%   Wrapper function for autogen_fk_viz.m
%   Computes the locations of the pivot and the spine tip of the jumping robot.
%
% Inputs:
%   x: the state vector, x = [q; q_dot];
%   params: a struct with many elements, generated by calling init_params.m
%
% Outputs:
%   FK_pivot = x and z locations of pivot
%   FK_tip = x and z locations of spine tip

function [FK_pivot,FK_tip] = fk_viz(x,params)

x_f = x(1);
z_f = x(2);
theta_f = x(3);
theta_s = x(4);

l_f = params.model.geom.foot.htop;
l_spine = params.model.geom.spine.h;

[FK_pivot,FK_tip] = autogen_fk_viz(l_f,l_spine,theta_f,theta_s,x_f,z_f);

end
