%% fk_com.m
%
% Description:
%   Wrapper function for autogen_fk_com.m
%   Computes the locations of the CoMs of the spine and body.
%
% Inputs:
%   x: the state vector, x = [q; q_dot];
%   params: a struct with many elements, generated by calling init_params.m
%
% Outputs:
%   FK_com_s = x and z locations of spine CoM
%   FK_com_b = x and z locations of body CoM

function [FK_com_s,FK_com_b] = fk_com(x,params)

x_f = x(1);
z_f = x(2);
theta_f = x(3);
theta_m = x(5);
theta_s = x(4);

l_f = params.model.geom.foot.htop;
l_s = params.model.geom.spine.l;
p_b = params.model.geom.body.pb;
r = params.model.geom.body.r;

[FK_com_s,FK_com_b] = autogen_fk_com(l_f,l_s,p_b,r,theta_f,theta_m,theta_s,x_f,z_f);

end

