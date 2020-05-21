%% energy.m
%
% Description:
%   Wrapper function for autogen_energy.m
%   Computes the kinetic and potential energy of the jumping robot.
%
% Inputs:
%   x: the state vector, x = [q; q_dot];
%   params: a struct with many elements, generated by calling init_params.m
%
% Outputs:
%   KE: kinetic energy
%   PE: potential energy

function [KE, PE] = energy(x,params)

z_f = x(2);
theta_f = x(3);
theta_m = x(5);
theta_s = x(4);
x_dot_f = x(6);
z_dot_f = x(7);
theta_dot_f = x(8);
theta_dot_m = x(10);
theta_dot_s = x(9);

J_b = params.model.dyn.body.J;
J_f = params.model.dyn.foot.J;
J_m = params.model.dyn.body.Jm;
K = params.model.dyn.K;
g = params.model.dyn.g;
h_f = params.model.geom.foot.hbot;
J_s = params.model.dyn.spine.J;
l_f = params.model.geom.foot.htop;
l_s = params.model.geom.spine.l;
m_b = params.model.dyn.body.m;
m_f = params.model.dyn.foot.m;
m_s = params.model.dyn.spine.m;
p_b = params.model.geom.body.pb;
r = params.model.geom.body.r;

[KE,PE] = autogen_energy(J_b,J_f,J_m,J_s,K,g,h_f,l_f,l_s,m_b,m_f,m_s,p_b,r,theta_f,theta_m,theta_s,theta_dot_f,theta_dot_m,theta_dot_s,x_dot_f,z_dot_f,z_f);

end

