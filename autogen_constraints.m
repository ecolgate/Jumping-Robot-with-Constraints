function C_all = autogen_constraints(h_f,theta_f,w_f,x_f,z_f)
%AUTOGEN_CONSTRAINTS
%    C_ALL = AUTOGEN_CONSTRAINTS(H_F,THETA_F,W_F,X_F,Z_F)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    20-Apr-2020 15:40:42

t2 = cos(theta_f);
t3 = sin(theta_f);
t7 = -z_f;
t4 = h_f.*t2;
t5 = h_f.*t3;
t6 = t3.*w_f;
t8 = t2-1.0;
t9 = t8.*w_f;
C_all = [t5-t9+x_f;t4+t6+t7;t5+t9+x_f;t4-t6+t7];
