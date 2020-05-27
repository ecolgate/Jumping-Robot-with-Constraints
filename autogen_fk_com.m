function [FK_com_f,FK_com_s,FK_com_b] = autogen_fk_com(h_f,l_f,l_s,m_b,m_f,m_s,r,theta_f,theta_m,theta_s,w_f)
%AUTOGEN_FK_COM
%    [FK_COM_F,FK_COM_S,FK_COM_B] = AUTOGEN_FK_COM(H_F,L_F,L_S,M_B,M_F,M_S,R,THETA_F,THETA_M,THETA_S,W_F)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    27-May-2020 10:08:28

t2 = cos(theta_f);
t3 = sin(theta_f);
t4 = l_s.*m_s;
t5 = m_b+m_s;
t6 = theta_f+theta_s;
t11 = m_b.*r.*theta_m;
t7 = h_f.*t2;
t8 = l_f.*t2;
t9 = m_f+t5;
t10 = cos(t6);
t12 = sin(t6);
t13 = t3.*w_f;
t15 = t4+t11;
t14 = 1.0./t9;
t16 = t5.*t14;
t18 = t14.*t15;
FK_com_f = [t12.*t18+l_f.*t3.*t16;t7+t13];
if nargout > 1
    t17 = t16-1.0;
    t19 = -t18;
    t20 = l_f.*t3.*t17;
    FK_com_s = [t20-t12.*(l_s+t19);t7+t8+t13+l_s.*t10];
end
if nargout > 2
    FK_com_b = [t20+t12.*(t18-r.*theta_m);t7+t8+t13+r.*t10.*theta_m];
end
