function H = autogen_H_eom(K,g,l_f,l_s,m_b,m_f,m_s,p_b,r,theta_f,theta_m,theta_s,theta_dot_f,theta_dot_m,theta_dot_s)
%AUTOGEN_H_EOM
%    H = AUTOGEN_H_EOM(K,G,L_F,L_S,M_B,M_F,M_S,P_B,R,THETA_F,THETA_M,THETA_S,THETA_DOT_F,THETA_DOT_M,THETA_DOT_S)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    20-Apr-2020 15:47:35

t2 = cos(theta_f);
t3 = cos(theta_s);
t4 = sin(theta_f);
t5 = sin(theta_s);
t6 = theta_f+theta_s;
t7 = r.^2;
t8 = theta_dot_f.^2;
t9 = theta_dot_s.^2;
t10 = cos(t6);
t11 = sin(t6);
t15 = m_b.*t7.*theta_m.*theta_dot_f.*theta_dot_m.*2.0;
t16 = m_b.*t7.*theta_m.*theta_dot_m.*theta_dot_s.*2.0;
t12 = g.*m_b.*p_b.*t10;
t13 = g.*l_s.*m_s.*t11;
t14 = g.*m_b.*r.*t11.*theta_m;
t17 = -t12;
t18 = -t13;
t19 = -t14;
H = [l_f.*m_b.*t4.*t8+l_f.*m_s.*t4.*t8+l_s.*m_s.*t8.*t11+l_s.*m_s.*t9.*t11+m_b.*p_b.*t8.*t10+m_b.*p_b.*t9.*t10+l_s.*m_s.*t11.*theta_dot_f.*theta_dot_s.*2.0+m_b.*p_b.*t10.*theta_dot_f.*theta_dot_s.*2.0+m_b.*r.*t8.*t11.*theta_m+m_b.*r.*t9.*t11.*theta_m-m_b.*r.*t10.*theta_dot_f.*theta_dot_m.*2.0-m_b.*r.*t10.*theta_dot_m.*theta_dot_s.*2.0+m_b.*r.*t11.*theta_m.*theta_dot_f.*theta_dot_s.*2.0;g.*m_b+g.*m_f+g.*m_s-l_f.*m_b.*t2.*t8-l_f.*m_s.*t2.*t8-l_s.*m_s.*t8.*t10-l_s.*m_s.*t9.*t10+m_b.*p_b.*t8.*t11+m_b.*p_b.*t9.*t11-l_s.*m_s.*t10.*theta_dot_f.*theta_dot_s.*2.0+m_b.*p_b.*t11.*theta_dot_f.*theta_dot_s.*2.0-m_b.*r.*t8.*t10.*theta_m-m_b.*r.*t9.*t10.*theta_m-m_b.*r.*t11.*theta_dot_f.*theta_dot_m.*2.0-m_b.*r.*t11.*theta_dot_m.*theta_dot_s.*2.0-m_b.*r.*t10.*theta_m.*theta_dot_f.*theta_dot_s.*2.0;t15+t16+t17+t18+t19-g.*l_f.*m_b.*t4-g.*l_f.*m_s.*t4-l_f.*l_s.*m_s.*t5.*t9-l_f.*m_b.*p_b.*t3.*t9-l_f.*l_s.*m_s.*t5.*theta_dot_f.*theta_dot_s.*2.0-l_f.*m_b.*p_b.*t3.*theta_dot_f.*theta_dot_s.*2.0-l_f.*m_b.*r.*t5.*t9.*theta_m+l_f.*m_b.*r.*t3.*theta_dot_f.*theta_dot_m.*2.0+l_f.*m_b.*r.*t3.*theta_dot_m.*theta_dot_s.*2.0-l_f.*m_b.*r.*t5.*theta_m.*theta_dot_f.*theta_dot_s.*2.0;t15+t16+t17+t18+t19+K.*theta_s+l_f.*l_s.*m_s.*t5.*t8+l_f.*m_b.*p_b.*t3.*t8+l_f.*m_b.*r.*t5.*t8.*theta_m;-m_b.*r.*(-g.*t10+l_f.*t3.*t8+r.*t8.*theta_m+r.*t9.*theta_m+r.*theta_m.*theta_dot_f.*theta_dot_s.*2.0)];
