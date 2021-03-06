function M = autogen_mass_matrix(J_b,J_f,J_m,J_s,l_f,l_s,m_b,m_f,m_s,p_b,r,theta_f,theta_m,theta_s)
%AUTOGEN_MASS_MATRIX
%    M = AUTOGEN_MASS_MATRIX(J_B,J_F,J_M,J_S,L_F,L_S,M_B,M_F,M_S,P_B,R,THETA_F,THETA_M,THETA_S)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    20-Apr-2020 15:43:35

t2 = cos(theta_f);
t3 = cos(theta_s);
t4 = sin(theta_f);
t5 = sin(theta_s);
t6 = theta_f+theta_s;
t7 = l_f.^2;
t8 = l_s.^2;
t9 = p_b.^2;
t10 = r.^2;
t11 = theta_m.^2;
t12 = m_b+m_f+m_s;
t14 = m_b.*p_b.*r;
t13 = cos(t6);
t15 = sin(t6);
t16 = m_s.*t8;
t17 = m_b.*t9;
t18 = l_f.*t2.*2.0;
t19 = l_f.*t4.*2.0;
t20 = l_f.*l_s.*m_s.*t3;
t21 = l_f.*m_b.*p_b.*t5;
t22 = l_f.*m_b.*r.*t5;
t23 = -t14;
t31 = l_f.*m_b.*r.*t3.*theta_m;
t35 = m_b.*t10.*t11;
t24 = l_s.*t13.*2.0;
t25 = p_b.*t13.*2.0;
t26 = l_s.*t15.*2.0;
t27 = p_b.*t15.*2.0;
t28 = l_s.*m_s.*t13;
t29 = m_b.*r.*t13;
t30 = l_s.*m_s.*t15;
t32 = m_b.*r.*t15;
t33 = r.*t15.*theta_m.*2.0;
t34 = J_m+t23;
t37 = -t21;
t38 = r.*t13.*theta_m.*2.0;
t36 = -t27;
t39 = -t28;
t40 = -t38;
t41 = -t30;
t42 = -t32;
t43 = t18+t24;
t44 = t19+t26;
t45 = t22+t34;
t46 = t25+t33;
t63 = J_b+J_m+J_s+t16+t17+t20+t31+t35+t37;
t47 = t27+t40;
t48 = (m_s.*t43)./2.0;
t49 = (m_s.*t44)./2.0;
t52 = t19+t46;
t53 = (m_b.*t46)./2.0;
t54 = t18+t36+t38;
t50 = -t48;
t51 = -t49;
t55 = -t53;
t56 = (m_b.*t47)./2.0;
t57 = (m_b.*t52)./2.0;
t58 = (m_b.*t54)./2.0;
t59 = -t57;
t60 = -t58;
t61 = t39+t56;
t62 = t41+t55;
t64 = t51+t59;
t65 = t50+t60;
M = reshape([t12,0.0,t65,t61,t42,0.0,t12,t64,t62,t29,t65,t64,J_b+J_f+J_m+J_s+t16+t17+t20.*2.0-t21.*2.0+t31.*2.0+t35+m_b.*t7+m_s.*t7,t63,t45,t61,t62,t63,J_b+J_m+J_s+t16+t17+t35,t34,t42,t29,t45,t34,J_m+m_b.*t10],[5,5]);
