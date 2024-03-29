function M = mymassmatrixDAE(I1,I2,I3,I4,d1,d2,d3,d4,l1,l2,l3,l4,m1,m2,m3,m4,th1,th2,th3,th4)
%MYMASSMATRIXDAE
%    M = MYMASSMATRIXDAE(I1,I2,I3,I4,D1,D2,D3,D4,L1,L2,L3,L4,M1,M2,M3,M4,TH1,TH2,TH3,TH4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    21-Nov-2018 19:30:44

t2 = sin(th1);
t3 = d1.*t2;
t4 = cos(th1);
t5 = d1.*t4;
t6 = cos(th2);
t7 = sin(th2);
t8 = d2.*t7;
t9 = cos(th3);
t10 = sin(th3);
t11 = d3.*t10;
t12 = cos(th4);
t13 = sin(th4);
t14 = d4.*t13;
t15 = t5-l1.*t4;
t16 = d2.*t6;
t17 = l1.*t2;
t18 = t16-l2.*t6;
t19 = d3.*t9;
t20 = l2.*t7;
t21 = t19-l3.*t9;
t22 = d4.*t12;
t23 = l3.*t10;
t24 = t22-l4.*t12;
t25 = l4.*t13;
M = reshape([-m1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-I1,0.0,0.0,0.0,-t5,-t3,t15,t3-t17,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,-m2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-I2,0.0,0.0,0.0,0.0,-t16,-t8,t18,t8-t20,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-I3,0.0,0.0,0.0,0.0,0.0,-t19,-t11,t21,t11-t23,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-I4,0.0,0.0,0.0,0.0,0.0,0.0,-t22,-t14,t24,t14-t25,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-d1.*t4,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,t3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,t3-l1.*t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,t15,-d2.*t6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,0.0,-t3+t17,t8,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,0.0,t18,-d3.*t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,0.0,-t8+t20,t11,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,1.0,0.0,0.0,0.0,t21,-d4.*t12,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,-1.0,0.0,0.0,-t11+t23,t14,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,-t5,0.0,0.0,t24,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,1.0,t3,0.0,0.0,-t14+t25,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[23,23]);
