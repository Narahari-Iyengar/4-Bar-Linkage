% 4-bar linkage DAE
clear all;
clc;

syms m1 m2 m3 m4
syms d1 d2 d3 d4
syms l1 l2 l3 l4
syms I1 I2 I3 I4
syms x1 x2 x3 x4 y1 y2 y3 y4 th1 th2 th3 th4
syms x1d x2d x3d x4d y1d y2d y3d y4d th1d th2d th3d th4d
syms x1dd y1dd th1dd x2dd y2dd th2dd x3dd y3dd th3dd x4dd y4dd th4dd
syms Rhox Rhoy Rh1y R12x R12y R23x R23y R34x R34y R4ox R4oy
syms g
%Define the axis
i=[1 0 0];
j=[0 1 0];
k=cross(i,j);

er1=i*sin(th1)-j*cos(th1);
er2=i*sin(th2)-j*cos(th2);
er3=i*sin(th3)-j*cos(th3);
er4=i*sin(th4)-j*cos(th4);

eth1=cross(k,er1);
eth2=cross(k,er2);
eth3=cross(k,er3);
eth4=cross(k,er4);

%Define lengths
rg1o=d1*er1;
rao=l1*er1;
rg2a=d2*er2;
rba=l2*er2;
rg3b=d3*er3;
rcb=l3*er3;
rg4c=d4*er4;
rdc=l4*er4;

%Define lengths wrt origin
rg1o=d1*er1;
rao=l1*er1;
rg2o=rg2a+rao;
rbo=rba+rao;
rg3o=rg3b+rbo;
rco=rcb+rbo;
rg4o=rg4c+rco;
rdo=rdc+rco;

rag1=rao-rg1o;
rbg2=rbo-rg2o;
rcg3=rco-rg3o;
rdg4=rdo-rg4o;

% %Define the accelerations
% ag1o=-d1*th1d.^2*er1+d1*th1dd*eth1;
% % ag1o=0;
% aao=-l1*th1d.^2*er1+l1*th1dd*eth1;
% % aao=0;
% ag2o=-d2*th2d.^2*er2+d2*th2dd*eth2+aao;
% abo=-l2*th2d.^2*er2+l2*th2dd*eth2+aao;
% ag3o=-d3*th3d.^2*er3+d3*th3dd*eth3+abo;
% aco=-l3*th3d.^2*er3+l3*th3dd*eth3+abo;
% ag4o=-d4*th4d.^2*er4+d4*th4dd*eth4+aco;
% ado=-l4*th4d.^2*er4+l4*th4dd*eth4+aco;

%Linear Moment Balance
L1L=(Rhox+R4ox-R12x)*i...
    +(-Rhoy-R4oy-Rh1y+R12y-m1*g)*j;
L1R=m1*(x1dd*i+y1dd*j);
% L1R=0;

L2L=(R12x-R23x)*i+(-R12y+R23y-m2*g)*j;
L2R=m2*(x2dd*i+y2dd*j);

L3L=(R23x-R34x)*i+(-R23y+R34y-m3*g)*j;
L3R=m3*(x3dd*i+y3dd*j);

L4L=(R34x-R4ox)*i+(-R34y+R4oy-m4*g)*j;
L4R=m4*(x4dd*i+y4dd*j);

%Angular Moment Balance about Gi
Mg1=cross(-rg1o,((Rhox+R4ox)*i+(-Rhoy-R4oy)*j))...
    +cross(rag1,((-R12x)*i+(-Rh1y+R12y)*j));
Hdg1=I1*th1dd*k;
% Hdg1=0;

Mg2=cross(-rg2a,(R12x*i-R12y*j))+cross(rbg2,(-R23x*i+R23y*j));
Hdg2=I2*th2dd*k;

Mg3=cross(-rg3b,(R23x*i-R23y*j))+cross(rcg3,(-R34x*i+R34y*j));
Hdg3=I3*th3dd*k;

Mg4=cross(-rg4c,(R34x*i-R34y*j))+cross(rdg4,(-R4ox*i+R4oy*j));
Hdg4=I4*th4dd*k;

%Constraint equations
C1L=x1dd*i+y1dd*j;
% C1L=0;
C1R=cross(th1dd*k,rg1o)-th1d.^2*rg1o;
% C1R=0;

C2L=x2dd*i+y2dd*j-cross(th2dd*k,rg2a)+th2d.^2*rg2a;
C2R=x1dd*i+y1dd*j+cross(th1dd*k,rag1)-th1d.^2*rag1;

C3L=x3dd*i+y3dd*j-cross(th3dd*k,rg3b)+th3d.^2*rg3b;
C3R=x2dd*i+y2dd*j+cross(th2dd*k,rbg2)-th2d.^2*rbg2;

C4L=x4dd*i+y4dd*j-cross(th4dd*k,rg4c)+th4d.^2*rg4c;
C4R=x3dd*i+y3dd*j+cross(th3dd*k,rcg3)-th3d.^2*rcg3;

C5L=0;
C5R=x4dd*i+y4dd*j+cross(th4dd*k,rdg4)-th4d.^2*rdg4;

C6L=th1dd*k;
C6R=0;

%Putting in equations
eqn1=L1L-L1R;
eqn2=L2L-L2R;
eqn3=L3L-L3R;
eqn4=L4L-L4R;
eqn5=Mg1-Hdg1;
eqn6=Mg2-Hdg2;
eqn7=Mg3-Hdg3;
eqn8=Mg4-Hdg4;
eqn9=C1L-C1R;
eqn10=C2L-C2R;
eqn11=C3L-C3R;
eqn12=C4L-C4R;
eqn13=C5L-C5R;
eqn14=C6L-C6R;

%Take compponents of equations
eqn1i=eqn1(1);
eqn1j=eqn1(2);
eqn2i=eqn2(1);
eqn2j=eqn2(2);
eqn3i=eqn3(1);
eqn3j=eqn3(2);
eqn4i=eqn4(1);
eqn4j=eqn4(2);
eqn5k=eqn5(3);
eqn6k=eqn6(3);
eqn7k=eqn7(3);
eqn8k=eqn8(3);
eqn9i=eqn9(1);
eqn9j=eqn9(2);
eqn10i=eqn10(1);
eqn10j=eqn10(2);
eqn11i=eqn11(1);
eqn11j=eqn11(2);
eqn12i=eqn12(1);
eqn12j=eqn12(2);
eqn13i=eqn13(1);
eqn13j=eqn13(2);
eqn14k=eqn14(3);

%Define equations array and unknowns array
eqns=[eqn1i; eqn1j;eqn2i; eqn2j;eqn3i; eqn3j;eqn4i; eqn4j;
    eqn5k;eqn6k;eqn7k;eqn8k;
    eqn9i; eqn9j;eqn10i; eqn10j; eqn11i; eqn11j;eqn12i; eqn12j;
    eqn13i; eqn13j;eqn14k];
unknowns=[x1dd, y1dd, th1dd, x2dd, y2dd, th2dd, x3dd, y3dd, th3dd, x4dd, y4dd, th4dd, Rhox, Rhoy, Rh1y, R12x, R12y, R23x, R23y, R34x, R34y, R4ox, R4oy];

%Create MATLAB function files
[M,b]=equationsToMatrix(eqns,unknowns);
matlabFunction(M,'file','mymassmatrixDAE');
matlabFunction(b,'file','myrhsstuffDAE');