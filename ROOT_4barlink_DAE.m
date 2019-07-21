%4-bar linkage- DAE method
%2 plots -  1. Simulation  2. Energy Conservation Check

tic
clear all;
clc;

%Input parameters
p.d1=0.5; p.d2 =0.5; p.d3 =0.5; p.d4=0.5; %Distance of center of mass from hinge
p.l1=sqrt(2); p.l2 =2; p.l3 =1; p.l4=1; %Length
p.m1=1; p.m2=1; p.m3=1; p.m4=1; %Mass
p.I1=0.2;p.I2=0.2;p.I3=0.2;p.I4=0.2; %Moment of Inertia
p.g=1; %acceleration due to gravity

rate = 10;
dur   = 40;
ntimes = dur*rate;
tspan = linspace(0,dur,ntimes);

%Pack the variables
d1=p.d1; d2=p.d2; d3=p.d3; d4=p.d4;
l1=p.l1; l2=p.l2; l3=p.l3; l4=p.l4;
m1=p.m1; m2=p.m2; m3=p.m3; m4=p.m4;
I1=p.I1; I2=p.I2; I3=p.I3; I4=p.I4;
g=p.g;

%Define intial angles (orientation)
th10=pi/4; th20=pi; th30=-pi/2; th40=pi/100;
th1d0=0; th2d0=0; th3d0=0; th4d0=0;

%Derive cartesian co-ordinates of CoG in fixed frame
x10=d1*sin(th10);
x20=l1*sin(th10)+ d2*sin(th20);
x30=l1*sin(th10)+l2*sin(th20)+d3*sin(th30);
x40=l1*sin(th10)+l2*sin(th20)+l3*sin(th30)+d4*sin(th40);
y10=-d1*cos(th10);
y20=-l1* cos(th10)-d2*cos(th20);
y30=-l1* cos(th10)-l2*cos(th20)-d3*cos(th30);
y40=-l1* cos(th10)-l2*cos(th20)-l3*cos(th30)-d4*cos(th40);

%Define initial velocities
x1d0=0; y1d0=0;
x2d0=0; y2d0=0;
x3d0=0; y3d0=0;
x4d0=0; y4d0=0;

z0= [x10; y10; th10; x20; y20; th20; x30; y30; th30; x40; y40; th40; x1d0; y1d0; th1d0; x2d0; y2d0; th2d0; x3d0; y3d0; th3d0; x4d0; y4d0; th4d0];

small = 1e-8;
options = odeset('RelTol', small , 'AbsTol', small );
f = @(t,z)my4barlinkrhsDAE(t,z,p);
[tarray, zarray] = ode45(f,tspan, z0,options);

x1array=zarray(:,1); y1array=zarray(:,2);
th1array = zarray(:,3);
x2array=zarray(:,4); y2array=zarray(:,5);
th2array = zarray(:,6);
x3array=zarray(:,7); y3array=zarray(:,8);
th3array = zarray(:,9);
x4array=zarray(:,10); y4array=zarray(:,11);
th4array = zarray(:,12);
x1darray=zarray(:,13); y1darray=zarray(:,14);
th1darray=zarray(:,15);
x2darray=zarray(:,16); y2darray=zarray(:,17);
th2darray=zarray(:,18);
x3darray=zarray(:,19); y3darray=zarray(:,20);
th3darray=zarray(:,21);
x4darray=zarray(:,22); y4darray=zarray(:,23);
th4darray=zarray(:,24);

%Unpack variables
l1 = p.l1;  l2 = p.l2; l3=p.l3; l4=p.l4;

%Convert to cartesian co-ordinates
x= [zeros(ntimes,1),  l1*sin(th1array),  l1*sin(th1array)  + l2*sin(th2array), l1*sin(th1array)  + l2*sin(th2array) + l3*sin(th3array), l1*sin(th1array)  + l2*sin(th2array) + l3*sin(th3array)+l4*sin(th4array)];
y= [zeros(ntimes,1), -l1*cos(th1array), -l1* cos(th1array) - l2*cos(th2array), -l1* cos(th1array) - l2*cos(th2array)- l3*cos(th3array), -l1* cos(th1array) - l2*cos(th2array)- l3*cos(th3array)-l4*cos(th4array)];

%Plot 4 bar linkage
for i = 1:ntimes
    figure(1)
    plot(x(i,:),y(i,:),'LineWidth', 3)
    xlabel('x-axis');
    ylabel('y-axis');
    title('4-bar linkage');
    axis equal
    axis([-2 2 -2 2])
    pause(.01)
    shg
end

%Store as array
n=4; %Number of links
d=[d1 d2 d3 d4];
l=[l1 l2 l3 l4];
I=[I1 I2 I3 I4];
m=[m1 m2 m3 m4];

%Derivation of Total Energy of System
total_energy=[];
for f=1:ntimes
    th=zarray(f,3:3:3*n);
    thd=zarray(f,15:3:6*n);
    r=sym('er',[n,3]); eth=sym('eth',[n,3]);
    rgo=sym('rgo',[n,3]);
    rao=sym('rao',[n,3]);
    
    i=[1 0 0];
    j=[0 1 0];
    k=cross(i,j);
    
    templ=0; rao(1,:)=0;
    tempv=0; KE=0; PE=0;
    tempvel=0; tempPE=0; vel=zeros(n,1); v=zeros(n,1);
    tempKE=0;v=vel;
    
    for s=1:n
        er(s,:)=i*sin(th(s))-j*cos(th(s));
        eth(s,:)=cross(k,er(s,:));
        
        rga(s,:)=d(s)*er(s,:);
        rba(s,:)=l(s)*er(s,:);
        rbg=rba-rga;
        
        templ=templ+rba(s,:);
        rgo(s,:)=rga(s,:)+templ-rba(s,:);
        rao(s,:)=templ;
        
        tempv=tempv+rba(s,:).*thd(s);
        vgo(s,:)=rga(s,:).*thd(s)+tempv-rba(s,:).*thd(s);
        vao(s,:)=tempv;
        
        for sz=1:2
            tempvel=vgo(s,sz).*vgo(s,sz);
            vel(s)=vel(s)+tempvel;
        end
        vel(s)=vel(s)-v(s);
        tempKE=0.5*m(s)*vel(s)+0.5*I(s)*thd(s).^2;
        KE=KE+tempKE;
        
        tempPE=m(s)*g*rgo(s,2);
        PE=PE+tempPE;
        
    end
    
    total_energy(f)=KE+PE;
end

figure;
hold on
plot(tspan,total_energy-total_energy(1));
xlabel('Time');
ylabel('Change in total energy - (E-E(1))');
title('Energy conservation check');
toc
