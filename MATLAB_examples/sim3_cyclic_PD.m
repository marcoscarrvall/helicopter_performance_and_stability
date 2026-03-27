%SIMULATION OF A CYCLIC PITCH INPUT THETA_C=1 DEG GIVEN FROM HOVER 
%0.5 SEC<T<1 SEC. Now from the 15th second a PD controller becomes active 
clear
%INITIAL DATA HELICOPTER
g=9.81;	
cla=5.7; %NACA 0012
volh=.075;	%blade solidity	
lok=6;
cds=1.5;
mass=2200;
rho=1.225;
vtip=200;
diam=2*7.32;
iy=10615;
mast=1;
omega=vtip/(diam/2);
area=pi/4*diam^2;
tau=.1;		%time constant in dynamiCs inflow!!!
collect(1)=6*pi/180;
longit(1)=0*pi/180;

%initial values;
t0=0;
u0=0;
w0=0;
q0=0;
pitch0=0*pi/180;
x0=0;
labi0=sqrt(mass*g/(area*2*rho))/vtip;

t(1)=t0;
u(1)=u0;
w(1)=w0;
q(1)=q0;
pitch(1)=pitch0;
x(1)=x0;
labi(1)=labi0;
z(1)=0;

%INTEGRATION 
aantal=800;
teind=80;
stap=(teind-t0)/aantal;

for i=1:aantal 
   if t(i)>=0.5 & t(i)<=1 longit(i)=1*pi/180;
   else longit(i)=0*pi/180;
   end
    
   if t(i)>=15 longitgrd(i)=.2*pitch(i)*180/pi+.2*q(i)*180/pi;%PD in deg
       longit(i)=longitgrd(i)*pi/180;	%in rad
   end    
   %longit(i)=longitgrd(i)*pi/180;	%in rad
   
%NO LAW FOR COLLECTIVE

c(i)=u(i)*sin(pitch(i))-w(i)*cos(pitch(i));
h(i)=-z(i);
collect(i)=collect(1);

%Defining the differential equations

%defining the nondimensional notations
qdiml(i)=q(i)/omega;
vdiml(i)=sqrt(u(i)^2+w(i)^2)/vtip;
if u(i)==0 	if w(i)>0 	phi(i)=pi/2;
        else phi(i)=-pi/2;end
else
phi(i)=atan(w(i)/u(i));
end
if u(i)<0
phi(i)=phi(i)+pi;
end
alfc(i)=longit(i)-phi(i);

mu(i)=vdiml(i)*cos(alfc(i));
labc(i)=vdiml(i)*sin(alfc(i));

%a1 Flapping calculi
teller(i)=-16/lok*qdiml(i)+8/3*mu(i)*collect(i)-2*mu(i)*(labc(i)+labi(i));
a1(i)=teller(i)/(1-.5*mu(i)^2);

%the thrust coefficient
ctelem(i)=cla*volh/4*(2/3*collect(i)*(1+1.5*mu(i)^2)-(labc(i)+labi(i)));
%Thrust coefficient from Glauert
alfd(i)=alfc(i)-a1(i);
ctglau(i)=2*labi(i)*sqrt((vdiml(i)*cos(alfd(i)))^2+(vdiml(i)*...
sin(alfd(i))+labi(i))^2);

%Equations of motion
labidot(i)=ctelem(i); 
thrust(i)=labidot(i)*rho*vtip^2*area;
helling(i)=longit(i)-a1(i);
vv(i)=vdiml(i)*vtip; 		%it is 1/sqrt(u^2+w^2)

udot(i)=-g*sin(pitch(i))-cds/mass*.5*rho*u(i)*vv(i)+...
thrust(i)/mass*sin(helling(i))-q(i)*w(i);

wdot(i)=g*cos(pitch(i))-cds/mass*.5*rho*w(i)*vv(i)-...
thrust(i)/mass*cos(helling(i))+q(i)*u(i);

qdot(i)=-thrust(i)*mast/iy*sin(helling(i));

pitchdot(i)=q(i);

xdot(i)=u(i)*cos(pitch(i))+w(i)*sin(pitch(i));

zdot(i)=-c(i);

labidot(i)=(ctelem(i)-ctglau(i))/tau;
%corrdot(i)=uwens-u(i);
%corrcdot(i)=cwens(i)-c(i);

u(i+1)=u(i)+stap*udot(i);
w(i+1)=w(i)+stap*wdot(i);
q(i+1)=q(i)+stap*qdot(i);
pitch(i+1)=pitch(i)+stap*pitchdot(i);
x(i+1)=x(i)+stap*xdot(i);
labi(i+1)=labi(i)+stap*labidot(i);
z(i+1)=z(i)+stap*zdot(i);
t(i+1)=t(i)+stap;
end;

plot(t,u),xlabel('t (s)'),ylabel('u(m)'),grid,pause;
plot(t,pitch*180/pi),xlabel('t (s)'),ylabel('pitch(deg)'),grid,pause;
plot(t,x),xlabel('t (s)'),ylabel('x(m)'),grid,pause;
plot(t,w),xlabel('t (s)'),ylabel('w(m)'),grid,pause;
plot(t,q),xlabel('t (s)'),ylabel('q(m)'),grid,pause; 
plot(t,labi),xlabel('t (s)'),ylabel('labi(m)'),grid,pause;
plot(t,-z),xlabel('t (s)'),ylabel('h(m)'),grid,pause;
plot(t(1:800),longit*180/pi),xlabel('t (s)'),ylabel('longit grd'),grid,pause;


