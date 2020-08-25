clear all
close all

syms m M J Ih Jh x0 y0 t

Ih=[1;0;0]; Jh=[0;1;0]; Kh=[0;0;1];

x1  = str2sym("x1(t)");
y1  = str2sym("y1(t)");
th  = str2sym("th(t)");

ih  = Ih*cos(th)+Jh*sin(th);
jh  = -Ih*sin(th)+Jh*cos(th);

rmG = (x1*ih+y1*jh+x0*Ih+y0*Jh); % position vector of m in G
rMG = (x0*Ih+y0*Jh);             % position vector of c.m of M in G

equ_cm = simplify(m*(rmG)+M*(rMG));

x0y0=solve(equ_cm==0,{x0,y0}); % assume system center of mass is at the origin of G

x0sol  = x0y0.x0; % x-coordinate of center of mass of M in G
y0sol  = x0y0.y0; % y-coordinate of the center of mass of M in G

rmG   = subs(rmG,{x0,y0},{x0sol,y0sol}); 
rMG   = subs(rMG,{x0,y0},{x0sol,y0sol});

Hm = simplify(m*cross(rmG,diff(rmG,t)) ); %angular momentum of m about system c.m. in G
HM = simplify(J*diff(th,t)*Kh+M*cross(rMG,diff(rMG,t))); %angular momentum of M about system c.m. in G

syms x1 y1 th x1d y1d thetad
equ1=subs(simplify(Hm(3)+HM(3)),...
    {str2sym('x1(t)'),str2sym('y1(t)'),str2sym('th(t)'),...
     str2sym('diff(x1(t), t)'),str2sym('diff(y1(t), t)'),...
     str2sym('diff(th(t), t)')},[x1,y1,th,x1d,y1d,thetad])

collect(simplify(solve(equ1,thetad)),[x1d,y1d])

global m M J
m=1; M=10; J=10;

function xdot = myode(t,X)
global M m J
x1 = X(1);
y1 = X(2);
th = X(3);
x1d = 0.5*sin(t);
y1d = 0.5*sin(t+pi/4);
thdot=((M*m*y1)/(M*m*x1^2 + M*m*y1^2 + J*m + J*M))*x1d + (-(M*m*x1)/(M*m*x1^2 + M*m*y1^2 + J*m + J*M))*y1d; 
xdot = [x1d;y1d;thdot];

end

