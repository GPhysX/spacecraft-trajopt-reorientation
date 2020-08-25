clear all
clc
close all

syms m1 m2 M J Ih Jh x0 y0 t

Ih=[1;0;0]; Jh=[0;1;0]; Kh=[0;0;1];

x1  = str2sym("x1(t)");
y1  = str2sym("y1(t)");
x2  = str2sym("x2(t)");
y2  = str2sym("y2(t)");
th  = str2sym("th(t)");

ih  = Ih*cos(th)+Jh*sin(th);
jh  = -Ih*sin(th)+Jh*cos(th);

rm1G = (x1*ih+y1*jh+x0*Ih+y0*Jh); % position vector of m in G
rm2G = (x2*ih+y2*jh+x0*Ih+y0*Jh); % position vector of m in G
rMG = (x0*Ih+y0*Jh);             % position vector of c.m of M in G

equ_cm = simplify(m1*(rm1G)+M*(rMG)+m2*(rm2G));

x0y0=solve(equ_cm==0,{x0,y0}); % assume system center of mass is at the origin of G

x0sol  = x0y0.x0; % x-coordinate of center of mass of M in G
y0sol  = x0y0.y0; % y-coordinate of the center of mass of M in G

rm1G   = subs(rm1G,{x0,y0},{x0sol,y0sol}); 
rm2G = subs(rm2G,{x0,y0},{x0sol,y0sol}); 
rMG   = subs(rMG,{x0,y0},{x0sol,y0sol});

Hm = simplify(m1*cross(rm1G,diff(rm1G,t))+m2*cross(rm2G,diff(rm2G,t)) ); %angular momentum of m about system c.m. in G
HM = simplify(J*diff(th,t)*Kh+M*cross(rMG,diff(rMG,t))); %angular momentum of M about system c.m. in G

syms x1 y1 x2 y2 th x1d y1d x2d y2d thetad
equ1=subs(simplify(Hm(3)+HM(3)),...
    {str2sym('x1(t)'),str2sym('y1(t)'),...
     str2sym('x2(t)'),str2sym('y2(t)'),str2sym('th(t)'),...
     str2sym('diff(x1(t), t)'),str2sym('diff(y1(t), t)'),...
     str2sym('diff(x2(t), t)'),str2sym('diff(y2(t), t)')...
     str2sym('diff(th(t), t)')},[x1,y1,x2,y2,th,x1d,y1d,x2d,y2d,thetad]);

disp('Equation for thetadot:');
collect(simplify(solve(equ1,thetad)),[x1d,y1d,x2d,y2d])

