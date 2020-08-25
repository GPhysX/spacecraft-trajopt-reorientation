import casadi.*

%-------------------------------------------------------------------------
% Parameters to change for different cases
%-------------------------------------------------------------------------

shape = 1; % Shape?
% 1 - circle
% 2 - square
% 3 - rectangle
% 4 - bowtie

fix_tf = 1; % Fixed end time?
% 0 - End time not fixed
% 1 - Fixed end time

orig = 1; % Begin and return to origin?
% 1 - Yes
% 0 - No

num_meth = 1;
% 0 - FE
% 1 - Trapezoidal
% 2 - RK2
% 3 - RK4

restricted = 0;
% 1 - masses restricted to separate sides
% 0 - masses are allowed to both be on the same side

tf = 80; % End time (if fixed)
N = 400; % Number of intervals

thf = pi/6; % Terminal angle

maxcontrolmag = 1/2; % Max. magnitude of control
%-------------------------------------------------------------------------

min_control = 1;
if fix_tf == 0
   min_control = 0; 
end

global m1 m2 M J
if shape == 1
    R = 1/sqrt(pi); % Spacecraft radius/length for 1 & 2 
else
    R = 1;
end
width = sqrt(2)/2; % Spacecraft width for 3
length = sqrt(2); % Spacecraft length for 3
side = 1; % Spacecraft side length for 4

% Masses of internal mass and spacecraft
m1 = 1/2;
m2 = 1/2;
M = 10;
if shape == 1
   J = 0.5*M*R^2; % Moment of inertia of a disc
elseif shape == 2
   J = (1/6)*M*R^2; % Moment of inertia of a square
elseif shape == 3
   J = (1/6)*M*(length^2+width^2); % Moment of inertia of rectangle
elseif shape == 4
   J = (5/6)*M*side^2;
end
T = MX.sym('T');
dt = T/N;

% State and control vars.
nx = 5;
x1 = SX.sym('x1');
y1 = SX.sym('y1');
x2 = SX.sym('x2');
y2 = SX.sym('y2');
th = SX.sym('th');
x = [x1; y1; x2; y2; th];
x1d = SX.sym('x1d');
y1d = SX.sym('y1d');
x2d = SX.sym('x2d');
y2d = SX.sym('y2d');
u = [x1d; y1d; x2d; y2d];

xdot = spacecraft_model(x, u);

f = Function('f', {x, u}, {xdot});

% Initial conditions
x10 = 0;
y10 = 0;
x20 = 0;
y20 = 0;
th0 = 0;

% Initialize NLP
w = {};
w0 = [];
lbw = [];
ubw = [];
g = {};
lbg = [];
ubg = [];

% Final time
w = {w{:}, T};
if fix_tf == 0
    lbw = [lbw; 0];
    ubw = [ubw; inf];
else
    lbw = [lbw; tf];
    ubw = [ubw; tf];
end
w0 = [w0; tf];

% Setting initial conditions
Xk = MX.sym('X0', nx);
w = {w{:}, Xk};
if orig == 0
lbw = [lbw; -inf; -inf; -inf; -inf; th0];
ubw = [ubw; inf; inf; inf; inf; th0];
else
    lbw = [lbw; x10; y10; x20; y20; th0];
    ubw = [ubw; x10; y10; x20; y20; th0];    
end
w0 = [w0; 0; 0; 0; 0; 0];

O = 0;

% NLP formulation
for k = 0:N-1
    
   U1k = MX.sym(['U1_' num2str(k)]);
   U2k = MX.sym(['U2_' num2str(k)]);
   U3k = MX.sym(['U3_' num2str(k)]);
   U4k = MX.sym(['U4_' num2str(k)]);
   w = {w{:}, U1k, U2k, U3k, U4k};
   if min_control == 0
      lbw = [lbw; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag];
      ubw = [ubw; maxcontrolmag; maxcontrolmag; maxcontrolmag; maxcontrolmag];
   else
      lbw = [lbw; -inf; -inf; -inf; -inf];
      ubw = [ubw; inf; inf; inf; inf]; 
   end
   w0 = [w0; 0; 0; 0; 0];
   
   Xkp1 = MX.sym(['X_' num2str(k+1)], nx);
   w = [w, {Xkp1}];
   lbw = [lbw; -inf; -inf; -inf; -inf; 0];
   ubw = [ubw; inf; inf; inf; inf; 2*pi];
   w0 = [w0; 0; 0; 0; 0; 0];
   
   
   if num_meth == 0
   g = [g, {dt*f(Xk, [U1k; U2k; U3k; U4k])+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0];
   elseif num_meth == 1
   g = [g, {dt*( 0.5*f(Xk, [U1k; U2k; U3k; U4k]) + 0.5*f(Xkp1, [U1k; U2k; U3k; U4k]))+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0];   
   elseif num_meth == 2
   Xmid = Xk + 0.5*dt*f(Xk, [U1k; U2k; U3k; U4k]);
   g = [g, {dt*f(Xmid, [U1k; U2k; U3k; U4k])+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0];
   elseif num_meth == 3
   f0 = f(Xk, [U1k; U2k; U3k; U4k]);
   f1 = f(Xk+0.5*dt*f0, [U1k; U2k; U3k; U4k]);
   f2 = f(Xk+0.5*dt*f1, [U1k; U2k; U3k; U4k]);
   f3 = f(Xk+dt*f2, [U1k; U2k; U3k; U4k]);
   g = [g, {(1/6)*dt*(f0+2*f1+2*f2+f3)+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0];
   end
   
   if shape == 1
      g = [g, {(Xk(1)^2+Xk(2)^2)}];
      lbg = [lbg; 0];
      ubg = [ubg; R^2];
      g = [g, {(Xk(3)^2+Xk(4)^2)}];
      lbg = [lbg; 0];
      ubg = [ubg; R^2];
   elseif shape == 2
      g = [g, {Xk(1)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(2)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(3)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(4)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
   elseif shape == 3
      g = [g, {Xk(1)}];
      lbg = [lbg; -width/2];
      ubg = [ubg; width/2];
      g = [g, {Xk(2)}];
      lbg = [lbg; -length/2];
      ubg = [ubg; length/2];
      g = [g, {Xk(3)}];
      lbg = [lbg; -width/2];
      ubg = [ubg; width/2];
      g = [g, {Xk(4)}];
      lbg = [lbg; -length/2];
      ubg = [ubg; length/2];
   elseif shape == 4
      g = [g, {100*Xk(2)^2-(100*Xk(1)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      
      g = [g, {100*Xk(4)^2-(100*Xk(3)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      
      if restricted == 1
      
      g = [g, {Xk(1)}];
      lbg = [lbg; 0];
      ubg = [ubg; side*sqrt(3)/2];
      
      g = [g, {Xk(3)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; 0];
      
      else
      
      g = [g, {Xk(1)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; side*sqrt(3)/2];
      
      g = [g, {Xk(3)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; side*sqrt(3)/2];
          
      if k > 1
            scle = [0, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, 1]';
            slpe1 = Xk(1:2)-Xkn1(1:2);
            xp1 = slpe1(1)*scle+Xkn1(1);
            yp1 = slpe1(2)*scle+Xkn1(2);
            g = [g, {100*yp1.^2-100*(1/3)*xp1.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            
            slpe2 = Xk(3:4)-Xkn1(3:4);
            xp2 = slpe2(1)*scle+Xkn1(1);
            yp2 = slpe2(2)*scle+Xkn1(2);
            g = [g, {100*yp2.^2-100*(1/3)*xp2.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
      end
          
      end
   end
   
   if min_control == 1
      O = O+(U1k^2+U2k^2+U3k^2+U4k^2)*1/2;
   end
   
   Xkn1 = Xk;
   Xk = Xkp1;
    
end

% Terminal conditions
ubw(end-8:end-5) = [0; 0; 0; 0];
if orig == 1
lbw(end-4:end) = [0; 0; 0; 0; thf];
ubw(end-4:end) = [0; 0; 0; 0; thf];
else
    lbw(end) = thf;
    ubw(end) = thf;
end
g = [g, {f(Xkn1, [U1k; U2k; U3k; U4k])}];
lbg = [lbg; 0; 0; 0; 0; 0];
ubg = [ubg; 0; 0; 0; 0; 0];
g = [g, {[U1k; U2k; U3k; U4k]}];
lbg = [lbg; 0; 0; 0; 0];
ubg = [ubg; 0; 0; 0; 0];

% Cost function
if fix_tf == 0
   O = T; 
end

% Supply gradient (not used in this example)
% dJdw = gradient(O,vertcat(w{:}));
% dGdw = jacobian(vertcat(g{:}),vertcat(w{:}));

options = struct;
options.ipopt.max_iter = 50000;

% Create an NLP solver
prob = struct('f',  O, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob, options);

% Solve the NLP
sol = solver('x0', w0, ...
            'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%% Plots

T_opt = w_opt(1);
t = linspace(0, T_opt, N+1);
x1_opt = w_opt(2:9:end);
y1_opt = w_opt(3:9:end);
x2_opt = w_opt(4:9:end);
y2_opt = w_opt(5:9:end);
th_opt = w_opt(6:9:end);
x1d_opt = w_opt(7:9:end);
y1d_opt = w_opt(8:9:end);
x2d_opt = w_opt(9:9:end);
y2d_opt = w_opt(10:9:end);
subplot(3, 2, 1);
plot(x1_opt, y1_opt, 'b-');
hold on
plot(x2_opt, y2_opt, 'g-');
if shape == 1
    sgtitle('Circular spacecraft reorientation', 'Interpreter', 'latex',...
    'FontSize', 14);
    theta = linspace(0, 2*pi, 100);
    plot(R*sin(theta), R*cos(theta), 'r--');
elseif shape == 2
    sgtitle('Square spacecraft reorientation', 'Interpreter', 'latex',...
    'FontSize', 14);
    plot(-R/2*ones(2), [-R/2 R/2], 'r--');
    plot(R/2*ones(2), [-R/2 R/2], 'r--');
    plot([-R/2 R/2], -R/2*ones(2), 'r--');
    plot([-R/2 R/2], R/2*ones(2), 'r--');
elseif shape == 3
    sgtitle('Rectangular spacecraft reorientation', 'Interpreter',...
        'latex', 'FontSize', 14);
    plot(-width/2*ones(2), [-length/2 length/2], 'r--');
    plot(width/2*ones(2), [-length/2 length/2], 'r--');
    plot([-width/2 width/2], -length/2*ones(2), 'r--');
    plot([-width/2 width/2], length/2*ones(2), 'r--');
elseif shape == 4
    x = linspace(-side*sqrt(3)/2, side*sqrt(3)/2, 3);
    y1 = -1/sqrt(3) * x;
    y2 = 1/sqrt(3) * x;
    sgtitle('Bowtie spacecraft reorientation',...
        'Interpreter', 'latex', 'FontSize', 14);
    plot(x, y1, 'r--');
    plot(x, y2, 'r--');
    plot(-side*sqrt(3)/2*ones(2), [-side/2 side/2], 'r--');
    plot(side*sqrt(3)/2*ones(2), [-side/2 side/2], 'r--');
end
hold off
xlabel('$x_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$','$M$', 'Interpreter', 'latex', 'Location',...
    'northeastoutside');
axis equal;
subplot(3, 2, 2);
plot(t, th_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\theta$, rad', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 2, 3);
plot(t, x1_opt);
hold on
plot(t, x2_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', 'Interpreter', 'latex');
subplot(3, 2, 4);
plot(t, y1_opt);
hold on
plot(t, y2_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', 'Interpreter', 'latex');
subplot(3, 2, 5)
plot(t(1:end-1), x1d_opt);
hold on
plot(t(1:end-1), y1d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{1}$', '$\dot{y}_{1}$',...
    'Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 2, 6);
plot(t(1:end-1), x2d_opt);
hold on
plot(t(1:end-1), y2d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{2}$', '$\dot{y}_{2}$',...
    'Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%% Plots for Report

figure(1);
plot(x1_opt, y1_opt, 'b-');
hold on
plot(x2_opt, y2_opt, 'g-');
if shape == 1
    sgtitle('Circular spacecraft reorientation', 'Interpreter', 'latex',...
    'FontSize', 14);
    theta = linspace(0, 2*pi, 100);
    plot(R*sin(theta), R*cos(theta), 'b--');
elseif shape == 2
    sgtitle('Square spacecraft reorientation', 'Interpreter', 'latex',...
    'FontSize', 14);
    plot(-R/2*ones(2), [-R/2 R/2], 'b--');
    plot(R/2*ones(2), [-R/2 R/2], 'b--');
    plot([-R/2 R/2], -R/2*ones(2), 'b--');
    plot([-R/2 R/2], R/2*ones(2), 'b--');
elseif shape == 3
    sgtitle('Rectangular spacecraft reorientation', 'Interpreter',...
        'latex', 'FontSize', 14);
    plot(-width/2*ones(2), [-length/2 length/2], 'r--');
    plot(width/2*ones(2), [-length/2 length/2], 'r--');
    plot([-width/2 width/2], -length/2*ones(2), 'r--');
    plot([-width/2 width/2], length/2*ones(2), 'r--');
elseif shape == 4
    x = linspace(-side*sqrt(3)/2, side*sqrt(3)/2, 3);
    y1 = -1/sqrt(3) * x;
    y2 = 1/sqrt(3) * x;
    sgtitle('Bowtie spacecraft reorientation',...
        'Interpreter', 'latex', 'FontSize', 14);
    plot(x, y1, 'r--');
    plot(x, y2, 'r--');
    plot(-side*sqrt(3)/2*ones(2), [-side/2 side/2], 'r--');
    plot(side*sqrt(3)/2*ones(2), [-side/2 side/2], 'r--');
end
hold off
xlabel('$x_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', '$M$', 'Interpreter', 'latex', 'Location',...
    'northeastoutside');
axis equal;

figure(2);
plot(t, th_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\theta$, rad', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

figure(3);
subplot(2, 1, 1);
plot(t, x1_opt);
hold on
plot(t, x2_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', 'Interpreter', 'latex','Location',...
    'northeastoutside');
subplot(2, 1, 2);
plot(t, y1_opt);
hold on
plot(t, y2_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$, m', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-1 1]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', 'Interpreter', 'latex','Location',...
    'northeastoutside');

figure(4);
subplot(2, 1, 1);
plot(t(1:end-1), x1d_opt);
hold on
plot(t(1:end-1), y1d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{1}$', '$\dot{y}_{1}$',...
    'Location',...
    'northeastoutside','Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(2, 1, 2);
plot(t(1:end-1), x2d_opt);
hold on
plot(t(1:end-1), y2d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{2}$', '$\dot{y}_{2}$',...
    'Interpreter', 'latex','FontSize', 12,'Location',...
    'northeastoutside');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);


%% Model

function [dxdt] = spacecraft_model(x, u)
global m1 m2 M J
% State and control vars.
x1 = x(1);
y1 = x(2);
x2 = x(3);
y2 = x(4);
th = x(5);
x1d = u(1);
y1d = u(2);
x2d = u(3);
y2d = u(4);

% Dynamics
thd = ((M*m1*y1 + m1*m2*y1 - m1*m2*y2)/(J*m1 + J*m2 + J*M + m1*m2*x1^2 ...
      + m1*m2*x2^2 + m1*m2*y1^2 + m1*m2*y2^2 + M*m1*x1^2 + M*m2*x2^2 ...
      + M*m1*y1^2 + M*m2*y2^2 - 2*m1*m2*x1*x2 - 2*m1*m2*y1*y2))*x1d ...
      + (-(M*m1*x1 + m1*m2*x1 - m1*m2*x2)/(J*m1 + J*m2 + J*M + ...
      m1*m2*x1^2 + m1*m2*x2^2 + m1*m2*y1^2 + m1*m2*y2^2 + M*m1*x1^2 ...
      + M*m2*x2^2 + M*m1*y1^2 + M*m2*y2^2 - 2*m1*m2*x1*x2 - ...
      2*m1*m2*y1*y2))*y1d + ((M*m2*y2 - m1*m2*y1 + m1*m2*y2)/(J*m1 ...
      + J*m2 + J*M + m1*m2*x1^2 + m1*m2*x2^2 + m1*m2*y1^2 + m1*m2*y2^2 ...
      + M*m1*x1^2 + M*m2*x2^2 + M*m1*y1^2 + M*m2*y2^2 - 2*m1*m2*x1*x2 ...
      - 2*m1*m2*y1*y2))*x2d + (-(M*m2*x2 - m1*m2*x1 + m1*m2*x2)/(J*m1 + ...
      J*m2 + J*M + m1*m2*x1^2 + m1*m2*x2^2 + m1*m2*y1^2 + m1*m2*y2^2 ...
      + M*m1*x1^2 + M*m2*x2^2 + M*m1*y1^2 + M*m2*y2^2 - 2*m1*m2*x1*x2 ...
      - 2*m1*m2*y1*y2))*y2d;
dxdt = [x1d;y1d;x2d;y2d;thd];

end