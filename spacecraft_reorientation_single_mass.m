import casadi.*

%-------------------------------------------------------------------------
% Parameters to change for different cases
%-------------------------------------------------------------------------

shape = 1; % Shape?
% 1 - circle
% 2 - square
% 3 - rectangle
% 4 - bowtie

right = 2; % Side for bowtie?
% 0 - Left side
% 1 - Right side
% 2 - Both

orig = 1; % Begin and return to origin?
% 1 - Yes
% 0 - No

num_meth = 1;
% 0 - FE
% 1 - Trapezoidal
% 2 - RK2
% 3 - RK4

fix_tf = 1; % Fixed end time?
% 0 - End time not fixed
% 1 - Fixed end time

tf = 80; % End time (if fixed)
N = 200; % Number of intervals

thf = pi/6; % Terminal angle

maxcontrolmag = 1/2; % Max. magnitude of control
%-------------------------------------------------------------------------

min_control = 1;
if fix_tf == 0
   min_control = 0; 
end

global m M J
if shape == 1
    R = 1/sqrt(pi); % Spacecraft radius/length for 1 & 2 
else
    R = 1;
end
width = sqrt(2)/2; % Spacecraft width for 3
length = sqrt(2); % Spacecraft length for 3
side = 1; % Spacecraft side length for 4

% Masses of internal mass and spacecraft
m = 1;
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
nx = 3;
x1 = SX.sym('x1');
y1 = SX.sym('y1');
th = SX.sym('th');
x = [x1; y1; th];
x1d = SX.sym('x1d');
y1d = SX.sym('y1d');
u = [x1d; y1d];

xdot = spacecraft_model(x, u);

f = Function('f', {x, u}, {xdot});

% Initial conditions
x10 = 0;
y10 = 0;
th0 = 0;

% Constraints
x1df = 0;
y1df = 0;

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
lbw = [lbw; -inf; -inf; th0];
ubw = [ubw; inf; inf; th0];
else
    lbw = [lbw; x10; y10; th0];
    ubw = [ubw; x10; y10; th0];    
end
w0 = [w0; 0; 0; 0];

O = 0;

% NLP formulation
for k = 0:N-1
    
   U1k = MX.sym(['U1_' num2str(k)]);
   U2k = MX.sym(['U2_' num2str(k)]);
   w = {w{:}, U1k, U2k};
   if min_control == 0
      lbw = [lbw; -maxcontrolmag; -maxcontrolmag];
      ubw = [ubw; maxcontrolmag; maxcontrolmag];
   else
      lbw = [lbw; -inf; -inf];
      ubw = [ubw; inf; inf]; 
   end

   w0 = [w0; 0; 0];
   
   Xkp1 = MX.sym(['X_' num2str(k+1)], nx);
   w = [w, {Xkp1}];
   lbw = [lbw; -inf; -inf; 0];
   ubw = [ubw; inf; inf; 2*pi];
   w0 = [w0; 0; 0; 0];
   
   if num_meth == 0
   g = [g, {dt*f(Xk, [U1k; U2k])+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0];
   ubg = [ubg; 0; 0; 0];
   elseif num_meth == 1
   g = [g, {dt*( 0.5*f(Xk, [U1k; U2k]) + 0.5*f(Xkp1, [U1k; U2k]))+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0];
   ubg = [ubg; 0; 0; 0];    
   elseif num_meth == 2
   Xmid = Xk + 0.5*dt*f(Xk, [U1k; U2k]);
   g = [g, {dt*f(Xmid, [U1k; U2k])+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0];
   ubg = [ubg; 0; 0; 0];
   elseif num_meth == 3
   f0 = f(Xk, [U1k; U2k]);
   f1 = f(Xk+0.5*dt*f0, [U1k; U2k]);
   f2 = f(Xk+0.5*dt*f1, [U1k; U2k]);
   f3 = f(Xk+dt*f2, [U1k; U2k]);
   g = [g, {(1/6)*dt*(f0+2*f1+2*f2+f3)+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0];
   ubg = [ubg; 0; 0; 0];
   end
   
   if shape == 1
      g = [g, {(Xk(1)^2+Xk(2)^2)}];
      lbg = [lbg; 0];
      ubg = [ubg; R^2];
   elseif shape == 2
      g = [g, {Xk(1)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(2)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
   elseif shape == 3
      g = [g, {Xk(1)}];
      lbg = [lbg; -width/2];
      ubg = [ubg; width/2];
      g = [g, {Xk(2)}];
      lbg = [lbg; -length/2];
      ubg = [ubg; length/2];
   elseif shape == 4
      g = [g, {100*Xk(2)^2-100*(Xk(1)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      g = [g, {Xk(1)}];
      if right == 1
         lbg = [lbg; 0];
         ubg = [ubg; side*sqrt(3)/2];
      elseif right == 0
         lbg = [lbg; -side*sqrt(3)/2];
         ubg = [ubg; 0];
      else
         lbg = [lbg; -side*sqrt(3)/2];
         ubg = [ubg; side*sqrt(3)/2];
         if k > 1
            scle = [0, 1/6, 1/3, 1/2, 2/3, 5/6, 1]';
            slpe = Xk(1:2)-Xkn1(1:2);
            xp = slpe(1)*scle+Xkn1(1);
            yp = slpe(2)*scle+Xkn1(2);
            g = [g, {100*yp.^2-100*(1/3)*xp.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0];
         end
      end
   end
   
   if min_control == 1
      O = O+U1k^2+U2k^2;
   end
   
   Xkn1 = Xk;
   Xk = Xkp1;
    
end

% Terminal conditions
ubw(end-4:end-3) = [0; 0];
if orig == 1
lbw(end-2:end) = [0; 0; thf];
ubw(end-2:end) = [0; 0; thf];
else
    lbw(end) = thf;
    ubw(end) = thf;
end
g = [g, {f(Xkn1, [U1k; U2k])}];
lbg = [lbg; 0; 0; 0];
ubg = [ubg; 0; 0; 0];
g = [g, {[U1k; U2k]}];
lbg = [lbg; 0; 0];
ubg = [ubg; 0; 0];

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
x1_opt = w_opt(2:5:end);
y1_opt = w_opt(3:5:end);
th_opt = w_opt(4:5:end);
x1d_opt = w_opt(5:5:end);
y1d_opt = w_opt(6:5:end);
subplot(3, 2, 1);
plot(x1_opt, y1_opt, 'b-');
hold on
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
    if right == 1
    sgtitle('Bowtie spacecraft reorientation (restricted to right side)',...
        'Interpreter', 'latex', 'FontSize', 14);
    elseif right == 0
    sgtitle('Bowtie spacecraft reorientation (restricted to left side)',...
    'Interpreter', 'latex', 'FontSize', 14);
    else
    sgtitle('Bowtie spacecraft reorientation (no restriction)',...
    'Interpreter', 'latex', 'FontSize', 14);    
    end
    plot(x, y1, 'r--');
    plot(x, y2, 'r--');
    plot(-side*sqrt(3)/2*ones(2), [-side/2 side/2], 'r--');
    plot(side*sqrt(3)/2*ones(2), [-side/2 side/2], 'r--');
end
hold off
xlabel('$x_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m$', '$M$', 'Interpreter', 'latex', 'Location',...
    'northeastoutside');
axis equal;
subplot(3, 2, 2);
plot(t, th_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\theta$, rad', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 2, 3);
plot(t, x1_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 2, 4);
plot(t, y1_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y_{1}$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 2, 5:6)
plot(t(1:end-1), x1d_opt);
hold on
plot(t(1:end-1), y1d_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{1}$', '$\dot{y}_{1}$', 'Interpreter', 'latex',...
        'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%% Plots for Report

figure(1);
plot(x1_opt, y1_opt, 'b-');
hold on
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
legend('$m$', '$M$', 'Interpreter', 'latex', 'Location',...
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
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m$', 'Interpreter', 'latex','Location',...
    'northeastoutside');
subplot(2, 1, 2);
plot(t, y1_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$, m', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-1 1]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m$', 'Interpreter', 'latex','Location',...
    'northeastoutside');

figure(4);
plot(t(1:end-1), x1d_opt);
hold on
plot(t(1:end-1), y1d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}$', '$\dot{y}$',...
    'Location',...
    'northeastoutside','Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%% Model

function [dxdt] = spacecraft_model(x, u)
global m M J
% State and control vars.
x1 = x(1);
y1 = x(2);
th = x(3);
x1d = u(1);
y1d = u(2);

% Dynamics
thd = ((M*m*y1)/(M*m*x1^2 + M*m*y1^2 + J*m + J*M))*x1d...
       + (-(M*m*x1)/(M*m*x1^2 + M*m*y1^2 + J*m + J*M))*y1d;
dxdt = [x1d;y1d;thd];

end