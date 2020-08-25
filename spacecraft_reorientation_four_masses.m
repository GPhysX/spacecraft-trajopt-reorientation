import casadi.*

%-------------------------------------------------------------------------
% Parameters to change for different cases
%-------------------------------------------------------------------------

shape = 4; % Shape?
% 1 - circle
% 2 - square
% 3 - rectangle
% 4 - bowtie

fix_tf = 1; % Fixed end time?
% 0 - End time not fixed
% 1 - Fixed end time

orig = 0; % Begin and return to origin?
% 1 - Yes
% 0 - No

num_meth = 1;
% 0 - FE
% 1 - Trapezoidal
% 2 - RK2
% 3 - RK4

tf = 80; % End time (if fixed)
N = 300; % Number of intervals

thf = pi/6; % Terminal angle

maxcontrolmag = 1/2; % Max. magnitude of control
%-------------------------------------------------------------------------

min_control = 1;
if fix_tf == 0
   min_control = 0; 
end

global m1 m2 m3 m4 M J
if shape == 1
    R = 1/sqrt(pi); % Spacecraft radius/length for 1 & 2 
else
    R = 1;
end
width = sqrt(2)/2; % Spacecraft width for 3
length = sqrt(2); % Spacecraft length for 3
side = 1; % Spacecraft side length for 4

% Masses of internal mass and spacecraft
m1 = 1/4;
m2 = 1/4;
m3 = 1/4;
m4 = 1/4;
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
nx = 9;
x1 = SX.sym('x1');
y1 = SX.sym('y1');
x2 = SX.sym('x2');
y2 = SX.sym('y2');
x3 = SX.sym('x3');
y3 = SX.sym('y3');
x4 = SX.sym('x4');
y4 = SX.sym('y4');
th = SX.sym('th');
x = [x1; y1; x2; y2; x3; y3; x4; y4; th];
x1d = SX.sym('x1d');
y1d = SX.sym('y1d');
x2d = SX.sym('x2d');
y2d = SX.sym('y2d');
x3d = SX.sym('x3d');
y3d = SX.sym('y3d');
x4d = SX.sym('x4d');
y4d = SX.sym('y4d');
u = [x1d; y1d; x2d; y2d; x3d; y3d; x4d; y4d];

xdot = spacecraft_model(x, u);

f = Function('f', {x, u}, {xdot});

% Initial conditions
x10 = 0;
y10 = 0;
x20 = 0;
y20 = 0;
x30 = 0;
y30 = 0;
x40 = 0;
y40 = 0;
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
lbw = [lbw; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; th0];
ubw = [ubw; inf; inf; inf; inf; inf; inf; inf; inf; th0];
else
    lbw = [lbw; x10; y10; x20; y20; x30; y30; x40; y40; th0];
    ubw = [ubw; x10; y10; x20; y20; x30; y30; x40; y40; th0];    
end
w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

O = 0;

% NLP formulation
for k = 0:N-1
    
   U1k = MX.sym(['U1_' num2str(k)]);
   U2k = MX.sym(['U2_' num2str(k)]);
   U3k = MX.sym(['U3_' num2str(k)]);
   U4k = MX.sym(['U4_' num2str(k)]);
   U5k = MX.sym(['U5_' num2str(k)]);
   U6k = MX.sym(['U6_' num2str(k)]);
   U7k = MX.sym(['U7_' num2str(k)]);
   U8k = MX.sym(['U8_' num2str(k)]);
   w = {w{:}, U1k, U2k, U3k, U4k, U5k, U6k, U7k, U8k};
   if min_control == 0
      lbw = [lbw; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag; -maxcontrolmag];
      ubw = [ubw; maxcontrolmag; maxcontrolmag; maxcontrolmag; maxcontrolmag; maxcontrolmag; maxcontrolmag; maxcontrolmag; maxcontrolmag];
   else
      lbw = [lbw; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
      ubw = [ubw; inf; inf; inf; inf; inf; inf; inf; inf]; 
   end
   w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0];
   
   Xkp1 = MX.sym(['X_' num2str(k+1)], nx);
   w = [w, {Xkp1}];
   lbw = [lbw; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; 0];
   ubw = [ubw; inf; inf; inf; inf; inf; inf; inf; inf; 2*pi];
   w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   
   if num_meth == 0
   g = [g, {dt*f(Xk, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k])+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   elseif num_meth == 1
   g = [g, {dt*( 0.5*f(Xk, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]) + 0.5*f(Xkp1, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]))+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];  
   elseif num_meth == 2
   Xmid = Xk + 0.5*dt*f(Xk, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]);
   g = [g, {dt*f(Xmid, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k])+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   elseif num_meth == 3
   f0 = f(Xk, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]);
   f1 = f(Xk+0.5*dt*f0, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]);
   f2 = f(Xk+0.5*dt*f1, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]);
   f3 = f(Xk+dt*f2, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]);
   g = [g, {(1/6)*dt*(f0+2*f1+2*f2+f3)+Xk-Xkp1}];
   lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
   end
   
   if shape == 1
      g = [g, {(Xk(1)^2+Xk(2)^2)}];
      lbg = [lbg; 0];
      ubg = [ubg; R^2];
      g = [g, {(Xk(3)^2+Xk(4)^2)}];
      lbg = [lbg; 0];
      ubg = [ubg; R^2];
      g = [g, {(Xk(5)^2+Xk(6)^2)}];
      lbg = [lbg; 0];
      ubg = [ubg; R^2];
      g = [g, {(Xk(7)^2+Xk(8)^2)}];
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
      g = [g, {Xk(5)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(6)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(7)}];
      lbg = [lbg; -R/2];
      ubg = [ubg; R/2];
      g = [g, {Xk(8)}];
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
      g = [g, {Xk(5)}];
      lbg = [lbg; -width/2];
      ubg = [ubg; width/2];
      g = [g, {Xk(6)}];
      lbg = [lbg; -length/2];
      ubg = [ubg; length/2];
      g = [g, {Xk(7)}];
      lbg = [lbg; -width/2];
      ubg = [ubg; width/2];
      g = [g, {Xk(8)}];
      lbg = [lbg; -length/2];
      ubg = [ubg; length/2];
   elseif shape == 4
      g = [g, {100*Xk(2)^2-100*(Xk(1)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      
      g = [g, {100*Xk(4)^2-100*(Xk(3)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      
      g = [g, {100*Xk(6)^2-100*(Xk(5)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      
      g = [g, {100*Xk(8)^2-100*(Xk(7)^2)/3}];
      lbg = [lbg; -inf];
      ubg = [ubg; 0];
      
      g = [g, {Xk(1)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; side*sqrt(3)/2];
      
      g = [g, {Xk(3)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; side*sqrt(3)/2];
      
      g = [g, {Xk(5)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; side*sqrt(3)/2];
      
      g = [g, {Xk(7)}];
      lbg = [lbg; -side*sqrt(3)/2];
      ubg = [ubg; side*sqrt(3)/2];
          
      if k > 1
            scle = [0, 1/11, 2/11, 3/11, 4/11, 5/11, 6/11, 7/11, 8/11, 9/11, 10/11, 1]';
            slpe1 = Xk(1:2)-Xkn1(1:2);
            xp1 = slpe1(1)*scle+Xkn1(1);
            yp1 = slpe1(2)*scle+Xkn1(2);
            g = [g, {100*yp1.^2-100*(1/3)*xp1.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            
            slpe2 = Xk(3:4)-Xkn1(3:4);
            xp2 = slpe2(1)*scle+Xkn1(1);
            yp2 = slpe2(2)*scle+Xkn1(2);
            g = [g, {100*yp2.^2-100*(1/3)*xp2.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            
            slpe3 = Xk(5:6)-Xkn1(5:6);
            xp3 = slpe3(1)*scle+Xkn1(1);
            yp3 = slpe3(2)*scle+Xkn1(2);
            g = [g, {100*yp3.^2-100*(1/3)*xp3.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            
            slpe4 = Xk(7:8)-Xkn1(7:8);
            xp4 = slpe4(1)*scle+Xkn1(1);
            yp4 = slpe4(2)*scle+Xkn1(2);
            g = [g, {100*yp4.^2-100*(1/3)*xp4.^2}];
            lbg = [lbg; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf; -inf];
            ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
      end
          
   end
   
   if min_control == 1
      O = O+(U1k^2+U2k^2+U3k^2+U4k^2+U5k^2+U6k^2+U7k^2+U8k^2)*1/4;
   end
   
   Xkn1 = Xk;
   Xk = Xkp1;
    
end

% Terminal conditions
ubw(end-16:end-9) = [0; 0; 0; 0; 0; 0; 0; 0];
if orig == 1
lbw(end-8:end) = [0; 0; 0; 0; 0; 0; 0; 0; thf];
ubw(end-8:end) = [0; 0; 0; 0; 0; 0; 0; 0; thf];
else
    lbw(end) = thf;
    ubw(end) = thf;
end
g = [g, {f(Xkn1, [U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k])}];
lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0; 0];
g = [g, {[U1k; U2k; U3k; U4k; U5k; U6k; U7k; U8k]}];
lbg = [lbg; 0; 0; 0; 0; 0; 0; 0; 0];
ubg = [ubg; 0; 0; 0; 0; 0; 0; 0; 0];

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
x1_opt = w_opt(2:17:end);
y1_opt = w_opt(3:17:end);
x2_opt = w_opt(4:17:end);
y2_opt = w_opt(5:17:end);
x3_opt = w_opt(6:17:end);
y3_opt = w_opt(7:17:end);
x4_opt = w_opt(8:17:end);
y4_opt = w_opt(9:17:end);
th_opt = w_opt(10:17:end);
x1d_opt = w_opt(11:17:end);
y1d_opt = w_opt(12:17:end);
x2d_opt = w_opt(13:17:end);
y2d_opt = w_opt(14:17:end);
x3d_opt = w_opt(15:17:end);
y3d_opt = w_opt(16:17:end);
x4d_opt = w_opt(17:17:end);
y4d_opt = w_opt(18:17:end);
subplot(3, 4, 1:4);
plot(x1_opt, y1_opt, 'b-');
hold on
plot(x2_opt, y2_opt, 'g-');
plot(x3_opt, y3_opt, 'c-');
plot(x4_opt, y4_opt, 'm-');
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
legend('$m_{1}$', '$m_{2}$', '$m_{3}$', '$m_{4}$', '$M$', 'Interpreter', 'latex', 'Location',...
    'northeastoutside');
axis equal;
subplot(3, 4, 5:6);
plot(t, th_opt);
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$\theta$, rad', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 4, 7);
plot(t, x1_opt);
hold on
plot(t, x2_opt);
plot(t, x3_opt);
plot(t, x4_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', '$m_{3}$','$m_{4}$', 'Interpreter', 'latex');
subplot(3, 4, 8);
plot(t, y1_opt);
hold on
plot(t, y2_opt);
plot(t, y3_opt);
plot(t, y4_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', '$m_{3}$','$m_{4}$', 'Interpreter', 'latex');
subplot(3, 4, 9)
plot(t(1:end-1), x1d_opt);
hold on
plot(t(1:end-1), y1d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{1}$', '$\dot{y}_{1}$',...
    'Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 4, 10);
plot(t(1:end-1), x2d_opt);
hold on
plot(t(1:end-1), y2d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{2}$', '$\dot{y}_{2}$',...
    'Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 4, 11);
plot(t(1:end-1), x3d_opt);
hold on
plot(t(1:end-1), y3d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{3}$', '$\dot{y}_{3}$',...
    'Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(3, 4, 12);
plot(t(1:end-1), x4d_opt);
hold on
plot(t(1:end-1), y4d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{4}$', '$\dot{y}_{4}$',...
    'Interpreter', 'latex','FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

%% Plots for Report

figure(1);
plot(x1_opt, y1_opt, 'b-');
hold on
plot(x2_opt, y2_opt, 'g-');
plot(x3_opt, y3_opt, 'c-');
plot(x4_opt, y4_opt, 'm-');
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
legend('$m_{1}$', '$m_{2}$', '$m_{3}$', '$m_{4}$', '$M$', 'Interpreter', 'latex', 'Location',...
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
plot(t, x3_opt);
plot(t, x4_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$x$, m', 'Interpreter', 'latex', 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', '$m_{3}$','$m_{4}$', 'Interpreter', 'latex','Location',...
    'northeastoutside');
subplot(2, 1, 2);
plot(t, y1_opt);
hold on
plot(t, y2_opt);
plot(t, y3_opt);
plot(t, y4_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$y$, m', 'Interpreter', 'latex', 'FontSize', 12);
ylim([-1 1]);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
legend('$m_{1}$', '$m_{2}$', '$m_{3}$','$m_{4}$', 'Interpreter', 'latex','Location',...
    'northeastoutside');

figure(4);
subplot(2, 2, 1);
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
subplot(2, 2, 2);
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
subplot(2, 2, 3);
plot(t(1:end-1), x3d_opt);
hold on
plot(t(1:end-1), y3d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{3}$', '$\dot{y}_{3}$',...
    'Interpreter', 'latex','FontSize', 12,'Location',...
    'northeastoutside');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
subplot(2, 2, 4);
plot(t(1:end-1), x4d_opt);
hold on
plot(t(1:end-1), y4d_opt);
hold off
xlabel('$t$, s', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('$u$, m/s', 'Interpreter', 'latex', 'FontSize', 12);
legend('$\dot{x}_{4}$', '$\dot{y}_{4}$',...
    'Interpreter', 'latex','FontSize', 12,'Location',...
    'northeastoutside');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);


%% Model

function [dxdt] = spacecraft_model(x, u)
global m1 m2 m3 m4 M J
% State and control vars.
x1 = x(1);
y1 = x(2);
x2 = x(3);
y2 = x(4);
x3 = x(5);
y3 = x(6);
x4 = x(7);
y4 = x(8);
th = x(9);
x1d = u(1);
y1d = u(2);
x2d = u(3);
y2d = u(4);
x3d = u(5);
y3d = u(6);
x4d = u(7);
y4d = u(8);

% Dynamics
thd = -(M*m1*x1*y1d - M*m1*x1d*y1 + M*m2*x2*y2d - M*m2*x2d*y2 +...
    M*m3*x3*y3d - M*m3*x3d*y3 + M*m4*x4*y4d - M*m4*x4d*y4 +...
    m1*m2*x1*y1d - m1*m2*x1d*y1 - m1*m2*x1*y2d - m1*m2*x2*y1d +...
    m1*m2*x1d*y2 + m1*m2*x2d*y1 + m1*m3*x1*y1d - m1*m3*x1d*y1 +...
    m1*m2*x2*y2d - m1*m2*x2d*y2 + m1*m4*x1*y1d - m1*m4*x1d*y1 -...
    m1*m3*x1*y3d - m1*m3*x3*y1d + m1*m3*x1d*y3 + m1*m3*x3d*y1 +...
    m2*m3*x2*y2d - m2*m3*x2d*y2 + m1*m3*x3*y3d - m1*m3*x3d*y3 -...
    m1*m4*x1*y4d - m1*m4*x4*y1d + m1*m4*x1d*y4 + m1*m4*x4d*y1 -...
    m2*m3*x2*y3d - m2*m3*x3*y2d + m2*m3*x2d*y3 + m2*m3*x3d*y2 +...
    m2*m4*x2*y2d - m2*m4*x2d*y2 + m2*m3*x3*y3d - m2*m3*x3d*y3 -...
    m2*m4*x2*y4d - m2*m4*x4*y2d + m2*m4*x2d*y4 + m2*m4*x4d*y2 +...
    m1*m4*x4*y4d - m1*m4*x4d*y4 + m3*m4*x3*y3d - m3*m4*x3d*y3 +...
    m2*m4*x4*y4d - m2*m4*x4d*y4 - m3*m4*x3*y4d - m3*m4*x4*y3d +...
    m3*m4*x3d*y4 + m3*m4*x4d*y3 + m3*m4*x4*y4d - m3*m4*x4d*y4)/(J*m1 +...
    J*m2 + J*m3 + J*m4 + J*M + m1*m2*x1^2 + m1*m2*x2^2 + m1*m3*x1^2 +...
    m1*m4*x1^2 + m1*m3*x3^2 + m2*m3*x2^2 + m2*m3*x3^2 + m2*m4*x2^2 +...
    m1*m4*x4^2 + m2*m4*x4^2 + m3*m4*x3^2 + m3*m4*x4^2 + m1*m2*y1^2 +...
    m1*m2*y2^2 + m1*m3*y1^2 + m1*m4*y1^2 + m1*m3*y3^2 + m2*m3*y2^2 +...
    m2*m3*y3^2 + m2*m4*y2^2 + m1*m4*y4^2 + m2*m4*y4^2 + m3*m4*y3^2 +...
    m3*m4*y4^2 + M*m1*x1^2 + M*m2*x2^2 + M*m3*x3^2 + M*m4*x4^2 +...
    M*m1*y1^2 + M*m2*y2^2 + M*m3*y3^2 + M*m4*y4^2 - 2*m1*m2*x1*x2 -...
    2*m1*m3*x1*x3 - 2*m1*m4*x1*x4 - 2*m2*m3*x2*x3 - 2*m2*m4*x2*x4 -...
    2*m3*m4*x3*x4 - 2*m1*m2*y1*y2 - 2*m1*m3*y1*y3 - 2*m1*m4*y1*y4 -...
    2*m2*m3*y2*y3 - 2*m2*m4*y2*y4 - 2*m3*m4*y3*y4);
dxdt = [x1d;y1d;x2d;y2d;x3d;y3d;x4d;y4d;thd];

end