%==========================================================================
% AE4314 - HELICOPTER FLIGHT MECHANICS
% ADS-33 Sustained Hover - Phase 2: PID Controller
%
% DESCRIPTION:
%   Simulates a 5-minute (300 s) stabilised hover at a fixed target point.
%   Two simultaneous feedback loops maintain position within +/-3.05 m:
%
%   1. LONGITUDINAL (Cyclic pitch, theta_c) - Cascaded PD controller
%      Outer loop : desired fuselage pitch from position + velocity error
%      Inner loop : cyclic pitch command from attitude error + pitch rate
%
%   2. VERTICAL (Collective pitch, theta_0) - PID controller
%      P  : proportional on climb-rate error
%      I  : integral state (deltatheta) eliminates steady-state height error
%      Outer: height error -> desired climb rate
%
% AIRCRAFT DATA: same helicopter as sim3_cyclic.m (mass = 2200 kg)
%
% INTEGRATION: Euler forward, dt = 0.1 s, 3000 steps (300 s total)
%==========================================================================

clear; clc;

%--------------------------------------------------------------------------
% 1.  HELICOPTER PARAMETERS  (identical to instructor scripts, have to adapt to Apache)
%--------------------------------------------------------------------------
g     = 9.81;
cla   = 5.7;          % NACA 0012 lift-curve slope
volh  = 0.075;        % blade solidity
lok   = 6;            % Lock number
cds   = 1.5;          % equivalent flat-plate drag area (m^2)
mass  = 2200;         % kg
rho   = 1.225;        % kg/m^3
vtip  = 200;          % m/s  blade tip speed
diam  = 2 * 7.32;     % m    rotor diameter
iy    = 10615;        % kg m^2  pitch moment of inertia
mast  = 1;            % m    hub-to-CG distance (positive nose-up moment arm)
omega = vtip / (diam/2);
area  = pi/4 * diam^2;
tau   = 0.1;          % dynamic inflow time constant (s)

%--------------------------------------------------------------------------
% 2.  TRIM COLLECTIVE  (hover equilibrium)
%--------------------------------------------------------------------------
% Trim inflow from momentum theory: labi0 = sqrt(T/(2*rho*A)) / vtip
labi0        = sqrt(mass*g / (area * 2 * rho)) / vtip;
collect_trim = 6 * pi/180;   % 6 deg as in lecture / previous scripts

%--------------------------------------------------------------------------
% 3.  MISSION TARGETS
%--------------------------------------------------------------------------
xtarget = 0;      % desired horizontal position (m)
htarget = 0;      % desired altitude, defined as h = -z = 0 at start

%--------------------------------------------------------------------------
% 4.  CONTROLLER GAINS
%
%   Horizontal (cascaded PD on cyclic):
%     Outer:  pitch_wens = -Kx*(xtarget - x)  +  Ku*u
%     Inner:  longit     =  Kp*(pitch - pitch_wens)  +  Kd*q
%
%   Vertical (PID on collective):
%     cdes           = K3*(htarget - h)
%     deltatheta_dot = cdes - c
%     collect        = thetagen + K1*(cdes - c) + K2*deltatheta
%
%   Tuning approach:
%     - Start with inner loop (attitude hold) before outer loop (position)
%     - K3 sets how aggressively altitude error drives a climb rate demand
%     - K1 provides fast damping on climb-rate error; K2 kills static error
%--------------------------------------------------------------------------

% --- Cyclic / horizontal ---
Kx = 0.010;    % position error  -> desired pitch (rad per m)
Ku = 0.060;    % velocity feed-forward (rad per m/s)
Kp = 0.800;    % attitude error  -> cyclic (rad per rad)
Kd = 0.400;    % pitch rate damp -> cyclic (rad per rad/s)

% --- Collective / vertical ---
K1      = 0.012;           % climb-rate error -> collective (rad per m/s)
K2      = 0.005;           % integral state   -> collective (rad per m)
K3      = 0.15;            % height error     -> desired climb rate (1/s)
thetagen = collect_trim;   % feed-forward baseline collective

%--------------------------------------------------------------------------
% 5.  INTEGRATION SETUP
%--------------------------------------------------------------------------
antal  = 3000;
teind  = 300;          % 5 minutes = 300 seconds
t0     = 0;
stap   = (teind - t0) / antal;   % 0.1 s

%--------------------------------------------------------------------------
% 6.  INITIAL CONDITIONS
%--------------------------------------------------------------------------
t(1)          = t0;
u(1)          = 1;
w(1)          = 10;
q(1)          = 0;
pitch(1)      = 0;
x(1)          = 0;
z(1)          = 0;
labi(1)       = labi0;
deltatheta(1) = 0;    % integral state for altitude PID

%--------------------------------------------------------------------------
% 7.  MAIN SIMULATION LOOP
%--------------------------------------------------------------------------
for i = 1:antal

    %----------------------------------------------------------------------
    % A.  CYCLIC PITCH COMMAND  (cascaded PD – horizontal controller)
    %----------------------------------------------------------------------
    pitch_wens = -Kx * (xtarget - x(i)) + Ku * u(i);
    longit(i)  =  Kp * (pitch(i) - pitch_wens) + Kd * q(i);

    %----------------------------------------------------------------------
    % B.  CLIMB RATE AND HEIGHT
    %----------------------------------------------------------------------
    c(i) = u(i)*sin(pitch(i)) - w(i)*cos(pitch(i));
    h(i) = -z(i);

    %----------------------------------------------------------------------
    % C.  COLLECTIVE PITCH COMMAND  (PID – vertical controller)
    %----------------------------------------------------------------------
    cdes              = K3 * (htarget - h(i));
    deltatheta_dot(i) = cdes - c(i);
    collect(i)        = thetagen + K1*(cdes - c(i)) + K2*deltatheta(i);

    % Clamp collective to physically reasonable limits (3 to 12 deg)
    collect(i) = max(3*pi/180, min(12*pi/180, collect(i)));

    %----------------------------------------------------------------------
    % D.  NON-DIMENSIONAL VELOCITY AND FLOW ANGLE
    %----------------------------------------------------------------------
    qdiml(i)  = q(i) / omega;
    vdiml(i)  = sqrt(u(i)^2 + w(i)^2) / vtip;

    % Quadrant-safe flight-path angle (handles reverse flight u < 0)
    if u(i) == 0
        if w(i) > 0
            phi(i) =  pi/2;
        else
            phi(i) = -pi/2;
        end
    else
        phi(i) = atan(w(i) / u(i));
        if u(i) < 0
            phi(i) = phi(i) + pi;
        end
    end

    alfc(i) = longit(i) - phi(i);
    mu(i)   = vdiml(i) * cos(alfc(i));
    labc(i) = vdiml(i) * sin(alfc(i));

    %----------------------------------------------------------------------
    % E.  FLAPPING ANGLE  a1  (blade element theory)
    %----------------------------------------------------------------------
    teller = -16/lok * qdiml(i) ...
             + 8/3   * mu(i) * collect(i) ...
             - 2     * mu(i) * (labc(i) + labi(i));
    a1(i) = teller / (1 - 0.5*mu(i)^2);

    %----------------------------------------------------------------------
    % F.  THRUST COEFFICIENT  (BEM element and Glauert)
    %----------------------------------------------------------------------
    ctelem(i) = cla*volh/4 * ( 2/3*collect(i)*(1 + 1.5*mu(i)^2) ...
                                - (labc(i) + labi(i)) );

    alfd(i)   = alfc(i) - a1(i);
    ctglau(i) = 2*labi(i) * sqrt( (vdiml(i)*cos(alfd(i)))^2 ...
                                 + (vdiml(i)*sin(alfd(i)) + labi(i))^2 );

    %----------------------------------------------------------------------
    % G.  THRUST AND ROTOR TILT
    %----------------------------------------------------------------------
    thrust(i)  = ctelem(i) * rho * vtip^2 * area;
    helling(i) = longit(i) - a1(i);
    vv(i)      = vdiml(i) * vtip;

    %----------------------------------------------------------------------
    % H.  EQUATIONS OF MOTION
    %----------------------------------------------------------------------
    udot(i) = -g*sin(pitch(i)) ...
              - cds/mass * 0.5*rho * u(i)*vv(i) ...
              + thrust(i)/mass * sin(helling(i)) ...
              - q(i)*w(i);

    wdot(i) =  g*cos(pitch(i)) ...
              - cds/mass * 0.5*rho * w(i)*vv(i) ...
              - thrust(i)/mass * cos(helling(i)) ...
              + q(i)*u(i);

    qdot(i) = -thrust(i) * mast / iy * sin(helling(i));

    pitchdot(i) = q(i);

    xdot(i) = u(i)*cos(pitch(i)) + w(i)*sin(pitch(i));

    zdot(i) = -c(i);

    % Quasi-dynamic inflow (Peters-He momentum theory)
    labidot(i) = (ctelem(i) - ctglau(i)) / tau;

    %----------------------------------------------------------------------
    % I.  EULER INTEGRATION  – advance all states one step
    %----------------------------------------------------------------------
    u(i+1)           = u(i)           + stap * udot(i);
    w(i+1)           = w(i)           + stap * wdot(i);
    q(i+1)           = q(i)           + stap * qdot(i);
    pitch(i+1)       = pitch(i)       + stap * pitchdot(i);
    x(i+1)           = x(i)           + stap * xdot(i);
    z(i+1)           = z(i)           + stap * zdot(i);
    labi(i+1)        = labi(i)        + stap * labidot(i);
    deltatheta(i+1)  = deltatheta(i)  + stap * deltatheta_dot(i);

    t(i+1) = t(i) + stap;

end % for i

%--------------------------------------------------------------------------
% 8.  ADS-33 COMPLIANCE CHECK
%--------------------------------------------------------------------------
x_max = max(abs(x));
h_vec = -z;
h_max = max(abs(h_vec));

fprintf('\n=== ADS-33 Sustained Hover – Compliance Report ===\n');
fprintf('Max horizontal drift  |x|_max = %6.3f m  (limit 3.05 m) --> %s\n', ...
        x_max, pass_fail(x_max, 3.05));
fprintf('Max altitude variation|h|_max = %6.3f m  (limit 3.05 m) --> %s\n', ...
        h_max, pass_fail(h_max, 3.05));

%--------------------------------------------------------------------------
% 9.  PLOTS
%--------------------------------------------------------------------------
t_plot = t(1:antal);

figure(1);
subplot(2,1,1);
plot(t_plot, x(1:antal), 'b', 'LineWidth', 1.5); hold on;
yline( 3.05,'r--','LineWidth',1.2);
yline(-3.05,'r--','LineWidth',1.2);
xlabel('t (s)'); ylabel('x (m)');
title('Horizontal Position  –  ADS-33 limit \pm3.05 m');
legend('x(t)','ADS-33 limit'); grid on;

subplot(2,1,2);
plot(t_plot, -z(1:antal), 'b', 'LineWidth', 1.5); hold on;
yline( 3.05,'r--','LineWidth',1.2);
yline(-3.05,'r--','LineWidth',1.2);
xlabel('t (s)'); ylabel('h (m)');
title('Altitude  –  ADS-33 limit \pm3.05 m');
legend('h(t)','ADS-33 limit'); grid on;

figure(2);
subplot(2,1,1);
plot(t_plot, u(1:antal), 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('u (m/s)'); title('Longitudinal velocity u'); grid on;

subplot(2,1,2);
plot(t_plot, w(1:antal), 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('w (m/s)'); title('Vertical velocity w'); grid on;

figure(3);
subplot(2,1,1);
plot(t_plot, pitch(1:antal)*180/pi, 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('\theta_f (deg)'); title('Fuselage pitch attitude'); grid on;

subplot(2,1,2);
plot(t_plot, q(1:antal)*180/pi, 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('q (deg/s)'); title('Pitch rate q'); grid on;

figure(4);
subplot(2,1,1);
plot(t_plot, longit(1:antal)*180/pi, 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('\theta_c (deg)'); title('Cyclic pitch command \theta_c'); grid on;

subplot(2,1,2);
plot(t_plot, collect(1:antal)*180/pi, 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('\theta_0 (deg)'); title('Collective pitch command \theta_0'); grid on;

figure(5);
plot(t_plot, labi(1:antal), 'LineWidth', 1.5);
xlabel('t (s)'); ylabel('\lambda_i (-)'); title('Non-dimensional inflow \lambda_i'); grid on;

%--------------------------------------------------------------------------
% HELPER FUNCTION
%--------------------------------------------------------------------------
function s = pass_fail(val, lim)
    if val <= lim
        s = 'PASS';
    else
        s = 'FAIL';
    end
end