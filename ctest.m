clear; clc;

% 1. HELICOPTER PARAMETERS
g     = 9.81;
cla   = 5.7;    volh  = 0.092;   lok   = 2.2819;
mass  = 8005.89;   rho   = 1.225;   vtip  = 200;
diam  = 14.63;  iy    = 69087.10;   mast  = 1.6;
vtip  = 200;    omega = vtip / (diam/2); cds = 1.5;
area  = pi/4 * diam^2;           tau   = 0.1;

% 2. DYNAMIC TRIM CALCULATION 
ct_trim      = (mass * g) / (rho * area * vtip^2);
labi0        = sqrt(ct_trim / 2);
collect_trim = 1.5 * ( (4 * ct_trim) / (volh * cla) + labi0 );

% 3. INTEGRATION SETUP
antal = 4000;
teind = 800; 
stap  = teind / antal;

% 4. INITIAL CONDITIONS
t(1) = 0; u(1) = 0; w(1) = 0; q(1) = 0; pitch(1) = 0; x(1) = 0; z(1) = 0; 
labi(1) = labi0;
% I-ACTION for altitude control and forward velocity control:
u_des = 0; w_des = 0; h_des = 0;
delta_collect(1) = 0;
delta_u(1) = 0;
delta_h(1) = 0;

% Dynamic altitude and position control
h_des_dynamic = 0;
theta_des_dynamic = 0;


% 5. MAIN SIMULATION LOOP
for i = 1:antal
    test17 = true;
    if test17 == true
        u_des = 70;
        u(1) = 0;

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.4;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);
    end



    % --- PHYSICS ENGINE (Unchanged BET/Glauert Model) ---
    qdiml = q(i) / omega;
    vdiml = sqrt(u(i)^2 + w(i)^2) / vtip;
    
    % Flow angles
    if u(i) == 0, phi = sign(w(i))*pi/2; else, phi = atan(w(i)/u(i)); if u(i)<0, phi=phi+pi; end; end
    alfc = longit(i) - phi;
    mu   = vdiml * cos(alfc);
    labc = vdiml * sin(alfc);

    % Flapping a1
    a1(i) = (-16/lok * qdiml + 8/3*mu*collect(i) - 2*mu*(labc + labi(i))) / (1 - 0.5*mu^2);

    % Thrust
    ctelem = cla*volh/4 * ( 2/3*collect(i)*(1 + 1.5*mu^2) - (labc + labi(i)) );
    alfd   = alfc - a1(i);
    ctglau = 2*labi(i) * sqrt( (vdiml*cos(alfd))^2 + (vdiml*sin(alfd) + labi(i))^2 );
    
    thrust = ctelem * rho * vtip^2 * area;
    tilt   = longit(i) - a1(i);

    % Equations of Motion
    udot = -g*sin(pitch(i)) - cds/mass*0.5*rho*u(i)*abs(u(i)) + thrust/mass*sin(tilt) - q(i)*w(i);
    wdot =  g*cos(pitch(i)) - cds/mass*0.5*rho*w(i)*abs(w(i)) - thrust/mass*cos(tilt) + q(i)*u(i);
    qdot = -thrust * mast / iy * sin(tilt);

    delta_collect_dot = (w_des - w(i));
    delta_u_dot = (u_des - u(i));
    delta_h_dot = (h_des - z(i));
    
    % Euler Integration
    u(i+1)     = u(i) + stap * udot;
    w(i+1)     = w(i) + stap * wdot;
    q(i+1)     = q(i) + stap * qdot;
    pitch(i+1) = pitch(i) + stap * q(i);
    x(i+1)     = x(i) + stap * (u(i)*cos(pitch(i)) + w(i)*sin(pitch(i)));
    z(i+1)     = z(i) + stap * (w(i)*cos(pitch(i)) - u(i)*sin(pitch(i)));
    labi(i+1)  = labi(i) + stap * (ctelem - ctglau)/tau;
    t(i+1)     = t(i) + stap;

    delta_collect(i+1) = delta_collect(i) + stap * delta_collect_dot;
    delta_u(i+1) = delta_u(i) + stap * delta_u_dot;
    delta_h(i+1) = delta_h(i) + stap * delta_h_dot;
end

collect(i+1) = collect(i); % Hold last value for plotting
longit(i+1) = longit(i); % Hold last value for plotting

% 6. PLOTTING
close all; 
% White background

figure;
plot(t, pitch*180/pi, 'k-', 'LineWidth', 1.5); 
hold on; grid on;
plot(t, u, 'k--', 'LineWidth', 1.5);
xlabel('Time (s)');
legend("\theta_f (deg)", 'u (m/s)');

figure;
plot(t, u, 'k--', 'LineWidth', 1.5);
hold on; grid on;
plot(t, w, 'k-.', 'LineWidth', 1.5);
xlabel('Time (s)');
legend('u (m/s)', "w (m/s)");

figure;
subplot(2,1,1); plot(t, x, 'k-', 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('x-Position (m)');
subplot(2,1,2); plot(t, z, 'k-', 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Altitude (m)');

figure;
plot(t, pitch*180/pi, 'k-', 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Pitch Angle (deg)');

figure;
subplot(2,1,1); plot(t, longit*180/pi, 'k-', 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Longitudinal Cyclic (deg)');
subplot(2,1,2); plot(t, collect*180/pi, 'k-', 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Collective (deg)');

% Force theme to light
set(gcf,'Color','w');