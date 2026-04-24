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
teind = 900; 
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

transition = false;


% 5. MAIN SIMULATION LOOP
for i = 1:antal

    % --- OPEN LOOP INPUTS ---
    % Default: Stay at trim
    longit(i)  = 0; 
    collect(i) = collect_trim;

    % RANDOM/PULSE INPUT: 
    test1 = false;
    if test1 == true
        if t(i) > 0.5 && t(i) < 1.5
            longit(i) = 1 * pi/180; 
        end
    end

    % PROPORTIONAL FEEDBACK: 
    test2 = false;
    if test2 == true
        % Initial disturbance
        if t(i) > 0.5 && t(i) < 1
            longit(i) = 1 * pi/180; 
        end

        if t(i) > 15 
            longit(i) = 0.2 * pitch(i); 
        end
    end

    % PD FEEDBACK:
    test3 = false;
    if test3 == true
        % Initial disturbance
        if t(i) > 0.5 && t(i) < 1
            longit(i) = 1 * pi/180; 
        end

        if t(i) > 15 
            longit(i) = 0.2 * pitch(i) + 0.2 * q(i); 
        end
    end

    % SAS System:
    test4 = false;
    if test4 == true
        % Initial disturbance
        if t(i) > 0.5 && t(i) < 1
            longit(i) = 2 * pi/180; 
        end

        if t(i) > 15
            longit(i) = -0.1 * q(i); 
        end
    end

    % -5º Attitude:
    test5 = false;
    if test5 == true
        longit(i) = 0.2 * (pitch(i)+5*pi/180) + 0.2 * q(i);
    end

    % Desired vertical velocity, P control action:
    test6 = false;
    if test6 == true
        w_des = -1;
        collect_gen = collect_trim;
        K1 = -0.1;
        collect(i) = collect_gen + K1 * (w_des - w(i));
    end

    % Desired vertical velocity, PI control action:
    test7 = false;
    if test7 == true
        w_des = -1;
        collect_gen = collect_trim;
        K1 = -0.1;
        K2 = -0.01;
        collect(i) = collect_gen + K1 * (w_des - w(i)) + K2 * (delta_collect(i));
    end

    % Altitude Setter: 
    test8 = false;
    if test8 == true
        h_des = 40;
        Kcollect = -0.1;
        Kw = -0.01;
        Kh = 0.05;

        w_des = Kh * (h_des - z(i));
        collect_gen = collect_trim;

        collect(i) = collect_gen + Kw * (w_des - w(i)) + Kcollect * (delta_collect(i));
    end

    % Altitude Setter - Smooth: 
    test9 = false;
    if test9 == true && i > 1
        
        Kw = -0.1;
        Kcollect = -0.01;
        Kh = 0.1;

        max_collect_rate = 2 * pi/180; % Limit climb rate
        max_collect_change = max_collect_rate * stap;

        h_des = 50;
        max_climb_rate = 2/(h_des+0.0000001); 
        if h_des_dynamic < h_des
            h_des_dynamic = h_des_dynamic + (h_des - h_des_dynamic) * stap;
        else
            h_des_dynamic = h_des;
        end

        w_des = Kh * (h_des_dynamic - z(i));
        collect_gen = collect_trim;
        collect_controller = collect_gen + Kw * (w_des - w(i)) + Kcollect * (delta_collect(i));
        collect_change = collect_controller - collect(i-1);
        collect(i) = collect(i-1) + max(min(collect_change, max_collect_change), -max_collect_change);
        collect(i) = max(min(collect(i), 18*pi/180), 2*pi/180);
    end

    % Altitude hold setter to 0m from xm, P control
    test9b = false;
    if test9b == true
        z(1) = 50;
        h_des = 0;
        K_h = -0.001;

        collect(i) = collect_trim + K_h * (h_des - z(i));
    end

    %Altitude hold setter to 0m from xm, PD control
    test9c = false;
    if test9c == true
        z(1) = 50;
        h_des = 0;
        K_h = -0.001;
        K_w = 0.01;

        collect(i) = collect_trim + K_h * (h_des - z(i)) + K_w * w(i);
    end

    %Altitude hold setter to 0m from xm, PD control, filtered
    test9d = false;
    if test9d == true
        z(1) = 50;
        h_des = 0;
        K_h = -0.001;
        K_w = 0.01;

        collect(i) = collect_trim + K_h * (h_des - z(i)) + K_w * w(i);
    end
    

    % Position hold, set position to x_des:
    test10 = false;
    if test10 == true
        u(1) = 0; % m/s,
        x_des = 100;
        K_pos = 0.05;
        K_vel = -0.01;
        u_des = K_pos * (x_des - x(i));
        theta_des = K_vel * (u_des - u(i));

        K_theta = 0.25;
        K_q = 0.1;
        longit(i) = K_theta * (pitch(i)-theta_des) + K_q * q(i); 
    end

    % Position hold with low-pass filter:
    test11 = false;
    if test11 == true
        u(1) = 0; % m/s,
        x_des = 100;
        K_pos = 0.05;
        K_vel = -0.01;
        u_des = K_pos * (x_des - x(i));
        theta_des = K_vel * (u_des - u(i));
        alpha = 0.003;

        theta_des_dynamic = (1-alpha) * theta_des_dynamic + alpha * theta_des;

        K_theta = 0.25;
        K_q = 0.1;
        longit(i) = K_theta * (pitch(i)-theta_des_dynamic) + K_q * q(i); 
    end

    % --- VELOCITY REDUCTION SIMULATIONS ---
    % Trim at 40 kts (20.58 m/s), OPEN LOOP:
    test12 = false;
    if test12 == true
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight
        longit(i) = 1.57 * pi/180;
        collect(i) = 8.08 * pi/180;
    end

    % Test proportional feedback
    test13 = false;
    if test13 == true
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight
        pitch_des = -0.29 * pi/180;
        K_theta = -0.1;
        longit(i) = 1.57 * pi/180 + K_theta * (pitch_des-pitch(i));
    end

    % Test proportional-derivative feedback
    test14 = false;
    if test14 == true
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight
        pitch_des = -0.29 * pi/180;
        K_theta = -0.1;
        K_q = 0.5;
        longit(i) = 1.57 * pi/180 + K_theta * (pitch_des-pitch(i)) + K_q * q(i);
    end

    % PID feedback
    test15 = false;
    if test15 == true
        u_des = 20.58;
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight

        K_vel_error = -0.001;
        pitch_des = K_vel_error * delta_u(i);

        K_theta = -0.5;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);
    end

    % PID feedback with damped outer loop
    test16 = false;
    if test16 == true
        u_des = 20.58;
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.5;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);
    end

    test16a = false;
    if test16a == true
        u_des = 20.58;
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.5;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);

        h_des = 0;

        K_h = -0.0059;
        K_w = 0.01;
        collect(i) = collect_trim + K_h * (h_des - z(i)) + K_w * w(i);
    end

    test16b = false;
    if test16b == true
        u_des = 20.58;
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.5;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);

        h_des = 0;

        K_h = -0.005;
        K_w = 0.01;
        K_deltah = -0.0001;
        collect(i) = collect_trim + K_h * (h_des - z(i)) + K_w * w(i) + K_deltah * delta_h(i);
    end

    % PID feedback with damped outer loop, 0 to 70 m/s and smooth pitch
    test17 = false;
    if test17 == true
        u_des_final = 70;
        u(1) = 0;

        u_des = u_des + 0.01 * (u_des_final - u_des) * stap; % Smoothly increase desired velocity
        
        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.4;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);
    end

    test18 = false;  % From 40 kts to hover
    if test18 == true
        if i == 1
            u_des_final = 0;
            u(1) = 20.58;
            u_des = u(1);
            pitch(1) = -0.283 * pi/180; % 
            labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight
        end

        u_des = u_des + 0.02 * (u_des_final - u_des) * stap; % Smoothly increase desired velocity

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.5;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);
    end
    

    test19 = false; % From 40kts to hover with altitude hold
    if test19 == true
        if i == 1
            h_des = 0;
            pitch_des = 0;
            u_des_final = 0;
            u(1) = 20.58;
            u_des = u(1);
            pitch(1) = -0.283 * pi/180; % 
            labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight
            colect(1) = collect_trim;
        end

        u_des = u_des + 0.02 * (u_des_final - u_des) * stap; % Smoothly increase desired velocity

        K_h = -0.1;
        K_w = 0.05;
        K_deltah = -0.001;

        K_theta = -0.5;
        K_q = 0.5;

        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);
        collect(i) = K_h * (h_des - z(i)) + K_w * w(i) + K_deltah * delta_h(i);
    end

    test20 = true; % Full maneuver from trim at 40kts to standtill, with altitude hold
    if test20 == true && t(i) < 300 % Trim aircraft at 40kts
        u_des = 20.58;
        u(1) = 20.58;
        pitch(1) = -0.29 * pi/180; % 
        labi(1) = 8.54 / (omega * diam/2); % Adjusted for forward flight

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));

        K_theta = -0.5;
        K_q = 0.5;
        longit(i) = K_theta * (pitch_des-pitch(i)) + K_q * q(i);

        h_des = 0;

        K_h = -0.005;
        K_w = 0.01;
        K_deltah = -0.0001;
        collect(i) = collect(i) + K_h * (h_des - z(i)) + K_w * w(i) + K_deltah * delta_h(i);
    end
    if test20 == true && t(i) >= 300 % At t=300s, start reducing velocity to 0

        if transition == false
            fprintf("Transitioning to hover control at t=%.2f seconds\n", t(i));
            h_des = z(i);
            pitch_des_final = 0;
            u_des_final = 0;
            pitch_des = pitch(i);
            u_des = u(i);
            delta_h(i) = 0;
            delta_collect(i) = 0;
            delta_u(i) = 0;
            collect_at_transition = collect(i-1);
            longit_at_transition = longit(i-1);
            transition = true;
        end

        K_vel_error = -0.001;
        K_vel = -0.01;
        pitch_des = K_vel_error * delta_u(i) + K_vel * (u_des - u(i));
        u_des = u_des + 0.02 * (u_des_final - u_des) * stap; % Smoothly increase desired velocity

        K_h = -0.1;
        K_w = 0.05;
        K_deltah = -0.001;

        K_theta = -0.5;
        K_q = 0.5;

        longit(i) = longit_at_transition + K_theta * (pitch_des-pitch(i)) + K_q * q(i);
        collect(i) = collect_at_transition + K_h * (h_des - z(i)) + K_w * w(i) + K_deltah * delta_h(i);
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
set(gcf,'Color','w');  