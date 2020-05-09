
% continuous time plant model
c = 5;
agp_ct = [-1 -3 c; 1 -1 0; 0 0 0];
bgp_ct = [0; 0; 1];
cgp_ct = [1 0 0];
dgp_ct = 0;
gp_ct = ss(agp_ct, bgp_ct, cgp_ct, dgp_ct);


% start with no controller
gc_ct = 1;
sys_ct = feedback(gp_ct, 1);
sys_ct_uncomp = sys_ct;

% plot uncompensated step response
figure
step(sys_ct)
grid
title("Continuous Time Step Response")
ylabel("Pitch Angle (rad)")


% plot multiple step responses with varying sample period
ts_arr = [0.3 0.1 0.03 0.01 0.003];
t_end = 5;

figure
hold on

step(sys_ct, t_end)

for ts = ts_arr
    gp_dt = c2d(gp_ct, ts);
    step(feedback(gp_dt, 1), t_end)
end

grid
title("Continuous and Discrete Time Step Responses")
ylabel("Pitch Angle (rad)")
legend(cat(2, "Continuous", ts_arr+" s"), "Location", "southeast")
hold off


% setup system for controller design
ts = 0.01;
gp_dt = c2d(gp_ct, ts);
sys_dt_uncomp = feedback(gp_dt, 1);


% design PID controller
gp_dt_poles = pole(gp_dt);
k_arr = [0.01 0.005 0.002];
t_end = 2;

figure
hold on

for k = k_arr
    a = [1 ts/2 1/ts -1; -1 ts/2 -2/ts (gp_dt_poles(1)+gp_dt_poles(2)); 0 0 1/ts -(gp_dt_poles(1)*gp_dt_poles(2)); 0 0 0 k];
    b = [0; 0; 0; 1];
    x = linsolve(a, b);

    kp = x(1);
    ki = x(2);
    kd = x(3);

    a0 = kp + ki*ts/2 + kd/ts;
    a1 = -kp + ki*ts/2 - 2*kd/ts;
    a2 = kd/ts;
    gc_dt = tf([a0 a1 a2], [1 -1 0], ts);

    sys_dt = feedback(gp_dt*gc_dt, 1);
    step(sys_dt, t_end)
end

grid
title("PID Step Responses")
ylabel("Pitch Angle (rad)")
legend(k_arr+"", "Location", "southeast")
hold off

sys_dt_pid = sys_dt;


% design pole-placement controller
[agp_dt, bgp_dt, cgp_dt, dgp_dt] = ssdata(gp_dt);
[~, ~, ~, order] = minfo(gp_dt); % = 3
controllability = rank(ctrb(agp_dt, bgp_dt)); % = 3, full rank

r = 0.98315;
om = pi/1000;
sys_poles = [r r*exp(1i*om) r*exp(-1i*om)];
k = place(agp_dt, bgp_dt, sys_poles);
sys_dt = ss(agp_dt - bgp_dt*k, bgp_dt, cgp_dt, dgp_dt, 0.001); % change ts for this controller

figure
hold on
step(sys_dt, t_end);
step(sys_dt_pid, t_end);
grid
title("Pole Placement versus PID Step Response")
ylabel("Pitch Angle (rad)")
legend(["Pole Placement" "PID"], "Location", "southeast")
hold off

sys_dt_sf = sys_dt;


% design one-step deadbeat controller
[z, p, k] = zpkdata(gp_dt);
gc_dt = zpk(p, z, 1/k, ts) * tf(1, [1 -1], ts);
sys_dt = feedback(gp_dt*gc_dt, 1);

[y, t, x] = step(ss(sys_dt), t_end);
figure
plot(t, y)
grid
title("One-Step Deadbeat Step Response")
xlabel("Time (seconds)")
ylabel("Pitch Angle (rad)")

figure
plot(t, x(:,3))
grid
title("One-Step Deadbeat Step Response Control Signal")
xlabel("Time (seconds)")
ylabel("Rate Command (rad/s)")


% design one-step deadbeat controller following Dahlin's method
lambda_arr = [0.02 0.04 0.08];

figure
hold on

for lambda = lambda_arr
    [z, p, k] = zpkdata(gp_dt);
    gc_dt = zpk(p, z, 1/k, ts) * tf(1-exp(-ts/lambda), [1 -1], ts);
    sys_dt = feedback(gp_dt*gc_dt, 1);
    
    step(ss(sys_dt), t_end);
end

grid
title("Dahlin Algorithm One-Step Deadbeat Step Response")
ylabel("Pitch Angle (rad)")
legend(lambda_arr+"", "Location", "southeast")
hold off

figure
hold on

for lambda = lambda_arr
    [z, p, k] = zpkdata(gp_dt);
    gc_dt = zpk(p, z, 1/k, ts) * tf(1-exp(-ts/lambda), [1 -1], ts);
    sys_dt = feedback(gp_dt*gc_dt, 1);
    
    [~, t, x] = step(ss(sys_dt), t_end);
    plot(t, x(:,3))
end

grid
title("Dahlin Algorithm One-Step Deadbeat Step Response Control Signal")
xlabel("Time (seconds)")
ylabel("Rate Command (rad/s)")
legend(lambda_arr+"", "Location", "southeast")
hold off

lambda = 0.04;
[z, p, k] = zpkdata(gp_dt);
gc_dt = zpk(p, z, 1/k, ts) * tf(1-exp(-ts/lambda), [1 -1], ts);
sys_dt_db = feedback(gp_dt*gc_dt, 1);


% plot all step responses
t_end = 8;
figure
hold on
step(sys_ct_uncomp, t_end);
step(sys_dt_uncomp, t_end);
step(sys_dt_pid, t_end);
step(sys_dt_sf, t_end);
step(sys_dt_db, t_end);

grid
title("Step Responses of All Systems")
xlabel("Time (seconds)")
ylabel("Pitch Angle (rad)")
legend(["Continuous-Time, Uncompensated" "Uncompensated" "PID" "Pole Placement" "Dahlin Deadbeat"], "Location", "southeast")
hold off


% plot all controlled step responses
t_end = 1;
figure
hold on
step(sys_dt_pid, t_end);
step(sys_dt_sf, t_end);
step(sys_dt_db, t_end);

grid
title("Step Responses of All Controlled Systems")
xlabel("Time (seconds)")
ylabel("Pitch Angle (rad)")
legend(["PID" "Pole Placement" "Dahlin Deadbeat"], "Location", "southeast")
hold off


% plot all controlled step responses
t_end = 1;
figure
hold on
[~, t, x] = step(sys_dt_pid, t_end);
plot(t, x(:,3))
[~, t, x] = step(sys_dt_sf, t_end);
plot(t, x(:,3))
[~, t, x] = step(sys_dt_db, t_end);
plot(t, x(:,3))

grid
title("Step Response Control Signals of All Controlled Systems")
xlabel("Time (seconds)")
ylabel("Pitch Angle (rad)")
legend(["PID" "Pole Placement" "Dahlin Deadbeat"], "Location", "southeast")
hold off
