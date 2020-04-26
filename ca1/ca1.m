
%% setup

fc_z = 100;
fr_z = 150;
omc_z = 2*pi * fc_z; % = 628 rad/s
omr_z = 2*pi * fr_z; % = 942 rad/s

a_sb_db = 10; % minimum attenuation in stopband
a_sb = 1/power(10, -a_sb_db/20);

fs = 1000;
t = 1/fs;

% prewarp if necessary
om_warp = 0.3/t; % = 300 rad/s, both omc and omr need prewarping
omc_s = 2/t * tan(omc_z * t / 2);
omr_s = 2/t * tan(omr_z * t / 2);


%% part 1, butterworth, indirect approach

order = ceil(log10(a_sb^2 - 1) / (2 * log10(omr_s/omc_s))); % = 3
[z, p, k] = buttap(order);
[num_s, den_s] = zp2tf(z, p, k);
[num_s, den_s] = lp2lp(num_s, den_s, omc_s);

[num_z, den_z] = bilinear(num_s, den_s, fs);
figure
freqz(num_z, den_z, 256, fs);


%% part 1, butterworth, direct approach

[num, den] = butter(order, fc_z / (fs/2));
figure
freqz(num, den, 256, fs)


%% part 2, chebyshev

a_pb_db = 1; % maximum attenuation in passband
ep = sqrt(power(10, a_pb_db / 10) - 1);
order = ceil(acosh(power(10, a_sb/20 - log10(ep^2))) / acosh(omr_s / omc_s));

[num, den] = cheby1(3, a_pb_db, fc_z / (fs/2));
figure
freqz(num, den, 256, fs)


%% part 3, bessel

%d0 = power(omc_s, order); % same order as part 2 (3)
num_s = 15;
den_s = [1 6 15 15]; % 3rd order bessel
[num_s, den_s] = lp2lp(num_s, den_s, omc_s);

[num_z, den_z] = bilinear(num_s, den_s, fs);
figure
freqz(num_z, den_z, 256, fs)
