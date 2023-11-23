clear;clc 
clearvars;
% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 20;                % [grid points]
pml_y_size = 20;                % [grid points]
Nx = 540 - 2 * pml_x_size;      % [grid points] 340
Ny = 540 - 2 * pml_y_size;      % [grid points] 340
x = 240e-3;                      % [m]240
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1540;                      % [m/s]1480
medium.density = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;
  
% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

% % create initial pressure distribution using makeDisc
% disc_x_pos = Nx/2;    % [grid points]
% disc_y_pos = Ny/2;    % [grid points]
% disc_radius = 15e-3/dx;    % [grid points]
% disc_1 =makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);
% 
% medium.sound_speed = medium.sound_speed * ones(Nx, Ny);   % [m/s]
% medium.sound_speed(disc_1==1) = 3000;       % [m/s]
% medium.density = 1000 * ones(Nx, Ny);       % [kg/m^3]
% medium.density(disc_1==1) = 3000;          % [kg/m^3]

% % define properties of the input signal
% source_strength = 1e6;          % [Pa]
% tone_burst_freq = 3e6;        % [Hz]
% tone_burst_cycles = 5;

% create the time array
t_end = 307e-6;   % [s]
cfl=0.3;
kgrid.makeTime(medium.sound_speed, cfl, t_end);
% % create the input signal using toneBurst 
% source.p =20*toneBurst(100e6, tone_burst_freq, tone_burst_cycles);
% define a centered Cartesian circular source
source_radius = 100e-3;     % [m]
source_angle = 2*pi;      % [rad]
source_pos = [0, 0];        % [m]
num_source_points = 256 ;
cart_source_mask = makeCartCircle(source_radius, num_source_points, source_pos, source_angle);
source.p_mask = cart2grid(kgrid, cart_source_mask);
% define a centered Cartesian circular sensor
sensor_radius = 100e-3;     % [m]
sensor_angle = 2*pi;      % [rad]
sensor_pos = [0, 0];        % [m]
num_sensor_points =256 ;
cart_sensor_mask = makeCartCircle(sensor_radius, num_sensor_points, sensor_pos, sensor_angle);
% assign to sensor structure
sensor.mask = cart_sensor_mask;

% set the input options
input_args = {'PMLInside', false, 'PlotPML', false, 'DataCast', 'gpuArray-double'};

% % run the simulation
% sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});


% add noise to the recorded sensor data
% signal_to_noise_ratio = 40;	% [dB]
% sensor_data = addNoise(sensor_data, signal_to_noise_ratio, 'peak');

% create a second computation grid for the reconstruction to avoid the
% inverse crime
% set the size of the perfectly matched layer (PML)
pml_x_size = 20;                % [grid points]
pml_y_size = 20;                % [grid points]
Nx = 540 - 2 * pml_x_size;      % [grid points]
Ny = 540 - 2 * pml_y_size;      % [grid points]
x = 240e-3;                      % [m]
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
kgrid_recon = kWaveGrid(Nx, dx, Ny, dy);%%%%wy%%%%
% use the same time array for the reconstruction
kgrid_recon.setTime(kgrid.Nt, kgrid.dt);
% sensor_data=gather(data_need_filter);
% % reset the initial pressure
% source.p0 = 0;
% assign the time reversal data
% load E:\user-MATLAB\Test\k-wave\自发自收\32发32收\sum_a_data_need_filter
load E:\user-MATLAB\Test\实测\256合成孔径间隔1\data2\breast_model_white2\image_data\sum_new_data_SA1
sum_new_data=circshift(sum_new_data,[1 -15]);
sensor.time_reversal_boundary_data = sum_new_data;
%% 高斯滤波filter the sensor data using a Gaussian filter
% Fs = 1/kgrid.dt;        % [Hz]
% center_freq = 3e6;      % [Hz]
% bandwidth = 75;        % [%]
% sum_a_data_need_filter_move_gaussian = gaussianFilter(sum_a_data_need_filter_move, Fs, center_freq, bandwidth);
% sensor.time_reversal_boundary_data = sum_a_data_need_filter_move_gaussian;
%% 高通滤波filter the sensor data using a high pass filter
% Fs = 1/kgrid.dt;        % [Hz]
% cutoff_freq = 3e6;      % [Hz]
% sum_a_data_need_filter_move_highpass = zeros(size(sum_a_data_need_filter_move));
% for index = 1:sum(sensor.mask(:))
%     sum_a_data_need_filter_move_highpass(index, :) = applyFilter(sum_a_data_need_filter_move(index, :), Fs, cutoff_freq, 'LowmPass', 'ZeroPhase', true);
% end
% sensor.time_reversal_boundary_data = sum_a_data_need_filter_move_highpass;
%%
% run the time-reversal reconstruction
% sensor.time_reversal_boundary_data = fliplr(sensor.time_reversal_boundary_data);
p0_recon = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor, input_args{:});
p0_recon1=gather(p0_recon);
%  Do logarithmic compression to a 60 dB dynamic range将对数压缩到60db动态范围

% log_p0_recon=20.*log10(p0_recon1);
% log_p0_recon=log_p0_recon-max(max(p0_recon)) + 50;
% log_p0_recon=127*log_p0_recon/60;
% % assign to sensor structure
% sensor_radius_grid_points = round(sensor_radius / kgrid_recon.dx);
% binary_sensor_mask = makeCircle(kgrid_recon.Nx, kgrid_recon.Ny, kgrid_recon.Nx/2 + 1, kgrid_recon.Ny/2 + 1, sensor_radius_grid_points, sensor_angle);
% 
% % assign to sensor structure
% sensor.mask = binary_sensor_mask;

% % interpolate data to remove the gaps and assign to sensor structure
% sensor.time_reversal_boundary_data = interpCartData(kgrid_recon, sensor_data, cart_sensor_mask, binary_sensor_mask);

% % run the time-reversal reconstruction 
% p0_recon_interp = kspaceFirstOrder2D(kgrid_recon, medium, source, sensor, input_args{:});

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the initial pressure and sensor distribution    p0是传感器，cart2grid原始图
figure;
% imagesc(kgrid.y_vec* 1e3, kgrid.x_vec * 1e3, 1000*cart2grid(kgrid, cart_sensor_mask) , [0, 1]);
image(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, cart2grid(kgrid, cart_sensor_mask));
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% % plot the simulated sensor data
% figure;
% imagesc(sensor_data, [-1, 1]);
% colormap(getColorMap);
% ylabel('Sensor Position');
% xlabel('Time Step'); 
% colorbar;
% load E:\user-MATLAB\Test\k-wave\自发自收\网格340_340\合成孔径间隔1阵元\4发4收\p0_recon1
% plot the reconstructed initial pressure 
% figure;
hold on
% imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon+cart2grid(kgrid, cart_sensor_mask), [0, 5]);
imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon1+cart2grid(kgrid, cart_sensor_mask));%);
% colormap(gray);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
colorbar
colormap(jet)
% figure;
% imshow ( p0_recon);
% 
% axis([-100 100 -100 100])
% ylabel('x-position [mm]');
% xlabel('y-position [mm]');
% axis image;

% % plot the reconstructed initial pressure using the interpolated重建的图
% figure;
% % imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon_interp, [-1, 1]);
% imagesc(kgrid_recon.y_vec * 1e3, kgrid_recon.x_vec * 1e3, p0_recon_interp, [-1, 1]);
% colormap(getColorMap);
% ylabel('x-position [mm]');
% xlabel('y-position [mm]');
% axis image;
% 
% % plot a profile for comparison
% slice_pos = 4.5e-3;  % [m] location of the slice from top of grid [m]
% figure;
% plot(kgrid.y_vec * 1e3, p0(round(slice_pos/kgrid.dx), :), 'k--', ...
%      kgrid_recon.y_vec * 1e3, p0_recon(round(slice_pos/kgrid_recon.dx), :), 'r-', ...
%      kgrid_recon.y_vec * 1e3, p0_recon_interp(round(slice_pos/kgrid_recon.dx), :), 'b-');
% xlabel('y-position [mm]');
% ylabel('Pressure');
% legend('Initial Pressure', 'Point Reconstruction', 'Interpolated Reconstruction');
% axis tight;
% set(gca, 'YLim', [0 2.1]);