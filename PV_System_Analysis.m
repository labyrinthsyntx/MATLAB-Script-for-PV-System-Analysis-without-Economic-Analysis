
% MATLAB Code for PV System Analysis without Economic Analysis
% Author: Derek / Team
% Location: Sacramento International Airport (38.7°N, -121.6°E)
% Purpose: Detailed irradiance calculations, PV array performance modeling,
%          environmental impact assessment, and additional performance metrics

%% Clear workspace and command window
clear all; close all; clc;

%% 1. Parameters Initialization

% Constants
G_sc = 1367;                 % Solar constant in W/m²
lat = 38.7;                  % Latitude in degrees (positive for Northern Hemisphere)
lon = -121.6;                % Longitude in degrees (negative for West)
tilt = 30;                   % Tilt angle of the PV panel in degrees
azimuth = 180;               % Surface azimuth angle in degrees (180° for south-facing in Northern Hemisphere)
albedo = 0.2;                % Ground reflectance (dimensionless)
day = 315;                   % Day of the year (November 10, 2024 is the 315th day)

% Time variables
timezone = -8;               % Time zone offset from UTC (Pacific Standard Time)
hours = 0:23;                % Hours of the day

% PV Panel Specifications
P_rated = 415;               % Rated power of each panel in watts (W)
gamma = -0.0045;             % Temperature coefficient (% per °C)
T_STC = 25;                  % Standard Test Condition temperature in degrees Celsius (°C)
NOCT = 45;                   % Nominal Operating Cell Temperature in degrees Celsius (°C)
panel_area = 1.94;           % Area of each panel in square meters (m²)
efficiency = 18;             % PV panel efficiency in percent (%), from the RETScreen report

% Environmental Conditions
ambient_temp = 16.4;         % Average annual temperature in degrees Celsius (°C) from RETScreen report

% System Configuration
installed_capacity = 721.61; % Installed capacity in kilowatts (kW) from RETScreen report
P_rated_total_W = installed_capacity * 1e3; % Total installed capacity in watts (W)
num_panels = ceil(P_rated_total_W / P_rated); % Total number of panels (rounded up)
total_area = num_panels * panel_area;        % Total PV array area in m²

% Correction Factors
CF_shading = 0.95;           % Shading loss factor (5% loss)
CF_soiling = 0.90;           % Soiling loss factor (10% loss)
eta_inverter = 0.96;         % Inverter efficiency (96%)
eta_wiring = 0.98;           % Wiring and mismatch efficiency (98%)

% Degradation Parameters
degradation_rate = 0.005;    % Annual degradation rate (0.5%)
years = 25;                  % Number of years for analysis

% Solar Position Parameters
n = day;                     % Day of the year
phi = lat;                   % Latitude (degrees North)
longitude = lon;             % Longitude (degrees East, negative for West)

% Actual Solar Radiation Data
H_actual = 4.68;             % Actual solar radiation (kWh/m²/day) from RETScreen report

% Environmental Parameters
emission_factor_grid = 0.5;  % kg CO₂ per kWh (average value)

% Preallocate arrays
G_oh = zeros(1,24); % Extraterrestrial horizontal irradiance
GHI = zeros(1,24);  % Global horizontal irradiance
DHI = zeros(1,24);  % Diffuse horizontal irradiance
DNI = zeros(1,24);  % Direct normal irradiance
theta_z = zeros(1,24); % Solar zenith angle
theta = zeros(1,24);   % Angle of incidence
G_bT = zeros(1,24);    % Beam component on tilted surface
G_dT = zeros(1,24);    % Diffuse component on tilted surface
G_rT = zeros(1,24);    % Ground-reflected component on tilted surface
G_T = zeros(1,24);     % Total irradiance on tilted surface

%% 2. Solar Irradiance Calculations

% Calculate declination angle (delta)
delta = 23.45 * sind(360*(284+day)/365);

% Equation of time (EoT) in minutes
B = 360*(day - 81)/364;
EoT = 9.87 * sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);

% Loop over each hour
for i = 1:length(hours)
    % Local standard time
    t = hours(i);

    % Solar time correction
    L_st = timezone * 15; % Reference longitude for the time zone
    TC = 4 * (lon - L_st) + EoT; % Time correction in minutes
    t_solar = t + TC/60; % Solar time in hours

    % Hour angle (H)
    H = 15 * (t_solar - 12); % Hour angle in degrees

    % Solar zenith angle (theta_z)
    theta_z(i) = acosd(sind(lat)*sind(delta) + cosd(lat)*cosd(delta)*cosd(H));

    % Solar azimuth angle (gamma_s)
    gamma_s = atand(sind(H) / (sind(lat)*cosd(H) - cosd(lat)*tand(delta)));
    if H > 0 % Afternoon
        gamma_s = gamma_s + 180;
    else % Morning
        gamma_s = gamma_s + 360;
    end

    % Angle of incidence (theta)
    theta(i) = acosd( sind(delta)*sind(lat)*cosd(tilt) - sind(delta)*cosd(lat)*sind(tilt)*cosd(azimuth) ...
                + cosd(delta)*cosd(lat)*cosd(tilt)*cosd(H) + cosd(delta)*sind(lat)*sind(tilt)*cosd(azimuth)*cosd(H) ...
                + cosd(delta)*sind(tilt)*sind(azimuth)*sin(H) );

    % Extraterrestrial normal irradiance (G_on)
    G_on = G_sc * (1 + 0.033 * cosd(360*day/365));

    % Extraterrestrial horizontal irradiance (G_oh)
    G_oh(i) = G_on * cosd(theta_z(i));
    if G_oh(i) < 0
        G_oh(i) = 0;
    end

    % Clearness index (K_t)
    % Use actual daily solar radiation data
    GHI_daily = H_actual * 1000 / 24; % Convert kWh/m²/day to W/m²/hour
    if G_oh(i) > 0
        K_t = GHI_daily / G_oh(i); % Calculate clearness index
        K_t = max(min(K_t, 0.82), 0); % Constrain K_t between 0 and 0.82
    else
        K_t = 0;
    end

    if G_oh(i) > 0
        % Global horizontal irradiance (GHI)
        GHI(i) = K_t * G_oh(i);

        % Estimate DHI using Erbs model
        if K_t <= 0.22
            K_d = 1 - 0.09*K_t;
        elseif K_t <= 0.80
            K_d = 0.9511 - 0.1604*K_t + 4.388*K_t^2 - 16.638*K_t^3 + 12.336*K_t^4;
        else
            K_d = 0.165;
        end
        DHI(i) = K_d * GHI(i);

        % Direct normal irradiance (DNI)
        DNI(i) = (GHI(i) - DHI(i)) / cosd(theta_z(i));

        % Beam component on tilted surface (G_bT)
        G_bT(i) = DNI(i) * cosd(theta(i));
        if G_bT(i) < 0
            G_bT(i) = 0;
        end

        % Diffuse component on tilted surface (G_dT)
        % Using the Hay-Davies model for anisotropic sky
        A = max(DNI(i) / G_on, 0); % Anisotropy index
        F = max(0, cosd(theta(i)) / cosd(theta_z(i))); % Circumsolar factor
        G_dT(i) = DHI(i) * (A * F + (1 - A) * (1 + cosd(tilt)) / 2);

        % Ground-reflected component (G_rT)
        G_rT(i) = GHI(i) * albedo * (1 - cosd(tilt)) / 2;

        % Total irradiance on tilted surface (G_T)
        G_T(i) = G_bT(i) + G_dT(i) + G_rT(i);
    else
        % Night time
        GHI(i) = 0;
        DHI(i) = 0;
        DNI(i) = 0;
        G_bT(i) = 0;
        G_dT(i) = 0;
        G_rT(i) = 0;
        G_T(i) = 0;
    end
end

% Plot the results
figure;
plot(hours, G_T, 'r', 'LineWidth', 2);
hold on;
plot(hours, G_bT, 'b--', 'LineWidth', 1.5);
plot(hours, G_dT, 'g--', 'LineWidth', 1.5);
plot(hours, G_rT, 'k--', 'LineWidth', 1.5);
xlabel('Hour of the Day');
ylabel('Irradiance (W/m^2)');
title('Hourly Irradiance Components on Tilted Surface');
legend('Total Irradiance (G_T)', 'Beam Component (G_bT)', 'Diffuse Component (G_dT)', 'Ground-reflected Component (G_rT)');
grid on;

% Plot solar zenith and incidence angles
figure;
plot(hours, theta_z, 'm-', 'LineWidth', 2);
hold on;
plot(hours, theta, 'c--', 'LineWidth', 2);
xlabel('Hour of the Day');
ylabel('Angle (Degrees)');
title('Solar Zenith Angle (θ_z) and Angle of Incidence (θ)');
legend('Solar Zenith Angle (θ_z)', 'Angle of Incidence (θ)');
grid on;

% Display results in a table
IrradianceTable = table(hours', G_T', G_bT', G_dT', G_rT', theta_z', theta', 'VariableNames', ...
    {'Hour', 'G_T', 'G_bT', 'G_dT', 'G_rT', 'Theta_z', 'Theta'});
disp('--- Irradiance Components and Angles ---');
disp(IrradianceTable);

%% 3. PV Array Performance Modeling

% Validate Irradiance Data Length
if length(G_T) ~= 24
    error('Irradiance data G_T must contain exactly 24 hourly values.');
end

% Calculate Sunset Hour Angle (H_s) in degrees
phi_rad = deg2rad(phi);
delta_rad = deg2rad(delta);
Hs_rad = acos(-tan(phi_rad) * tan(delta_rad));
Hs = rad2deg(Hs_rad);
fprintf('Sunset Hour Angle (H_s): %.2f°\n', Hs);

% Extraterrestrial Radiation (H_0) Calculation
H0_J = (24 * 3600 * G_sc / pi) * (1 + 0.033 * cosd((360 * n) / 365)) * ...
       (cosd(phi) * cosd(delta) * sind(Hs) + ...
       (Hs * sind(phi) * sind(delta)) * pi / 180);

H0_kWh = H0_J / 3.6e6; % Convert Joules to kWh
fprintf('Extraterrestrial Radiation (H_0): %.2f kWh/m²/day\n', H0_kWh);

% Clearness Index (K_t) Calculation
Kt = H_actual / H0_kWh;
Kt = max(min(Kt, 0.82), 0); % Constrain K_t between 0 and 0.82
fprintf('Clearness Index (K_t): %.2f\n', Kt);

% Calculate Actual Operating Temperature (T) for each hour
T = ambient_temp + (G_T / 800) * (NOCT - 20);

% Calculate Power Output for Each Hour (W)
P = P_rated * (G_T / 1000) .* (1 + gamma .* (T - T_STC));
P(P < 0) = 0; % Ensure no negative power output

% Total Power Output for All Panels (W)
P_total_W = P * num_panels;

% Convert Power Output to kWh for Each Hour
P_total_kWh = P_total_W / 1000; % kWh/hour

% Calculate Total Power Output over the Day (kWh)
sum_P = sum(P_total_kWh); % Total energy output in kWh/day

% Calculate Total Incident Energy over the Day (kWh)
Incident_energy_kWh = (G_T .* total_area) / 1000; % kWh/hour
sum_Incident = sum(Incident_energy_kWh); % Total incident energy in kWh/day

% Calculate Average Efficiency (Dimensionless)
eta_avg = sum_P / sum_Incident;
fprintf('Average PV Array Efficiency: %.2f%%\n', eta_avg * 100);

% Apply Shading and Soiling Corrections
P_corrected_kWh = P_total_kWh * CF_shading * CF_soiling;

% Sum Corrected Power Output
sum_P_corrected = sum(P_corrected_kWh);

% Updated Average Efficiency with Corrections
eta_avg_corrected = sum_P_corrected / sum_Incident;
fprintf('Average PV Array Efficiency after Shading and Soiling Corrections: %.2f%%\n', eta_avg_corrected * 100);

% Apply Inverter and Wiring Efficiencies
P_system_kWh = P_corrected_kWh * eta_inverter * eta_wiring;

% Sum System Power Output
sum_P_system = sum(P_system_kWh);

% Updated Average Efficiency with All Corrections
eta_avg_system = sum_P_system / sum_Incident;
fprintf('Average PV Array Efficiency after All Corrections: %.2f%%\n', eta_avg_system * 100);

% Apply Degradation Over Time
% Calculate Degraded Efficiency after each year
years_vector = 1:years;
eta_degraded_vector = eta_avg_system * (1 - degradation_rate) .^ years_vector;

% Plot Efficiency Degradation Over Time
figure;
plot(years_vector, eta_degraded_vector * 100, 'r-o', 'LineWidth', 2);
xlabel('Years');
ylabel('Efficiency (%)');
title('PV Array Efficiency Degradation Over Time');
grid on;

% Calculate Annual Energy Production Over Time
annual_energy_production = sum_P_system * 365 * (1 - degradation_rate) .^ (years_vector - 1); % kWh/year

% Plot Annual Energy Production Over Time
figure;
plot(years_vector, annual_energy_production / 1000, 'b-s', 'LineWidth', 2);
xlabel('Years');
ylabel('Annual Energy Production (MWh)');
title('Annual Energy Production Over 25 Years');
grid on;

%% 4. Environmental Impact Assessment

% Annual CO₂ Emissions Reduction
annual_CO2_reduction = annual_energy_production * emission_factor_grid; % kg CO₂/year

% Total CO₂ Emissions Reduction Over Lifetime
total_CO2_reduction = sum(annual_CO2_reduction);

% Display Environmental Impact Results
fprintf('\n--- Environmental Impact Assessment ---\n');
fprintf('Annual CO₂ Emissions Reduction (First Year): %.2f tons\n', annual_CO2_reduction(1) / 1000);
fprintf('Total CO₂ Emissions Reduction Over %d Years: %.2f tons\n', years, total_CO2_reduction / 1000);

% Plot Annual CO₂ Emissions Reduction
figure;
plot(years_vector, annual_CO2_reduction / 1000, 'k-*', 'LineWidth', 2);
xlabel('Years');
ylabel('Annual CO₂ Emissions Reduction (tons)');
title('Annual CO₂ Emissions Reduction Over 25 Years');
grid on;

%% 5. Additional Performance Metrics

% Total Irradiance on the Plane of Array over the Day (kWh/m²/day)
total_irradiance_POA = sum(G_T) / 1000; % Convert from W/m² to kWh/m²

% Theoretical Energy Output under STC (kWh)
theoretical_energy_STC = installed_capacity * total_irradiance_POA / (G_sc / 1000); % Considering STC irradiance

% Performance Ratio (PR)
PR = sum_P_system / theoretical_energy_STC;
fprintf('Performance Ratio (PR): %.2f%%\n', PR * 100);

% Capacity Factor (CF)
total_hours_year = 365 * 24;
CF = (annual_energy_production(1) / 1000) / (installed_capacity * total_hours_year);
fprintf('Capacity Factor (CF): %.2f%%\n', CF * 100);

% Specific Yield (kWh/kW/year)
specific_yield = annual_energy_production(1) / installed_capacity;
fprintf('Specific Yield: %.2f kWh/kW/year\n', specific_yield);

% Energy Yield (kWh/m²/year)
energy_yield = annual_energy_production(1) / total_area;
fprintf('Energy Yield: %.2f kWh/m²/year\n', energy_yield);

% Loss Analysis
% Irradiance Losses
irradiance_loss = sum_P - sum_P_corrected; % Losses due to shading and soiling
% Temperature and System Losses
temperature_loss = sum_P_corrected - sum_P_system; % Losses due to temperature, inverter, and wiring
% Total Losses
total_losses = sum_P - sum_P_system;

% Display Loss Analysis
fprintf('\n--- Loss Analysis ---\n');
fprintf('Irradiance Losses: %.2f kWh/day\n', irradiance_loss);
fprintf('Temperature and System Losses: %.2f kWh/day\n', temperature_loss);
fprintf('Total Losses: %.2f kWh/day\n', total_losses);

% Plot Loss Diagram (Pie Chart)
losses = [irradiance_loss, temperature_loss, sum_P_system];
labels = {'Irradiance Losses', 'Temperature & System Losses', 'Net Energy Output'};
figure;
pie(losses, labels);
title('Loss Analysis of PV System');

%% 6. Total Lifetime Energy Production

total_lifetime_energy = sum(annual_energy_production);
fprintf('\nTotal Lifetime Energy Production over %d years: %.2f MWh\n', years, total_lifetime_energy / 1000);

%% 7. Equivalent Environmental Benefits

% Assume one tree absorbs about 21.77 kg CO₂ per year
trees_planted_equivalent = total_CO2_reduction / (21.77 * years);

% Assume an average car emits about 4.6 metric tons CO₂ per year
cars_removed_equivalent = (total_CO2_reduction / 1000) / (4.6 * years);

fprintf('\nEquivalent Environmental Benefits:\n');
fprintf('- Equivalent Trees Planted: %.0f trees\n', trees_planted_equivalent);
fprintf('- Equivalent Cars Taken Off the Road: %.2f cars\n', cars_removed_equivalent);

%% 8. Summary of Results

% Create a summary table
Results = table({'Average Efficiency'; 'After Shading & Soiling'; 'After System Losses'}, ...
               [eta_avg; eta_avg_corrected; eta_avg_system] * 100, ...
               'VariableNames', {'Calculation_Stage', 'Efficiency_Percentage'});

disp('--- PV Array Efficiency Summary ---');
disp(Results);

% Create a table for Efficiency Degradation Over Time
DegradationTable = table(years_vector', eta_degraded_vector' * 100, annual_energy_production', ...
    annual_CO2_reduction', 'VariableNames', {'Year', 'Efficiency_Percentage', 'Annual_Energy_kWh', 'Annual_CO2_Reduction_kg'});

disp('--- PV Array Efficiency and Annual Energy Production Over Time ---');
disp(DegradationTable);

%% 9. Conclusions

% Display key findings
fprintf('\n--- Key Findings ---\n');
fprintf('Total Daily Energy Generated (After Losses): %.2f kWh\n', sum_P_system);
fprintf('Annual Energy Production (First Year): %.2f MWh\n', sum_P_system * 365 / 1000);
fprintf('Total Number of Panels Required: %d\n', num_panels);
fprintf('Total CO₂ Emissions Reduction Over %d Years: %.2f tons\n', years, total_CO2_reduction / 1000);
fprintf('Performance Ratio (PR): %.2f%%\n', PR * 100);
fprintf('Capacity Factor (CF): %.2f%%\n', CF * 100);
fprintf('Specific Yield: %.2f kWh/kW/year\n', specific_yield);
fprintf('Energy Yield: %.2f kWh/m²/year\n', energy_yield);
fprintf('Total Lifetime Energy Production: %.2f MWh\n', total_lifetime_energy / 1000);
fprintf('Equivalent Trees Planted: %.0f trees\n', trees_planted_equivalent);
fprintf('Equivalent Cars Taken Off the Road: %.2f cars\n', cars_removed_equivalent);

