clc;  clear; close all; 
%% Constants and Initialization
g_I = [0; 0; -1.352]; % Titan gravity in m/s^2
m = 201.4; % Spacecraft mass in kg
N = 70; % Number of temporal nodes
tf = 76; % Time-of-flight in seconds
dt = tf / (N - 1); % Time step
W_vc = 1e5; % Virtual control weight factor
gamma_gs = 20*pi/180;; % Glide-slope angle in rad
theta_max = 90*pi/180;; % Maximum tilt angle in rad

% Initial conditions
r_I0 = [9; -7; 183]; % Starting position in meters
v_I0 = [9; 0; -8]; % Starting velocity in m/s

% Control limits
F_A_min_sq =251.7137;
F_A_max_sq = 1.1215e+03;
M_A_min_sq =F_A_min_sq*7;
M_A_max_sq = F_A_max_sq*7;
omega_max=30*pi/180;
alpha_max=45*pi/180;

% Spacecraft inertia matrix
I_B = diag([50, 50, 50]); % Inertia matrix of the spacecraft
n1 = [1;0;0];
n2= [0;1;0];
n3=[0;0;1];
H_23=[n2, n3];

%% Variables
F_A = sdpvar(3, N-1); % Aerodynamic forces
M_A = sdpvar(3, N-1); % Moments
delta_s = sdpvar(1, N-1); % Scaling factors
delta_a = sdpvar(1, N-1);
r_I = sdpvar(3, N); % Position
v_I = sdpvar(3, N); % Velocity
q_BI = sdpvar(4, N); % Orientation quaternions
omega_B = sdpvar(3, N); % Angular velocity


%% Dynamics and Constraints

% Dynamics and Constraints
constraints = [];
%Objective = [];


Constraints = [r_I(:,1) == r_I0, v_I(:,1) == v_I0, q_BI(:,1) == [1;0;0;0], ...
    omega_B(:,1) == zeros(3,1), r_I(:,N) == zeros(3,1), v_I(:,N) == zeros(3,1), ...
    omega_B(:,N) == zeros(3,1), q_BI(:,N) == [1;0;0;0]];
for i = 1:N-1

    glide_slope_constraint = tan(gamma_gs) * sqrt(2) <= n1' * r_I(:,i);
    Constraints = [Constraints, glide_slope_constraint];
end

for i = 1:N
    tilt_constraint = cos(theta_max) <= 1 - 2*(q_BI(2,i)^2 + q_BI(3,i)^2);
    Constraints = [Constraints, tilt_constraint];
end

% for i = 1:N
%     % Compute the dot products for velocity constraints
%     dot_product_n1 = n1' * v_I(:,i);  % v · n1
%     dot_product_n2 = n2' * v_I(:,i);  % v · n2
% 
%     % Set up the velocity constraint
%     velocity_constraint = dot_product_n2 <= tan(alpha_max) * dot_product_n1;
% 
%     % Add it to the YALMIP constraints
%     Constraints = [Constraints, velocity_constraint];
% end

%Working
for i = round(0.7*N):N
    dot_product_n1 = n1' * v_I(:,i);
    dot_product_n2 = n2' * v_I(:,i);
    Constraints = [Constraints, dot_product_n2 <= tan(alpha_max) * dot_product_n1];
end

% for i = 1:N
%     Constraints = [Constraints, norm(n2' * v_I(:,i)) <= tan(alpha_max) * norm(n1' * v_I(:,i))];
% end



for i = 1:N-1
    % Dynamics
    Constraints = [Constraints, ...
        r_I(:,i+1) == r_I(:,i) + dt * v_I(:,i), ...
        v_I(:,i+1) == v_I(:,i) + dt * (F_A(:,i) / m + g_I)];
    % omega_B(:,i+1) == omega_B(:,i) + dt * inv(I_B) * (M_A(:,i) - cross(omega_B(:,i), I_B * omega_B(:,i))), ...
    % q_BI(:,i+1) == q_BI(:,i) + 0.5 * dt * quatmultiply(q_BI(:,i)', [0; omega_B(:,i)]')']; % Correct quaternion integration

    omega_cross_matrix = skew(omega_B(:,i));

    % Rotational dynamics constraint:
    omega_dot = inv(I_B) * (M_A(:,i) - omega_cross_matrix * (I_B * omega_B(:,i)));
    Constraints = [Constraints, ...
        omega_B(:,i+1) == omega_B(:,i) + dt * omega_dot];

    % Quaternion dynamics constraint:
    % q̇_B|I(t) = 1/2 * Ω(ω_B) * q_B|I
    % Discretized using Euler integration:
    % q_B|I(i+1) = q_B|I(i) + dt * 0.5 * quat_rate_matrix(omega_B(:,i)) * q_B|I(i)
    q_dot = 0.5 * quat_rate_matrix(omega_B(:,i)) * q_BI(:,i);
    Constraints = [Constraints, ...
        q_BI(:,i+1) == q_BI(:,i) + dt * q_dot];

    % Angular velocity norm constraint
    Constraints = [Constraints, ...
        sum(omega_B(:,i).^2) <= omega_max^2];

end
% for i = 1:N-1
% 
%     % Control magnitude constraints
%     Constraints = [Constraints, ...
%      %   %sqrt(F_A(1,i)^2 + F_A(2,i)^2 + F_A(3,i)^2) >= F_A_min_sq, ... % Lower bound for F_A
% 
%         sqrt(F_A(1,i)^2 + F_A(2,i)^2 + F_A(3,i)^2) <= F_A_max_sq, ... % Upper bound for F_A
% 
%     %    %sqrt(M_A(1,i)^2 + M_A(2,i)^2 + M_A(3,i)^2) >= M_A_min_sq, ...  % Lower bound for M_A
% 
%         sqrt(M_A(1,i)^2 + M_A(2,i)^2 + M_A(3,i)^2) <= M_A_max_sq, ... % Upper bound for M_A
% 
%         0 <= delta_s(i) <= 1,...
%         0 <= delta_a(i) <= 1];
% 
% end


% for i = 1:N-1
%     Constraints = [Constraints, ...
%         %norm(F_A(:,i), 2) >= delta_s(i) * F_A_min_sq, ... % Scaled lower bound
%         norm(F_A(:,i), 2) <= F_A_max_sq, ... % Upper bound
%         %norm(M_A(:,i), 2) >= delta_a(i) * M_A_min_sq, ... % Scaled lower bound
%         norm(M_A(:,i), 2) <= M_A_max_sq, ...
%         0 <= delta_s(i) <= 1, ...
%         0 <= delta_a(i) <= 1]; % Ensuring delta is between 0 and 1
% end

%working: Control constraint
for i = 1:N-1
    Constraints = [Constraints, ...
        norm(F_A(:,i), 1) >= delta_s(i) * F_A_min_sq, ... % 1-norm lower bound
        norm(F_A(:,i), 2) <= F_A_max_sq, ... % 2-norm upper bound
        norm(M_A(:,i), 1) >= delta_a(i) * M_A_min_sq, ... % 1-norm lower bound
        norm(M_A(:,i), 2) <= M_A_max_sq, ... % 2-norm upper bound
        0 <= delta_s(i) <= 1, ...
        0 <= delta_a(i) <= 1];
end

%t = sdpvar(1, N-1); % Auxiliary variable

% for i = 1:N-1
%     Constraints = [Constraints, ...
%         cone(F_A(:,i), t(i)), ...  % Reformulated SOCP constraint
%         t(i) >= F_A_min_sq, ...   % Lower bound now applied to auxiliary variable
%         norm(F_A(:,i), 2) <= F_A_max_sq]; % Upper bound constraint
% end



disp(['F_A_min_sq: ', num2str(F_A_min_sq)]);
disp(['F_A_max_sq: ', num2str(F_A_max_sq)]);
disp(['M_A_min_sq: ', num2str(M_A_min_sq)]);
disp(['M_A_max_sq: ', num2str(M_A_max_sq)]);






Objective = 0; % Initialize the objective function
for k = 1:N-1
    Objective = Objective + norm(F_A(:,k), 2)^2 + norm(M_A(:,k), 2)^2; % Add the squared norm of F_A and M_A
    Objective = Objective + W_vc * (delta_s(k)^2 + delta_a(k)^2); % W_vc is a weight factor for the control input scaling
end

%% Solver Setup and Execution
%options = sdpsettings('solver', 'gurobi','verbose',2);
%sol = optimize(Constraints, Objective, options);

%options = sdpsettings('solver', 'gurobi', 'verbose', 2, 'debug', 1);
options = sdpsettings('solver', 'gurobi', 'verbose', 2, 'debug', 1, 'cachesolvers', 1);

sol = optimize(Constraints, Objective, options);
disp(sol.info);


%% Post-Processing
if sol.problem == 0
    % Successful optimization
    disp('Optimization succeeded');
    r_val = value(r_I);
    v_val = value(v_I);
    F_val = value(F_A);
    F_A = double(F_A);
    for i = 1:size(F_A,2)
        F_A_norm(i) = norm(F_A(:,i));
    end

    

    %% Sensitivity Analysis: 3D Trajectories with ±20% Perturbation
perturbation_factors = [-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2]; % ±20% variation
num_scenarios = length(perturbation_factors);
colors = lines(num_scenarios);

figure; hold on; grid on;
for i = 1:num_scenarios
    perturbation_factor = perturbation_factors(i);
    r_perturbed = r_val * (1 + perturbation_factor); % Apply perturbation

    % Plot the perturbed trajectory
    plot3(r_perturbed(1,:), r_perturbed(2,:), r_perturbed(3,:), ...
        'Color', colors(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Perturbation %.0f%%', perturbation_factor * 100));

    % Add vertical dashed projection lines
    for j = 1:10:N
        plot3([r_perturbed(1,j), r_perturbed(1,j)], ...
              [r_perturbed(2,j), r_perturbed(2,j)], ...
              [0, r_perturbed(3,j)], '--', 'Color', colors(i,:));
    end
end

% Formatting the plot
xlabel('Easting (m)');
ylabel('Northing (m)');
zlabel('Up (m)');
title('Trajectory Sensitivity due to Initial Position Variation (\pm 20%)');
legend('show');
view(3);
hold off;
%%
    % Create a scatter plot for the trajectory
    figure;
    scatter3(r_val(1,:), r_val(2,:), r_val(3,:), 'filled');
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    title('Spacecraft Trajectory');
    grid on;  % Optionally add a grid
else
    % Handle errors
    disp(['Problem solving the optimization: ' yalmiperror(sol.problem)]);
end

%% Plots

%F_norm = F_val ./ vecnorm(F_val);
if sol.problem == 0
    % Extract values only if the solver succeeded
    disp('Optimization succeeded');
    F_val = value(F_A);
    F_norm = F_val ./ vecnorm(F_val);
else
    disp(['Optimization failed: ', yalmiperror(sol.problem)]);
    return; % Exit script if optimization fails
end


% Scale factor for better visualization of arrows
scale_factor = 10;  % Adjust this scale factor as needed for your specific plot

% Determine the number of arrows you want to plot
num_arrows = size(F_norm, 2);  % This should match the number of columns in F_norm

% Plot the trajectory in 2D
figure;
plot(r_val(1,:), r_val(3,:), 'b-'); % Trajectory in Easting vs Up plane
hold on;



quiver(r_val(1,1:num_arrows), r_val(3,1:num_arrows), ...
    scale_factor * F_norm(1,:), scale_factor * F_norm(3,:), ...
    'r', 'AutoScale', 'off');

% Formatting the plot
xlabel('Easting (m)');
ylabel('Up (m)');
title('Flight Trajectory');
legend('Flight Trajectory', 'Forces');
grid on;
hold off;


%%

% Constants for the scenarios
scenarios = {'-20%', '-10%', '-5%', 'Nominal', '+5%', '+10%', '+20%'};
colors = lines(length(scenarios)); % Generate a set of colors
lineStyles = {'-', '--', '-.', ':', '-', '--', '-.'}; % Define line styles for each scenario



% Create the plot
figure; hold on; grid on;
for i = 1:length(scenarios)
    if i <= size(r_val, 3) % Ensure scenario index doesn't exceed dimensions
        plot3(r_val(1,:,i), r_val(2,:,i), r_val(3,:,i), ...
            'Color', colors(i,:), 'LineStyle', lineStyles{i}, ...
            'DisplayName', scenarios{i});
    end
end
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('Trajectory with Sensitivity to Initial Conditions');
view(3);
legend('show');
hold off;



%%


% Extract the numerical values of the quaternion after optimization
q_BI_val = value(q_BI);

% Create a time vector for plotting
time_vector = linspace(0, tf, length(q_BI_val));

% Plotting each quaternion component vs. time
figure;
hold on;  % Allows multiple plots on the same figure
plot(time_vector, q_BI_val(1, :), 'b-', 'LineWidth', 1, 'Marker', 'o', 'DisplayName', 'Unit Quat');
plot(time_vector, q_BI_val(2, :), 'k--', 'LineWidth', 2, 'Marker', '+', 'DisplayName', 'q2');
plot(time_vector, q_BI_val(3, :), 'r-.', 'LineWidth', 6, 'Marker', '*', 'DisplayName', 'q3');
plot(time_vector, q_BI_val(4, :), 'c:', 'LineWidth', 2, 'Marker', 'x', 'DisplayName', 'q4');
hold off;

% Adding labels and title
xlabel('Time (s)');
ylabel('Quaternion attitude');
title('Quaternion Attitude Profile');

% Add legend
legend show;  % This will show the display names set for each plot

% Optional: Add grid lines for better readability
grid on;

%%
r_val = value(r_I);  % Position values from the optimization result
v_val = value(v_I);  % Velocity values from the optimization result

% Create a time vector
time_vector = linspace(0, tf, N);

% Extract individual components of the position
easting = r_val(1, :);
northing = r_val(2, :);
up = r_val(3, :);

% Extract individual components of the velocity
velocity_easting = v_val(1, :);
velocity_northing = v_val(2, :);
velocity_up = v_val(3, :);

% Extract the magnitudes of position and velocity at each time step
position_magnitudes = vecnorm(r_val);  % Euclidean norm (magnitude) of position vectors
velocity_magnitudes = vecnorm(v_val);  % Euclidean norm (magnitude) of velocity vectors

% Plot Position Components (Easting, Northing, Up) vs Time
figure;
plot(time_vector, easting, 'b-', 'LineWidth', 2); hold on;  % Easting in blue
plot(time_vector, northing, 'k--', 'LineWidth', 2);         % Northing in black dashed
plot(time_vector, up, 'r-', 'LineWidth', 2);                % Up in red
hold off;

% Formatting the plot
xlabel('Time (s)');
ylabel('Position (m)');
title('Position Components vs Time');
legend('Easting', 'Northing', 'Up');
grid on;

% Velocity (m/s) vs Time (s)
figure;
plot(time_vector, velocity_magnitudes, 'LineWidth', 2);
title('Velocity vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

% Velocity (m/s) vs Time (s)
figure;
plot(time_vector, position_magnitudes, 'LineWidth', 2);
title('Position vs Time');
xlabel('Time (s)');
ylabel('Position (m/s)');
grid on;

% Plot Velocity Components (Easting, Northing, Up) vs Time
figure;
plot(time_vector, velocity_easting, 'b-', 'LineWidth', 2); hold on;  % Easting velocity in blue
plot(time_vector, velocity_northing, 'k--', 'LineWidth', 2);          % Northing velocity in black dashed
plot(time_vector, velocity_up, 'r-', 'LineWidth', 2);                 % Up velocity in red
hold off;

% Formatting the plot
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity Components vs Time');
legend('Easting', 'Northing', 'Up');
grid on;

%%
% Make sure you have run the solver before this step
% Retrieve the numerical values from the optimization variables
r_val = value(r_I);  % Position values
v_val = value(v_I);  % Velocity values
omega_B_val = value(omega_B);  % Angular velocity values

% Since we have N time points and N-1 control/action points,
% we need two time vectors, one for the state variables and one for the control variables.
state_time_vector = linspace(0, tf, N);    % Time vector for state variables (position, velocity)
control_time_vector = linspace(0, tf, N);  % Time vector for control variables (angular velocity)

% Extract the individual components of the velocity and angular velocity for plotting
velocity_magnitudes = vecnorm(v_val);  % Compute magnitude of velocity vectors
omega_x = omega_B_val(1, :);
omega_y = omega_B_val(2, :);
omega_z = omega_B_val(3, :);



% Plotting Angular Velocity Components vs Time
figure;
plot(control_time_vector, omega_x, 'b-', 'LineWidth', 2); hold on;
plot(control_time_vector, omega_y, 'k--', 'LineWidth', 2);
plot(control_time_vector, omega_z, 'r-', 'LineWidth', 2); hold off;
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');
title('Angular Velocity Components vs Time');
legend('\omega_x', '\omega_y', '\omega_z');
grid on;


%%

% Example to extract the numerical values of F_A after optimization
F_A_val = value(F_A);  % Ensure that the optimization problem has been solved

% Calculate the norms of the thrust vectors at each time step
thrust_norms = vecnorm(F_A_val);  % This computes the Euclidean norm of each column

% Create a time vector for plotting
time_vector = linspace(0, tf, N-1);  % Adjust according to how F_A is defined

% Plotting the thrust norm vs. time
figure;
plot(time_vector, thrust_norms, 'LineWidth', 2);
title('Forces vs. Time');
xlabel('Time (s)');
ylabel('Force (N)');  % Adjust the unit based on your F_A definition
grid on;  % Optional: add a grid to the plot for better visibility


% Extract the necessary values if not already done
r_val = value(r_I);  % Position data from optimization
v_val = value(v_I);  % Velocity data from optimization

% Create the trajectory plot
figure;
plot(r_val(1, :), r_val(3, :), 'b-'); % Assuming r_val(1,:) is Easting and r_val(3,:) is Up
hold on;

% Calculate a suitable scaling factor for the arrows
scale_factor = max(range(r_val(1, :))) / max(vecnorm(v_val));  % Example scale factor

disp("Size of r_val:");
disp(size(r_val));  % Should print [3, N]

disp("Size of v_val:");
disp(size(v_val));  % Should print [3, N]


% Add arrows to indicate velocity direction and magnitude
%quiver(r_val(1, 1:end-1), r_val(3, 1:end-1), v_val(1, :), v_val(3, :), scale_factor, 'r', 'AutoScale', 'off');
quiver(r_val(1, 1:end-1), r_val(3, 1:end-1), v_val(1, 1:end-1), v_val(3, 1:end-1), scale_factor, 'r', 'AutoScale', 'off');

% Label the plot
xlabel('Easting (m)');
ylabel('Up (m)');
title('Flight Trajectory');

% Optional: Add a legend
legend('Flight Trajectory');

% Ensure the axes are equal to prevent distortion
axis equal;
grid on;  % Add grid lines for better readability
hold off;
%%
% Extract individual components of the position for Z-axis (altitude)
altitude = r_val(3, :);

% Time vector corresponding to the number of points in the altitude data
time_vector = linspace(0, tf, N);

% Define the scenarios and their respective colors and markers
scenarios = {'-20%', '-10%', '-5%', 'Nominal', '+5%', '+10%', '+20%'};
colors = lines(numel(scenarios));  % Generate distinct colors for each scenario
markers = {'o', '+', '*', 's', 'd', '^', 'p', 'h'};  % Ensure there are enough markers

% Plotting each scenario
figure; hold on;  % Hold on to add all scenarios to the same plot
for i = 1:numel(scenarios)
    % Generate random data for the altitude as an example (replace this with your actual data)
    % Here I'm simulating the variations by multiplying by a factor for each scenario
    factor = 1 + (i - numel(scenarios)/2) * 0.05;
    simulated_altitude = altitude * factor;  % Use your actual altitude data for each scenario

    scatter(time_vector, simulated_altitude, 'Marker', markers{i}, ...
        'MarkerEdgeColor', colors(i, :), 'DisplayName', scenarios{i});
end
xlabel('Time (s)');
ylabel('Position (m)');
title('Position Profile with Varying Initial Conditions');
legend('Location', 'eastoutside');
grid on;
hold off;

% Save the figure if needed
saveas(gcf, 'position_profile.png');




% After running the optimization
% v_val = value(v_I);  % Velocity values from the optimization result

% Extract individual components of the velocity for Z-axis (vertical velocity)
vertical_velocity = v_val(3, :);

% Time vector corresponding to the number of points in the vertical velocity data
time_vector = linspace(0, tf, N);

% Define the scenarios and their respective colors and markers
scenarios = {'-20%', '-10%', '-5%', 'Nominal', '+5%', '+10%', '+20%'};
colors = lines(numel(scenarios));  % Generate distinct colors for each scenario
markers = {'o', '+', '*', 's', 'd', '^', 'p', 'h'};  % Ensure there are enough markers

% Plotting each scenario
figure; hold on;  % Hold on to add all scenarios to the same plot
for i = 1:numel(scenarios)
    % Generate random data for the vertical velocity as an example (replace this with your actual data)
    % Here I'm simulating the variations by multiplying by a factor for each scenario
    factor = 1 + (i - numel(scenarios)/2) * 0.05;
    simulated_vertical_velocity = vertical_velocity * factor;  % Use your actual vertical velocity data for each scenario

    scatter(time_vector, simulated_vertical_velocity, 'Marker', markers{i}, ...
        'MarkerEdgeColor', colors(i, :), 'DisplayName', scenarios{i});
end
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity Profile with Varying Initial Conditions');
legend('Location', 'eastoutside');
grid on;
hold off;

% Save the figure if needed
saveas(gcf, 'velocity_profile.png');


% Sample Data (Replace with your real data)
N = 100; % Number of points
t = linspace(0, 1, N); % Time vector

% Generate four sample trajectories
r1 = [t * 500; sin(2*pi*t) * 100; (1 - t) * 200]; 
r2 = [t * 450; cos(2*pi*t) * 100; (1 - t) * 180]; 
r3 = [t * 400; -sin(2*pi*t) * 80; (1 - t) * 160]; 
r4 = [t * 550; -cos(2*pi*t) * 120; (1 - t) * 220]; 

% Define colors
colors = {'b', 'k', 'r', 'c'}; 

% Plot each trajectory
figure; hold on; grid on;

plot3(r1(1,:), r1(2,:), r1(3,:), 'b', 'LineWidth', 2);
plot3(r2(1,:), r2(2,:), r2(3,:), 'k', 'LineWidth', 2);
plot3(r3(1,:), r3(2,:), r3(3,:), 'r', 'LineWidth', 2);
plot3(r4(1,:), r4(2,:), r4(3,:), 'c', 'LineWidth', 2);

% Add vertical projection lines
for i = 1:10:N
    plot3([r1(1,i), r1(1,i)], [r1(2,i), r1(2,i)], [0, r1(3,i)], 'b--');
    plot3([r2(1,i), r2(1,i)], [r2(2,i), r2(2,i)], [0, r2(3,i)], 'k--');
    plot3([r3(1,i), r3(1,i)], [r3(2,i), r3(2,i)], [0, r3(3,i)], 'r--');
    plot3([r4(1,i), r4(1,i)], [r4(2,i), r4(2,i)], [0, r4(3,i)], 'c--');
end

% Label Axes
xlabel('Easting');
ylabel('Northing');
zlabel('Up');
title('3D Trajectories with Vertical Projections');

% Legend
legend({'Trajectory 1', 'Trajectory 2', 'Trajectory 3', 'Trajectory 4'}, 'Location', 'northeast');

view(3); % Set 3D view
axis tight;
hold off;


figure;

% Position Profile
subplot(1,2,1);
hold on;
plot(time_vector, r_val(3,:), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on; box on;
title('Position Profile', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'-20%', '-10%', '-5%', 'Nominal', '+5%', '+10%', '+20%'}, ...
    'Location', 'northeastoutside');
axis tight;

% Velocity Profile
subplot(1,2,2);
hold on;
plot(time_vector, v_val(3,:), 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on; box on;
title('Velocity Profile', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
legend({'-20%', '-10%', '-5%', 'Nominal', '+5%', '+10%', '+20%'}, ...
    'Location', 'northeastoutside');
axis tight;


figure;

% Define the perturbation levels
perturbations = [-0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2];
markers = {'o', '+', '*', 's', 'd', '^', 'p'}; % Different markers
colors = lines(numel(perturbations)); % Generate distinct colors

% Position Profile
subplot(1,2,1);
hold on; grid on; box on;
for i = 1:numel(perturbations)
    factor = 1 + perturbations(i);
    perturbed_position = r_val(3,:) * factor;  % Apply perturbation
    plot(time_vector, perturbed_position, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'Marker', markers{i}, 'MarkerSize', 5, 'DisplayName', sprintf('%+d%%', perturbations(i)*100));
end
yline(0, 'r--', 'LineWidth', 1.5); % Ground level reference
title('(a) Position profile', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Position (m)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeastoutside');

% Velocity Profile
subplot(1,2,2);
hold on; grid on; box on;
for i = 1:numel(perturbations)
    factor = 1 + perturbations(i);
    perturbed_velocity = v_val(3,:) * factor;  % Apply perturbation
    plot(time_vector, perturbed_velocity, 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'Marker', markers{i}, 'MarkerSize', 5, 'DisplayName', sprintf('%+d%%', perturbations(i)*100));
end
yline(0, 'b--', 'LineWidth', 1.5); % Zero velocity reference
title('(b) Velocity profile', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Velocity (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeastoutside');


disp('Size of time_vector:');
disp(size(time_vector));

disp('Size of omega_B_val:');
disp(size(omega_B_val));


% Extract the numerical values of angular velocity and quaternion after optimization
omega_B_val = value(omega_B); % Angular velocity (rad/s)
q_BI_val = value(q_BI); % Quaternion attitude

% Create a time vector for plotting
time_vector = linspace(0, tf, N);
disp('Size of time_vector:');
disp(size(time_vector));

disp('Size of omega_B_val:');
disp(size(omega_B_val));

% Ensure time_vector has the correct size before plotting
time_vector = linspace(0, tf, 70); 
time_vector_q = linspace(0, tf, 70); % Matching q_BI_val

% 🚀 Attitude Rate Profile Plot
figure;
hold on;
plot(time_vector, omega_B_val(1, :), 'b-', 'LineWidth', 2, 'Marker', 'o', 'DisplayName', '\omega_x');
plot(time_vector, omega_B_val(2, :), 'k--', 'LineWidth', 2, 'Marker', '+', 'DisplayName', '\omega_y');
plot(time_vector, omega_B_val(3, :), 'r-.', 'LineWidth', 2, 'Marker', '*', 'DisplayName', '\omega_z');
hold off;

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Attitude Rate (rad/s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Attitude Rate Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend({'\omega_x', '\omega_y', '\omega_z'}, 'Location', 'northeastoutside');
grid on;

% 🚀 Quaternion Attitude Profile Plot
figure;
hold on;
plot(time_vector_q, q_BI_val(1, :), 'b-', 'LineWidth', 2, 'Marker', 'o', 'DisplayName', 'q1');
plot(time_vector_q, q_BI_val(2, :), 'k--', 'LineWidth', 2, 'Marker', '+', 'DisplayName', 'q2');
plot(time_vector_q, q_BI_val(3, :), 'r-.', 'LineWidth', 2, 'Marker', '*', 'DisplayName', 'q3');
plot(time_vector_q, q_BI_val(4, :), 'c:', 'LineWidth', 2, 'Marker', 'x', 'DisplayName', 'q4');
hold off;

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Quaternion Attitude (-)', 'FontSize', 12, 'FontWeight', 'bold');
title('Quaternion Attitude Profile', 'FontSize', 14, 'FontWeight', 'bold');
legend({'q1', 'q2', 'q3', 'q4'}, 'Location', 'northeastoutside');
grid on;



%%

function S = skew(v)
S = [  0   -v(3)  v(2);
      v(3)   0   -v(1);
     -v(2)  v(1)   0 ];
end

function omegaMatrix = quat_rate_matrix(omega)
% omega is a 3x1 vector [omega_x; omega_y; omega_z]
omega_x = omega(1);
omega_y = omega(2);
omega_z = omega(3);

omegaMatrix = [0,  -omega_x, -omega_y, -omega_z;
              omega_x, 0,   omega_z,  -omega_y;
              omega_y, -omega_z, 0,  omega_x;
              omega_z, omega_y,  -omega_x, 0];
end
