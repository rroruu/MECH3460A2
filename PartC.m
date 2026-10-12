clc; close all;

%% SYSTEM LOADS
fprintf('========================================================================\n');
fprintf('        PART C: WELLED JOINT ANALYSIS & OPTIMIZATION (PLATE 1 TO BEAM)  \n');
fprintf('========================================================================\n\n');

% External Load & Bench Properties 
F_ext_max = 3500;       
V_bench = 0.06628;      
rho_steel = 7800;       
g = 9.81;               
W_bench = V_bench * rho_steel * g; 

x_W = 0.2955;           % Distance to bench CG along x-axis (m)
x_beam = 0.250;          % Distance to support beam CG along x-axis (m)
n_req = 1.5;            % Target minimum design safety factor

%% LOAD PART A DATA
mat_filename = 'PartA_selected_section.mat';

if exist(mat_filename, 'file') == 2
    partA_data = load(mat_filename);
    
    b         = partA_data.best(2);        % Section width (mm)
    d         = partA_data.best(3);        % Section depth/height (mm)
    t_beam    = partA_data.best(4);        % Wall thickness (mm)
    mass_beam = partA_data.best(5);        % Beam mass (kg)
    
    g         = partA_data.P.g;            % 9.81 m/s^2
    L1        = partA_data.P.L1;           % m
    L2        = partA_data.P.L2;           % m
    L_beam    = partA_data.P.L;            % 0.500 m
    W_beam    = mass_beam * g;             % Beam weight (N)
    
    fprintf('Loaded section properties from Part A:\n');
    fprintf('  Size: %.0f x %.0f x %.1f mm\n', b, d, t_beam);
    fprintf('  Mass: %.2f kg\n\n', mass_beam);
else
    fprintf('Warning: %s not found. Using standard default parameters.\n\n', mat_filename);
    b         = 100;      % mm
    d         = 100;      % mm
    t_beam    = 5;        % mm
    mass_beam = 7.41;     % kg
    g         = 9.81;     % m/s^2
    W_beam    = mass_beam * g;
    L1        = 0.150;    % m
    L2        = 0.830;    % m
    L_beam    = 0.500;    % m
end

x_F = 0.120 + 0.035 + L1 + 0.015; % Distance along x-axis to F_ext (m)[cite: 5]

%% WELD PROFILE PROPERTIES 
A_u = 2 * (b + d);                      % Unit weld area (mm)
I_u = (d^2 / 6) * (3*b + d);          % Unit second moment of area (mm^3)
J_u = (b + d)^3 / 6;                  % Unit polar moment of area (mm^3)

c = d / 2;                            % Distance for bending stress (mm)
r_corner = sqrt((b/2)^2 + (d/2)^2);   % Radius to corner for torsional shear (mm)

%% ELECTRODE MATERIAL 
S_ut = 482; % Ultimate tensile strength (MPa)
S_y = 393;  % Yield strength (MPa)[cite: 8]

ka = 272 * (S_ut^-0.995);             
kb = 0.85;                           
kc = 0.814;                           % Reliability factor (99% reliability)
kd = 1.0;                             % Temperature factor
ke = 1.0;                             % Miscellaneous effects

S_e_prime = 0.5 * S_ut;               
S_e = S_e_prime * ka * kb * kc * kd * ke; 

Vz_min = W_bench + W_beam;
My_min = (W_bench * x_W) + (W_beam * x_beam); % N-m
My_min_Nmm = My_min * 1000;
Tx_min_Nmm = 0;

Vz_max = Vz_min + F_ext_max;
My_max = My_min + (F_ext_max * x_F);
My_max_Nmm = My_max * 1000;
Tx_max = F_ext_max * L2;
Tx_max_Nmm = Tx_max * 1000;

h_candidates = [3, 4, 5, 6, 8, 10];   % Candidate weld leg sizes (mm)
results = [];                         % Results storage matrix

for i = 1:length(h_candidates)
    h = h_candidates(i);
    t_e = 0.707 * h;                  % Effective throat thickness (mm)
    
    tau_V_min = Vz_min / (t_e * A_u);
    sigma_M_min = (My_min_Nmm * c) / (t_e * I_u);
    vm_min = sqrt(sigma_M_min^2 + 3 * (tau_V_min^2));
    
    tau_V_max = Vz_max / (t_e * A_u);
    sigma_M_max = (My_max_Nmm * c) / (t_e * I_u);
    
    tau_T_max = (Tx_max_Nmm * r_corner) / (t_e * J_u);
    tau_Ty_max = tau_T_max * ((d/2) / r_corner); % Horizontal component (y)
    tau_Tz_max = tau_T_max * ((b/2) / r_corner); % Vertical component (z)
    
    tau_z_total = tau_V_max + tau_Tz_max;
    tau_y_total = tau_Ty_max;
    
    % Von Mises Equivalent Stress on Weld Throat (MPa)
    vm_max = sqrt(sigma_M_max^2 + 3 * (tau_y_total^2 + tau_z_total^2));
    
    sigma_a = (vm_max - vm_min) / 2;   % Alternating stress amplitude
    sigma_m = (vm_max + vm_min) / 2;   % Mean stress
    
    % Static Yield Safety Factor
    n_y = S_y / vm_max;
    
    % Modified Goodman Line
    n_goodman = 1 / ((sigma_a / S_e) + (sigma_m / S_ut));
    
    % Gerber Parabola
    a_gerb = sigma_a / S_e; 
    b_gerb = (sigma_m / S_ut)^2;
    n_gerber = 0.5 * (1 / b_gerb) * a_gerb * (-1 + sqrt(1 + (2 * sigma_m * S_e / (S_ut * sigma_a))^2));
    
    % ASME Ellipse
    n_asme = 1 / sqrt((sigma_a / S_e)^2 + (sigma_m / S_y)^2);
    
    results = [results; h, t_e, vm_min, vm_max, sigma_m, sigma_a, n_y, n_goodman, n_gerber, n_asme];
end

%% REQUIREMENT
h_max_limit = 5.0; 

valid_idx = find(results(:, 9) >= n_req & results(:, 1) <= h_max_limit);

if isempty(valid_idx)
    error('No weld sizes satisfy the minimum safety factor under the 5mm maximum limit.');
end

% Select the smallest valid size
opt_row = results(valid_idx(1), :);
h_opt = opt_row(1);

if h_opt <= 6.0
    pass_type = 'Single-Pass Fillet Weld';
    rec_electrode_dia = 3.2; 
else
    pass_type = 'Multi-Pass Fillet Weld (2-3 passes)';
    rec_electrode_dia = 4.0; 
end

%% OUTPUT RESULTS
fprintf('----------------------------------------------------------------------------------------------------\n');
fprintf(' Leg h | Throat | Max Stress | Mean Str | Alt Str | Yield SF | Goodman SF | Gerber SF | ASME Ell. SF \n');
fprintf(' (mm)  | t_e(mm)| vm_max(MPa)|  (MPa)   |  (MPa)  |   (n_y)  | (n_goodm)  | (n_gerber)|  (n_asme)   \n');
fprintf('----------------------------------------------------------------------------------------------------\n');
for i = 1:size(results, 1)
    marker = ' ';
    if results(i, 1) == h_opt, marker = '*'; end
    fprintf('%c %2d   |  %4.2f  |   %6.1f   |  %6.1f  |  %5.1f  |   %5.2f  |   %5.2f    |   %5.2f   |   %5.2f\n', ...
        marker, results(i,1), results(i,2), results(i,4), results(i,5), results(i,6), ...
        results(i,7), results(i,8), results(i,9), results(i,10));
end
fprintf('----------------------------------------------------------------------------------------------------\n');
fprintf(' * Denotes Selected Optimal Weld Leg Size (h = %d mm)\n\n', h_opt);
fprintf('========================================================================\n');
fprintf('      PART C(h): SMAW / STICK ELECTRODE SELECTION SPECIFICATION        \n');
fprintf('========================================================================\n');
fprintf(' Nominated Process:          Shielded Metal Arc Welding (SMAW / Stick)\n');
fprintf(' Recommended Electrode:     E7018 (AWS A5.1) / E4818 (AS/NZS 1553.1)\n');
fprintf(' Nominated Commercial Size:  %.1f mm Diameter Core Wire Rod\n', rec_electrode_dia);
fprintf(' Weld Pass Configuration:   %s\n', pass_type);
fprintf(' Maximum Allowable Leg:     %.1f mm (Hard Constraint)\n', h_max_limit);
fprintf(' Selected Leg Size (h):     %.1f mm (Effective Throat t_e = %.2f mm)\n', h_opt, opt_row(2));
fprintf(' Gerber Factor of Safety:    %.2f (>= %.1f Design Target)\n', opt_row(9), n_req);
fprintf(' Goodman Factor of Safety:   %.2f\n', opt_row(8));
fprintf(' Yield Factor of Safety:    %.2f\n', opt_row(7));
fprintf('========================================================================\n\n');

%% SAVE RESULTS
save('PartC_welded_joint.mat', 'h_opt', 'opt_row', 'rec_electrode_dia', ...
    'S_ut', 'S_y', 'S_e', 'results', 'F_ext_max', 'W_bench', 'W_beam');

%% PROFESSIONAL VISUALIZATION
% FATIGUE DIAGRAM 
figure('Name', 'Part C - Weld Fatigue Diagram', 'Color', 'w', 'Position', [100, 100, 800, 550]);
hold on; grid on; box on;

sm_plot = linspace(0, S_ut, 200);
sa_gerber = S_e .* (1 - (sm_plot ./ S_ut).^2);
sa_goodman = S_e .* (1 - (sm_plot ./ S_ut));
sa_asme = S_e .* sqrt(1 - (sm_plot ./ S_y).^2);
sa_asme(sm_plot > S_y) = NaN;
sa_langer = S_y - sm_plot; sa_langer(sa_langer < 0) = NaN;

plot(sm_plot, sa_gerber, 'g-', 'LineWidth', 2.0, 'DisplayName', 'Gerber Parabola');
plot(sm_plot, sa_goodman, 'b--', 'LineWidth', 1.8, 'DisplayName', 'Goodman Line');
plot(sm_plot, sa_asme, 'c-.', 'LineWidth', 1.8, 'DisplayName', 'ASME Ellipse');
plot(sm_plot, sa_langer, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Langer Yield Line');

for i = 1:size(results, 1)
    if results(i, 1) == h_opt
        plot(results(i,5), results(i,6), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 10, ...
            'DisplayName', sprintf('Optimal Weld (h = %d mm)', h_opt));
        plot([0, results(i,5)], [0, results(i,6)], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    else
        plot(results(i,5), results(i,6), 'ko', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 6, ...
            'HandleVisibility', 'off');
        text(results(i,5)+2, results(i,6)+1, sprintf('%d mm', results(i,1)), 'FontSize', 8);
    end
end

title('Part C: Fatigue Failure Criteria & Operating Points (E7018 Weld)');
xlabel('Mean Stress, \sigma_m (MPa)'); ylabel('Alternating Stress Amplitude, \sigma_a (MPa)');
xlim([0, S_ut + 20]); ylim([0, S_e + 20]);
legend('Location', 'northwest');

% SAFETY FACTORS VS WELD LEG SIZE
figure('Name', 'Part C - Safety Factor Evaluation', 'Color', 'w', 'Position', [920, 100, 750, 550]);
hold on; grid on; box on;

plot(results(:,1), results(:,8), 'b-s', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'Goodman SF');
plot(results(:,1), results(:,9), 'g-o', 'LineWidth', 2.0, 'MarkerSize', 8, 'DisplayName', 'Gerber SF');
plot(results(:,1), results(:,10), 'c-^', 'LineWidth', 1.8, 'MarkerSize', 7, 'DisplayName', 'ASME Ellipse SF');
plot(results(:,1), results(:,7), 'r--d', 'LineWidth', 1.5, 'MarkerSize', 7, 'DisplayName', 'Static Yield SF');

yline(n_req, 'k--', sprintf('Target n_{req} = %.1f', n_req), 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'bottom');

xline(h_opt, 'm:', sprintf('Optimal h = %d mm', h_opt), 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'top');

title('Part C: Design Safety Factors vs. Weld Leg Size (h)');
xlabel('Weld Leg Size, h (mm)'); ylabel('Factor of Safety, n');
xticks(h_candidates);
ylim([0, max(results(:,7)) + 1]);
legend('Location', 'northwest');

font_name = 'Times New Roman';
font_size = 12; 
line_w    = 3.0;

fig_pos = [100, 100, 600, 450]; 

% FATIGUE DIAGRAM 
c_blue = [0.00 0.45 0.74]; c_green = [0.00 0.55 0.30];
c_purple = [0.55 0.20 0.65]; c_red = [0.80 0.10 0.10]; c_pink = [0.80 0.10 0.60];

fig1 = figure('Name', 'Part C - Weld Fatigue Diagram', 'Color', 'w', ...
              'Units', 'pixels', 'Position', fig_pos, 'Visible', 'off');
fig1.PaperPositionMode = 'auto';
hold on; grid on; box on;

sm_plot = linspace(0, S_ut, 200);
sa_gerber = S_e .* (1 - (sm_plot ./ S_ut).^2);
sa_goodman = S_e .* (1 - (sm_plot ./ S_ut));
sa_asme = S_e .* sqrt(1 - (sm_plot ./ S_y).^2);
sa_asme(sm_plot > S_y) = NaN;
sa_langer = S_y - sm_plot; sa_langer(sa_langer < 0) = NaN;

plot(sm_plot, sa_gerber, '--', 'Color', c_green, 'LineWidth', line_w, 'DisplayName', 'Gerber Parabola');
plot(sm_plot, sa_goodman, '--', 'Color', c_blue, 'LineWidth', line_w, 'DisplayName', 'Goodman Line');
plot(sm_plot, sa_asme, '--', 'Color', c_purple, 'LineWidth', line_w, 'DisplayName', 'ASME Ellipse');
plot(sm_plot, sa_langer, '-', 'Color', c_red, 'LineWidth', line_w, 'DisplayName', 'Yield Line');

for i = 1:size(results, 1)
    if results(i, 1) == h_opt
        plot(results(i,5), results(i,6), 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 11, ...
            'DisplayName', sprintf('Optimal (h = %d mm)', h_opt));
        plot([0, results(i,5)], [0, results(i,6)], 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        text(results(i,5)+10, results(i,6), sprintf('h=%d mm', h_opt), ...
            'FontSize', 18, 'FontName', font_name, 'Interpreter', 'tex', ...
            'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    else
        plot(results(i,5), results(i,6), 'ko', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerSize', 7, ...
            'HandleVisibility', 'off');
        text(results(i,5)+8, results(i,6), sprintf('%d mm', results(i,1)), ...
            'FontSize', 16, 'FontName', font_name, 'Interpreter', 'tex', ...
            'VerticalAlignment', 'middle');
    end
end
xlabel('Mean Stress, \sigma_m (MPa)', 'Interpreter', 'tex'); 
ylabel('Alternating Stress, \sigma_a (MPa)', 'Interpreter', 'tex');
xlim([0, S_ut * 1.25]); ylim([0, S_e * 1.25]);

legend('Location', 'northeast', 'Interpreter', 'tex', 'FontSize', 16, 'NumColumns', 1);
set(gca, 'FontName', font_name, 'FontSize', 18, 'XMinorTick', 'on', 'YMinorTick', 'on');
exportgraphics(fig1, 'PartC_Fatigue_Diagram.png', 'Resolution', 300);
close(fig1);

% SAFETY FACTORS 
fig2 = figure('Name', 'Part C - Safety Factor Evaluation', 'Color', 'w', ...
              'Units', 'pixels', 'Position', fig_pos, 'Visible', 'off');
fig2.PaperPositionMode = 'auto';
hold on; grid on; box on;

plot(results(:,1), results(:,8), '-s', 'Color', c_blue, 'LineWidth', line_w, 'MarkerSize', 8, 'DisplayName', 'Goodman SF');
plot(results(:,1), results(:,9), '-o', 'Color', c_green, 'LineWidth', line_w, 'MarkerSize', 8, 'DisplayName', 'Gerber SF');
plot(results(:,1), results(:,10), '-^', 'Color', c_purple, 'LineWidth', line_w, 'MarkerSize', 8, 'DisplayName', 'ASME Ellipse');

yline(n_req, 'k--', sprintf('Target n = %.1f', n_req), 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 16, 'FontName', font_name);

xline(4, ':', 'h = 4 mm', 'Color', c_pink, 'LineWidth', 2.0, ...
    'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'right', 'FontSize', 16, 'FontName', font_name);

xlabel('Weld Leg Size, h (mm)', 'Interpreter', 'tex'); 
ylabel('Factor of Safety, n', 'Interpreter', 'tex');
xticks(h_candidates); 

xlim([2,11]);
max_sf = max(max(results(:, 8:10)));
ylim([0, max_sf * 1.25]);

legend('Goodman SF', 'Gerber SF', 'ASME Ellipse', ...
       'Location', 'north', ...
       'Interpreter', 'tex', ...
       'FontSize', 16, ...
       'NumColumns', 1);
set(gca, 'FontName', font_name, 'FontSize', 18, 'XMinorTick', 'on', 'YMinorTick', 'on');
exportgraphics(fig2, 'PartC_Safety_Factors.png', 'Resolution', 300);
close(fig2);

% WELD PROFILE
fig3 = figure('Name', 'Part C - Weld Profile Comparison', 'Color', 'w', ...
              'Units', 'pixels', 'Position', fig_pos, 'Visible', 'off');
fig3.PaperPositionMode = 'auto';
hold on; grid on; axis equal; axis padded; box on;
rectangle('Position', [-b/2, -d/2, b, d], 'EdgeColor', 'k', 'LineWidth', line_w);
rectangle('Position', [-b/2+t_beam, -d/2+t_beam, b-2*t_beam, d-2*t_beam], 'EdgeColor', 'k', 'LineStyle', '--', 'LineWidth', 1.2);
rectangle('Position', [-b/2-h_opt, -d/2-h_opt, b+2*h_opt, d+2*h_opt], 'EdgeColor', '#D95319', 'LineStyle', '-.', 'LineWidth', line_w);

p1 = plot(NaN, NaN, 'k-', 'LineWidth', line_w);
p2 = plot(NaN, NaN, 'k--', 'LineWidth', 1.2);
p3 = plot(NaN, NaN, 'Color', '#D95319', 'LineStyle', '-.', 'LineWidth', line_w);
title('Weld Profile (Plate 1 to Beam)', 'Interpreter', 'tex', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('Width (mm)', 'Interpreter', 'tex'); ylabel('Depth (mm)', 'Interpreter', 'tex');
legend([p1, p2, p3], {'Outer Boundary', 'Inner Boundary', sprintf('Weld (h = %.1f mm)', h_opt)}, ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 11, 'Interpreter', 'tex');
set(gca, 'FontName', font_name, 'FontSize', font_size);
exportgraphics(fig3, 'PartC_Weld_Profile.png', 'Resolution', 300);
close(fig3);

% FORCE COMPONENTS 
fig4 = figure('Name', 'Part C - Force Components on Weld Outline', 'Color', 'w', ...
              'Units', 'pixels', 'Position', fig_pos, 'Visible', 'off');
fig4.PaperPositionMode = 'auto';
hold on; grid on; axis equal; axis padded; view(35, 30); box on;
x_w = zeros(1,5);
y_w = [-b/2, b/2, b/2, -b/2, -b/2];
z_w = [d/2, d/2, -d/2, -d/2, d/2];
plot3(x_w, y_w, z_w, 'k-', 'LineWidth', line_w, 'DisplayName', 'Weld Outline');
scale_F = max(b, d) * 0.6 / max(Vz_max, 1); 
scale_M = max(b, d) * 0.6 / max(My_max, 1);
quiver3(0, 0, 0, 0, 0, -Vz_max * scale_F, 0, 'Color', 'b', 'LineWidth', line_w, 'MaxHeadSize', 0.5, 'DisplayName', 'V_z (Shear)');
text(0, 0, -Vz_max * scale_F - 15, 'V_z', 'Color', 'b', 'FontWeight', 'bold', 'FontName', font_name, 'FontSize', 13, 'Interpreter', 'tex');
quiver3(0, 0, 0, 0, My_max * scale_M, 0, 0, 'Color', 'g', 'LineWidth', line_w, 'MaxHeadSize', 0.5, 'DisplayName', 'M_y (Moment)');
text(0, My_max * scale_M + 15, 0, 'M_y', 'Color', 'g', 'FontWeight', 'bold', 'FontName', font_name, 'FontSize', 13, 'Interpreter', 'tex');
quiver3(0, 0, 0, Tx_max * scale_M * 5, 0, 0, 0, 'Color', '#D95319', 'LineWidth', line_w, 'MaxHeadSize', 0.5, 'DisplayName', 'T_x (Torsion)');
text(Tx_max * scale_M * 5 + 15, 0, 0, 'T_x', 'Color', '#D95319', 'FontWeight', 'bold', 'FontName', font_name, 'FontSize', 13, 'Interpreter', 'tex');
title('Critical Force and Moment Components', 'Interpreter', 'tex', 'FontSize', 15, 'FontWeight', 'bold');
xlabel('X (Long.)', 'Interpreter', 'tex'); ylabel('Y (Trans.)', 'Interpreter', 'tex'); zlabel('Z (Vert.)', 'Interpreter', 'tex');
legend('Location', 'northeast', 'Interpreter', 'tex', 'FontSize', 11, 'NumColumns', 1);
set(gca, 'FontName', font_name, 'FontSize', font_size);
exportgraphics(fig4, 'PartC_Force_Components.png', 'Resolution', 300);
close(fig4);

fprintf('\nSuccess: Exported 4 high-resolution, LaTeX-ready diagrams cleanly.\n');