clear; clc; close all;

%% LOADS
g = 9.81;  E = 210e9;
if isfile('PartA_selected_section.mat')
    PA = load('PartA_selected_section.mat');
    W_bench = PA.P.W;  x_W = PA.P.xW;  x_F = PA.P.xF;  L2 = PA.P.L2;  F_ext = PA.P.F;
    L_beam  = PA.P.L;  m_beam = PA.best(5);                    
    fprintf('Loads linked from Part A (beam %g x %g x %g mm, %.2f kg)\n', PA.best(2), PA.best(3), PA.best(4), m_beam);
else
    W_bench = 0.06628*7800*g;  x_W = 0.2955;  x_F = 0.650;  L2 = 0.830;  F_ext = 3500;
    L_beam  = 0.500;           m_beam = 7.41;                  
    fprintf('Part A file not found - using default loads.\n');
end
W_beam = m_beam*g;  x_beam = L_beam/2;

%% ASSUMPTIONS
B.t_plate  = 0.020;        % Plate 1 thickness (m)
B.w        = 0.300;        % Plate 1 width (m)
B.h        = 0.300;        % Plate 1 height (m)
B.edge     = 0.040;        
B.frames   = 0.200;        % building-frame thickness (m)
B.joint    = 'through';    
B.washer   = 0.002;        % washer thickness (m)
B.engage   = 2.0;          
B.L_stock  = 0.300;        % longest stocked bolt (m)
B.frac     = 0.75;         % preload Fi = frac * Sp * At
B.crit     = 'goodman';    
B.Sy_plate = 220.594;      % plate yield strength (MPa)
B.E        = E;
B.n_req    = 1.5;          
B.fs_def   = 'lecture';    
B.mu_thread = 0.4;         

n_list     = [4 6 8 10 12];
bolt_sizes = [10 58.0 1.5; 12 84.3 1.75; 14 115 2.0; 16 157 2.0; 20 245 2.5; 24 353 3.0];
classes    = [ 8.8  600  830  660  129;
              10.9  830 1040  940  162;
              12.9  970 1220 1100  190];

%% LOAD STATES 
LD.Vmin = W_bench + W_beam;
LD.Mmin = W_bench*x_W + W_beam*x_beam + LD.Vmin*B.t_plate;
LD.Tmin = 0;
LD.Vmax = LD.Vmin + F_ext;
LD.Mmax = W_bench*x_W + W_beam*x_beam + F_ext*x_F + LD.Vmax*B.t_plate;
LD.Tmax = F_ext*L2;

fprintf('\nJoint loads  min (F = 0):  V = %.0f N, M = %.1f N.m, T = %.1f N.m\n', LD.Vmin, LD.Mmin, LD.Tmin);
fprintf('Joint loads  max (F = %.0f): V = %.0f N, M = %.1f N.m, T = %.1f N.m\n', F_ext, LD.Vmax, LD.Mmax, LD.Tmax);

%% NUMBER OF BOLTS x DIAMETER x CLASS
S = [];  store = {};
for nb = n_list
    for i = 1:size(bolt_sizes,1)
        for j = 1:size(classes,1)
            r = eval_bolted(nb, bolt_sizes(i,1), classes(j,:), bolt_sizes(i,2), bolt_sizes(i,3), B, LD);
            S(end+1,:) = [nb, bolt_sizes(i,1), classes(j,1), r.C, r.Fi/1e3, r.nf, r.nL, r.n0, r.np, r.ns, r.nbr, r.n_gov, r.ok, r.geom_ok, r.sl_ok]; %#ok<SAGROW>
            store{end+1} = r; 
        end
    end
end

%% TOP CANDIDATES
idxf = find(S(:,13) == 1);
if isempty(idxf)
    fprintf('\nNo bolt configuration met every check. Relax assumptions or enlarge search.\n'); return;
end
[~, o] = sortrows(S(idxf,[1 2 3]));
idxf = idxf(o);

fprintf('\nFatigue factor: %s criterion, %s definition\n', B.crit, B.fs_def);
fprintf('TOP 5 FEASIBLE JOINTS (ranked: fewest bolts, smallest diameter, lowest class)\n');
fprintf(' Rank | Bolts | Size | Class | C     | Preload (kN) | n_fatigue | n_load | n_sep | n_proof | n_shear | n_bearing\n');
for k = 1:min(5, numel(idxf))
    s = S(idxf(k),:);
    fprintf('  #%d  |   %2d  | M%-3d | %4.1f  | %.3f | %10.1f   | %8.1f  | %6.1f | %5.1f | %6.2f  | %6.2f  | %8.1f\n', ...
        k, s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8), s(9), s(10), s(11));
end

sel = idxf(1);  R = store{sel};  s = S(sel,:);
cls = classes(classes(:,1) == s(3), :);
Sp = cls(2); Sut = cls(3); Sy = cls(4); Se = cls(5);

%% SELECTED JOINT
fprintf('\n=== SELECTED JOINT: %d x M%d, ISO class %.1f ===\n', s(1), s(2), s(3));
fprintf('Stress area At = %.1f mm^2 | Proof load Fp = %.1f kN | Preload Fi = %.1f kN (%.0f%% of proof)\n', R.At, R.Fp/1e3, R.Fi/1e3, B.frac*100);

%% FATIGUE DIAGRAM
fig1 = figure('Name', 'Bolt fatigue diagram', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 10, 9]);
t1 = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_main = nexttile(t1); hold on; grid on;

baseFontSize = 10;
labelFontSize = 11;
legendFontSize = 9;

Sm_line = linspace(0, Sut, 600);
switch lower(B.crit)
    case 'goodman', Sa_line = Se*(1 - Sm_line/Sut);
    case 'gerber',  Sa_line = Se*(1 - (Sm_line/Sut).^2);
    case 'asme',    Sm_line = linspace(0, Sy, 600); Sa_line = Se*sqrt(1 - (Sm_line/Sy).^2);
end

plot(ax_main, Sm_line, Sa_line, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Criterion Line');
plot(ax_main, [0 Sp], [Sp 0], 'r--', 'LineWidth', 1.2, 'DisplayName', 'Proof Line');
plot(ax_main, [R.si R.Sm_c], [0 R.Sa_c], 'k-', 'LineWidth', 1.2, 'DisplayName', 'Load Line');
plot(ax_main, R.si, 0, 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 6.5, 'DisplayName', '$\sigma_i$ (Preload)');

nbolt  = numel(R.sm);
mk1    = {'o', 's', '^', 'd', 'v', 'p', 'h', '>'};
col1   = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
          0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; ...
          0.6350 0.0780 0.1840; 0.2500 0.2500 0.2500];

for b = 1:nbolt
    c  = col1(mod(b-1, size(col1,1)) + 1, :);
    m  = mk1{mod(b-1, numel(mk1)) + 1};
    ms = 11 - 2*(b-1);  ms = max(ms, 4);
    plot(ax_main, R.sm(b), R.sa(b), m, 'Color', 'k', 'MarkerFaceColor', c, 'LineWidth', 0.8, ...
        'MarkerSize', ms, 'LineStyle', 'none', 'DisplayName', sprintf('Bolt %d', b));
end
plot(ax_main, R.Sm_c, R.Sa_c, 'k^', 'MarkerFaceColor', 'r', 'MarkerSize', 6.5, 'DisplayName', 'Intersection');

xlabel(ax_main, 'Mean Stress $\sigma_m$ (MPa)', 'Interpreter', 'latex', 'FontSize', labelFontSize);
ylabel(ax_main, 'Alternating Stress $\sigma_a$ (MPa)', 'Interpreter', 'latex', 'FontSize', labelFontSize);

set(ax_main, 'TickLabelInterpreter', 'latex', 'Box', 'off', 'FontSize', baseFontSize);
xlim(ax_main, [0 1.05*Sut]);
ylim(ax_main, [0 1.05*Se]);

lgd1 = legend(ax_main, 'Location', 'northeast', 'NumColumns', 1, 'Interpreter', 'latex', 'FontSize', legendFontSize);
lgd1.Box = 'off';

drawnow;                                   
ip = ax_main.InnerPosition;                
xl = xlim(ax_main);  yl = ylim(ax_main);

inX    = [0.08, 0.44] * xl(2);            
inYbot = 0.08 * yl(2);
ySafe  = interp1(Sm_line, Sa_line, inX(2), 'linear', Inf);  
inYtop = min(0.38*yl(2), 0.80*ySafe);      

toFig = @(xd, yd) [ip(1) + ip(3)*(xd - xl(1))/diff(xl), ...
                   ip(2) + ip(4)*(yd - yl(1))/diff(yl)];
p0 = toFig(inX(1), inYbot);
p1 = toFig(inX(2), inYtop);
ax_inset = axes('Parent', fig1, 'Position', [p0, p1 - p0], 'Color', 'w');
hold(ax_inset, 'on'); box(ax_inset, 'on'); grid(ax_inset, 'on');

xs = R.sm(:);
ys = R.sa(:);

xr   = max(xs) - min(xs);
yr   = max(ys) - min(ys);
xpad = max(1.0,  0.75*xr);
ypad = max(0.25, 0.60*yr);
x0   = min(xs) - xpad;
clusterW = xr + 2*xpad;
xlim(ax_inset, [x0, x0 + clusterW/0.30]);
ylim(ax_inset, [max(0, min(ys) - ypad), max(ys) + ypad]);

xl_in = xlim(ax_inset);  yl_in = ylim(ax_inset);
xLab  = xl_in(1) + 0.36*diff(xl_in);       

[~, ord] = sort(ys, 'descend');
slotY    = linspace(yl_in(2) - 0.14*diff(yl_in), yl_in(1) + 0.14*diff(yl_in), nbolt);

for k = 1:nbolt
    b = ord(k);
    plot(ax_inset, [xs(b), xLab - 0.02*diff(xl_in)], [ys(b), slotY(k)], '-', ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
end

for b = 1:nbolt
    c  = col1(mod(b-1, size(col1,1)) + 1, :);
    m  = mk1{mod(b-1, numel(mk1)) + 1};
    ms = 9 - 1.5*(b-1);  ms = max(ms, 4);
    plot(ax_inset, xs(b), ys(b), m, 'Color', 'k', 'MarkerFaceColor', c, ...
        'LineWidth', 0.8, 'MarkerSize', ms, 'LineStyle', 'none');
end

for k = 1:nbolt
    b = ord(k);
    text(ax_inset, xLab, slotY(k), sprintf('B%d (%.1f, %.2f)', b, xs(b), ys(b)), ...
        'FontSize', 7, 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
end

set(ax_inset, 'TickLabelInterpreter', 'latex', 'FontSize', 7);
title(ax_inset, 'Bolt Zoom ($\sigma_m,\ \sigma_a$)', 'Interpreter', 'latex', 'FontSize', 8);

exportgraphics(fig1, 'PartB_fatigue_diagram.pdf', 'ContentType', 'vector');
exportgraphics(fig1, 'PartB_fatigue_diagram.png', 'Resolution', 300);

%% SAFETY FACTOR
fig2 = figure('Name', 'Governing safety factor', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 9, 7.5]);
t2 = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t2); hold on; grid on;

line_styles = {'-', '--', '-.', ':', '-'};
markers     = {'o', 's', '^', 'd', 'v'};
colors      = [ ...
    0.0000 0.4470 0.7410;  % Blue
    0.8500 0.3250 0.0980;  % Red/Orange
    0.9290 0.6940 0.1250;  % Yellow/Gold
    0.4940 0.1840 0.5560;  % Purple
    0.4660 0.6740 0.1880]; % Green

for k = 1:numel(n_list)
    nb = n_list(k);
    m = S(:,1) == nb & S(:,3) == s(3);
    
    ls = line_styles{mod(k-1, numel(line_styles)) + 1};
    mk = markers{mod(k-1, numel(markers)) + 1};
    clr = colors(mod(k-1, size(colors,1)) + 1, :);
    
    plot(S(m,2), S(m,12), ...
        'LineStyle', ls, ...
        'Marker', mk, ...
        'Color', clr, ...
        'LineWidth', 1.3, ...
        'MarkerSize', 5, ...
        'MarkerFaceColor', clr, ...
        'DisplayName', sprintf('%d bolts', nb));
end

yline(B.n_req, 'r--', sprintf('Req. $n = %.1f$', B.n_req), ...
    'Interpreter', 'latex', ...
    'LineWidth', 1.2, ...
    'HandleVisibility', 'off', ...
    'LabelHorizontalAlignment', 'left', ...  
    'LabelVerticalAlignment', 'bottom');

set(gca, 'YScale', 'log', 'TickLabelInterpreter', 'latex', 'Box', 'off');

ylim([1.45, 1.85]); 

xlabel('Bolt Diameter (mm)', 'Interpreter', 'latex'); 
ylabel('Governing Safety Factor', 'Interpreter', 'latex');
title(sprintf('Class %.1f Safety Factors', s(3)), 'Interpreter', 'latex');

lgd2 = legend('Location', 'east', 'NumColumns', 1, 'Interpreter', 'latex', 'FontSize', 7);
lgd2.Box = 'off';

% Export
exportgraphics(fig2, 'PartB_governing_safety_factor.pdf', 'ContentType', 'vector');
exportgraphics(fig2, 'PartB_governing_safety_factor.png', 'Resolution', 300);

%% 9. FATIGUE DIAGRAM BEST CONFIG PER BOLT COUNT
nN = numel(n_list);
best_idx = zeros(1, nN);  feas9 = false(1, nN);
for k = 1:nN
    ik  = find(S(:,1) == n_list(k));
    ikf = ik(S(ik,13) == 1);
    if ~isempty(ikf)
        [~, o] = sortrows(S(ikf, [2 3]));
        best_idx(k) = ikf(o(1));  feas9(k) = true;
    else
        [~, o] = max(S(ik,12));
        best_idx(k) = ik(o);
    end
end

fig3 = figure('Name', 'Fatigue diagram vs number of bolts', 'Color', 'w', ...
    'Units', 'centimeters', 'Position', [2, 2, 15, 8]);
t3 = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile(t3); hold on; grid on;

cols9 = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; ...
         0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880];
mks9  = {'o', 's', '^', 'd', 'v'};
Sut_max = 0;  Se_max = 0;

for k = 1:nN
    id  = best_idx(k);  rk = store{id};  w = rk.worst;
    cl  = classes(classes(:,1) == S(id,3), :);
    Sp9 = cl(2);  Sut9 = cl(3);  Sy9 = cl(4);  Se9 = cl(5);
    Sut_max = max(Sut_max, Sut9);  Se_max = max(Se_max, Se9);
    c = cols9(mod(k-1, size(cols9,1)) + 1, :);

    switch lower(B.crit)
        case 'goodman', Sm9 = linspace(0, Sut9, 300); Sa9 = Se9*(1 - Sm9/Sut9);
        case 'gerber',  Sm9 = linspace(0, Sut9, 300); Sa9 = Se9*(1 - (Sm9/Sut9).^2);
        case 'asme',    Sm9 = linspace(0, Sy9, 300);  Sa9 = Se9*sqrt(1 - (Sm9/Sy9).^2);
    end
    plot(Sm9, Sa9, '-', 'Color', [c 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');

    % Load path from preload to the worst-case bolt
    plot([rk.si rk.sm(w)], [0 rk.sa(w)], '-', 'Color', c, 'LineWidth', 0.8, 'HandleVisibility', 'off');

    tag = ''; if ~feas9(k), tag = '*'; end
    plot(rk.sm(w), rk.sa(w), mks9{k}, 'Color', 'k', 'MarkerFaceColor', c, 'MarkerSize', 6.5, ...
        'LineStyle', 'none', ...
        'DisplayName', sprintf('%d bolts: M%d, %.1f%s', n_list(k), S(id,2), S(id,3), tag));
end

xlabel('Mean Stress $\sigma_m$ (MPa)', 'Interpreter', 'latex');
ylabel('Alternating Stress $\sigma_a$ (MPa)', 'Interpreter', 'latex');
title(sprintf('Worst-case bolt, best configuration per bolt count (%s)', B.crit), 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'Box', 'off');
xlim([0 1.05*Sut_max]);
ylim([0 1.05*Se_max]);

lgd3 = legend('Location', 'eastoutside', 'NumColumns', 1, 'Interpreter', 'latex', 'FontSize', 7);
lgd3.Box = 'off';
if any(~feas9)
    lgd3.Title.String = '* = no feasible joint';
    lgd3.Title.FontSize = 6.5;
end

exportgraphics(fig3, 'PartB_fatigue_vs_bolts.pdf', 'ContentType', 'vector');
exportgraphics(fig3, 'PartB_fatigue_vs_bolts.png', 'Resolution', 300);

%% FORCE-DEFLECTION DIAGRAMS
nbolt = numel(R.sm);
fig4 = figure('Name', 'Force-Deflection Diagrams', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 16, 13]);
t4 = tiledlayout(ceil(nbolt/2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_list = gobjects(nbolt, 1); 
for b = 1:nbolt
    ax_list(b) = nexttile(t4); hold on; grid on;
    
    ax_list(b).FontSize = 12; 
    
    delta_b_i = R.Fi / R.kb;
    delta_m_i = R.Fi / R.km;
    
    % Forces at max external load
    Pb_max = R.Pmax(b);
    Fb_max = R.Fi + R.C * Pb_max;
    Fm_min = R.Fi - (1 - R.C) * Pb_max;
    
    % Max deflections
    delta_b_max = Fb_max / R.kb;
    delta_m_zero = delta_b_i + delta_m_i; 
    
    
    plot([0, delta_b_max]*1e3, [0, Fb_max]/1e3, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Bolt Line ($F_b$)');
    
    
    plot([delta_m_zero, delta_b_i, delta_b_max]*1e3, [0, R.Fi, Fm_min]/1e3, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Member Line ($F_m$)');
    
    
    plot(delta_b_i*1e3, R.Fi/1e3, 'ks', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'DisplayName', 'Preload ($F_i$)');
    
    
    xlabel('Deflection $\delta$ (mm)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Force (kN)', 'Interpreter', 'latex', 'FontSize', 14);
    title(sprintf('Bolt %d ($P_{max} = %.1f$ N)', b, Pb_max), 'Interpreter', 'latex', 'FontSize', 14);
end

lgd4 = legend(ax_list(1), 'Orientation', 'horizontal', 'Interpreter', 'latex', 'FontSize', 14);
lgd4.Layout.Tile = 'north'; 
lgd4.Box = 'on';            

exportgraphics(fig4, 'PartB_Force_Deflection_Bolts.pdf', 'ContentType', 'vector');
exportgraphics(fig4, 'PartB_Force_Deflection_Bolts.png', 'Resolution', 300);


%% FORCE VECTOR DIAGRAM 
fig5 = figure('Name', 'Bolt Force Vectors', 'Color', 'w', 'Units', 'centimeters', 'Position', [2, 2, 16, 13]);
t5 = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax5 = nexttile(t5); hold(ax5, 'on'); grid(ax5, 'on'); axis(ax5, 'equal');

ax5.FontSize = 12;

Y_mm = R.Y * 1e3;
Z_mm = R.Z * 1e3;

V_prim_z = LD.Vmax / nbolt; 
zc = R.Z - mean(R.Z); 
yc = R.Y - mean(R.Y);
Jp = sum(yc.^2 + zc.^2);

s_vec = 0.012; 

txt_dy = [-30, -30,  30,  30]; 
txt_dz = [ 30,  30, -30, -30];

for b = 1:nbolt
    text(ax5, Y_mm(b) + txt_dy(b), Z_mm(b) + txt_dz(b), ...
        sprintf('\\textbf{B%d}', b), 'FontSize', 15, 'Interpreter', 'latex', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'Margin', 1.0);
end

% Plot Shear Vectors
for b = 1:nbolt
    % Secondary Shear (Torsional)
    V_sec_y = -LD.Tmax * zc(b) / Jp;
    V_sec_z =  LD.Tmax * yc(b) / Jp;
    
    % Plot Primary Shear Arrow
    q1 = quiver(ax5, Y_mm(b), Z_mm(b), 0, V_prim_z * s_vec, 0, 'b', 'LineWidth', 2.5, 'MaxHeadSize', 0.9);
    
    % Plot Secondary Shear Arrow 
    q2 = quiver(ax5, Y_mm(b), Z_mm(b), V_sec_y * s_vec, V_sec_z * s_vec, 0, 'r', 'LineWidth', 2.5, 'MaxHeadSize', 0.9);
    
    % Plot Resultant Shear Arrow 
    q3 = quiver(ax5, Y_mm(b), Z_mm(b), V_sec_y * s_vec, (V_prim_z + V_sec_z) * s_vec, 0, 'k', 'LineWidth', 2.8, 'MaxHeadSize', 0.9);
    
    if b == 1
        q1.DisplayName = 'Primary Shear ($V/n$)';
        q2.DisplayName = 'Secondary Shear ($V_{T}$)';
        q3.DisplayName = 'Resultant Shear';
    else
        q1.HandleVisibility = 'off';
        q2.HandleVisibility = 'off';
        q3.HandleVisibility = 'off';
    end
end

xlabel(ax5, 'Y Position (mm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel(ax5, 'Z Position (mm)', 'Interpreter', 'latex', 'FontSize', 14);

xlim(ax5, [-190, 290]);
ylim(ax5, [-180, 180]); 

lgd5 = legend(ax5, 'Location', 'east', 'Orientation', 'vertical', 'Interpreter', 'latex', 'FontSize', 10);
lgd5.Box = 'on'; 

exportgraphics(fig5, 'PartB_Force_Vectors.pdf', 'ContentType', 'vector');
exportgraphics(fig5, 'PartB_Force_Vectors.png', 'Resolution', 300);

%% FUNCTIONS
function r = eval_bolted(nb, d_mm, cls, At_mm2, pitch_mm, B, LD)
    d = d_mm/1e3;  At = At_mm2;
    Sp = cls(2); Sut = cls(3); Sy = cls(4); Se = cls(5);
    rows  = nb/2;
    y_top = B.h/2 - B.edge;  x_b = B.w/2 - B.edge;
    Zs = linspace(-y_top, y_top, rows)';
    Y  = [-x_b*ones(rows,1); x_b*ones(rows,1)];
    Z  = [Zs; Zs];
    r.Y = Y;  r.Z = Z;
    r.geom_ok = (2*y_top/(rows-1) >= 3*d) && (2*x_b >= 3*d) && (B.edge >= 1.5*d);
    
    if strcmpi(B.joint, 'tapped')
        l_g    = B.t_plate + d/2;
        L_req  = B.t_plate + B.washer + B.engage*d + 2*pitch_mm*1e-3;
        L_bolt = ceil(L_req/0.005 - 1e-9)*0.005;
        [C, kb, km] = joint_C(d, At*1e-6, l_g, B.E, L_bolt);  lw = l_g;
    else
        C = 0;
        for tf = B.frames
            [Ci, kbi, kmi] = joint_C(d, At*1e-6, B.t_plate + tf, B.E, B.t_plate + tf + 0.8*d);
            if Ci > C, C = Ci; kb = kbi; km = kmi; lw = B.t_plate + tf; end
        end
        L_req  = B.t_plate + max(B.frames) + 2*B.washer + 0.8*d + 3*pitch_mm*1e-3;
        L_bolt = ceil(L_req/0.005 - 1e-9)*0.005;
    end
    len_ok = L_bolt <= B.L_stock + 1e-9;
    kb_ap = (pi/4*d^2)*B.E/lw;  C_ap = kb_ap/(kb_ap + km);
    Fp = Sp*At;  Fi = B.frac*Fp;  si = Fi/At;
    
    [Pmax, Vmax] = bolt_group_loads(Y, Z, LD.Vmax, LD.Mmax, LD.Tmax, B.h);
    [Pmin, ~]    = bolt_group_loads(Y, Z, LD.Vmin, LD.Mmin, LD.Tmin, B.h);
    
    Fb_max = Fi + C*Pmax;  Fb_min = Fi + C*Pmin;
    sig_max = Fb_max/At;   sig_min = Fb_min/At;
    sa = (sig_max - sig_min)/2;  sm = (sig_max + sig_min)/2;
    nf_ll = zeros(nb,1); nf_lec = nf_ll; Sac = nf_ll; Smc = nf_ll;
    for k = 1:nb
        [nf_ll(k), Sac(k), Smc(k)] = fatigue_nf(B.crit, sa(k), sm(k), si, Se, Sut, Sy);
        nf_lec(k) = fatigue_FS_lecture(B.crit, sa(k), sm(k), Se, Sut, Sy);
    end
    if strcmpi(B.fs_def, 'lecture'), nf_all = nf_lec; else, nf_all = nf_ll; end
    nL  = (Fp - Fi) ./ (C*max(Pmax, 1e-9));
    n0  = Fi ./ ((1 - C)*max(Pmax, 1e-9));
    np  = Fp ./ Fb_max;
    ns  = 0.577*Sy ./ (Vmax/At);
    nbr = B.Sy_plate ./ (Vmax/(d*B.t_plate)/1e6);
    nf_min = min(nf_all);
    cand = find(nf_all <= nf_min*(1 + 1e-9));
    [~, ii] = max(Vmax(cand));  kw = cand(ii);
    r.worst = kw;  r.nf = nf_all(kw);  r.Sa_c = Sac(kw);  r.Sm_c = Smc(kw);
    r.nf_ll = nf_ll(kw);  r.nf_lec = nf_lec(kw);  r.FS_pt = nf_lec(kw)*[sm(kw), sa(kw)];
    r.nL = min(nL);  r.n0 = min(n0);  r.np = min(np);  r.ns = min(ns);  r.nbr = min(nbr);
    r.n_gov = min([r.nf r.nL r.n0 r.ns r.nbr]);
    r.pitch_d = pitch_mm/d_mm;
    r.sl_ok   = r.pitch_d <= pi*B.mu_thread;
    r.L_bolt = L_bolt;  r.len_ok = len_ok;
    r.ok = r.geom_ok && r.sl_ok && r.len_ok && all([r.nf r.nL r.n0 r.ns r.nbr] >= B.n_req) && r.np >= 1.0;
    r.C = C; r.C_ap = C_ap; r.kb = kb; r.km = km; r.Fi = Fi; r.Fp = Fp; r.si = si; r.At = At;
    r.Pmax = Pmax; r.Pmin = Pmin; r.Vmax = Vmax; r.sa = sa; r.sm = sm; r.nf_all = nf_all;
end

function [P, Vs] = bolt_group_loads(Y, Z, Vz, My_j, Tx, h)
    n = numel(Y);
    d_piv = Z + h/2;
    P_piv = My_j*d_piv/sum(d_piv.^2);
    zc = Z - mean(Z);  yc = Y - mean(Y);
    P_cen = My_j*zc/sum(zc.^2);
    P = max([P_piv, P_cen, zeros(n,1)], [], 2);
    Jp = sum(yc.^2 + zc.^2);
    Fy = -Tx*zc/Jp;
    Fz =  Vz/n + Tx*yc/Jp;
    Vs = sqrt(Fy.^2 + Fz.^2);
end

function [C, kb, km] = joint_C(d, At, l, E, L)
    Ad = pi/4*d^2;
    if L <= 0.125, LT = 2*d + 0.006; elseif L <= 0.200, LT = 2*d + 0.012; else, LT = 2*d + 0.025; end
    ld = min(max(L - LT, 0), l);
    lt = l - ld;
    kb = At*Ad*E/(Ad*lt + At*ld);
    km = 0.5774*pi*E*d/(2*log(5*(0.5774*l + 0.5*d)/(0.5774*l + 2.5*d)));
    C  = kb/(kb + km);
end

function [nf, Sa_c, Sm_c] = fatigue_nf(crit, sa, sm, si, Se, Sut, Sy)
    if sa <= 0, nf = Inf; Sa_c = NaN; Sm_c = NaN; return; end
    k = (sm - si)/sa;
    switch lower(crit)
        case 'goodman'
            Sa_c = (1 - si/Sut)/(1/Se + k/Sut);
        case 'gerber'
            a2 = k^2/Sut^2;  a1 = 2*si*k/Sut^2 + 1/Se;  a0 = si^2/Sut^2 - 1;
            Sa_c = solve_pos(a2, a1, a0);
        case 'asme'
            a2 = 1/Se^2 + k^2/Sy^2;  a1 = 2*si*k/Sy^2;  a0 = si^2/Sy^2 - 1;
            Sa_c = solve_pos(a2, a1, a0);
        otherwise
            error('Unknown fatigue criterion: %s', crit);
    end
    Sm_c = si + k*Sa_c;
    nf   = Sa_c/sa;
end

function x = solve_pos(a2, a1, a0)
    if abs(a2) < 1e-14, x = -a0/a1; else, x = (-a1 + sqrt(a1^2 - 4*a2*a0))/(2*a2); end
end

function FS = fatigue_FS_lecture(crit, sa, sm, Se, Sut, Sy)
    if sa <= 0 && sm <= 0, FS = Inf; return; end
    switch lower(crit)
        case 'goodman'
            FS = 1/(sa/Se + sm/Sut);
        case 'gerber'
            a = sa/Se;  b = (sm/Sut)^2;
            if b < 1e-14, FS = 1/a; else, FS = (-a + sqrt(a^2 + 4*b))/(2*b); end
        case 'asme'
            FS = 1/sqrt((sa/Se)^2 + (sm/Sy)^2);
        otherwise
            error('Unknown fatigue criterion: %s', crit);
    end
end