clear; clc; close all;

%% MATERIAL AND REQUIREMENTS
P.g = 9.81;  P.rho = 7800;
P.E = 210e9; P.G = 79e9;
P.Sy  = 220.594e6;            
P.Sut = 399.826e6;            
n_req = 1.5;                  
require_fatigue = true;       
n_top = 5;                    
P.use_wb   = false;         
P.thin_lim = 10;            

%% LOAD POSITION 
sid_d1 = 1;   % last digit of SID
sid_d2 = 4;   % 2nd to last digit
sid_d4 = 0;   % 4th to last digit 
tbl = [  0     0     1      100 730;
         1     0     1      100 780;
         0     1     0      100 830;
         0     1     1      150 730;
         1     1     0      150 780;
         1     0     0      150 830;
         1     1     1      200 730;
         0     0     0      200 780 ];
row = find(all(tbl(:,1:3) == [mod(sid_d1,2) mod(sid_d2,2) mod(sid_d4,2)], 2));
P.L1 = tbl(row,4)/1000;   P.L2 = tbl(row,5)/1000;

%% GEOMETRY AND LOADS
P.L  = 0.500;                          % beam length (m)
P.F  = 3500;                           % alternating point load 0-3500 N
P.xF = P.L + P.L1;                     % F from Plate 1 face (m) 
P.W  = 0.06628 * P.rho * P.g;          % bench weight (N)
% Bench COM
a_offset = 0.120;  t_b = 0.035;  H_leg = 0.650;  L_leg = 0.650;
A_v = t_b*H_leg;           xc_v = t_b/2;
A_h = t_b*(L_leg - t_b);   xc_h = t_b + (L_leg - t_b)/2;
P.xW = a_offset + (A_v*xc_v + A_h*xc_h)/(A_v + A_h);   % 0.2955 m

P.T    = P.F * P.L2;                                    % torque (N.m)
P.Ptip = P.W + P.F;                                     % tip force (N)
P.Mtip = P.F*(P.xF - P.L) + P.W*(P.xW - P.L);           % tip couple (N.m); W term negative

%% LIMITS
Delta_ratio = 1/125;                   % cantilever, total deflection (span beams: 1/250)
l_ref       = P.xF;                    % reference cantilever length l (m). 0.650 = support to line of action of F.

Delta_lim_mm = Delta_ratio * l_ref * 1000;

use_twist_limit = false;
tilt_N          = 250;                 
if use_twist_limit
    theta_lim_deg = rad2deg(atan(1/tilt_N));
else
    theta_lim_deg = Inf;
end
fprintf('Loads: W = %.1f N, F = %.0f N, xW = %.4f m, xF = %.3f m, L2 = %.3f m\n', P.W, P.F, P.xW, P.xF, P.L2);
fprintf('Root actions: V = %.0f N, M = %.1f N.m, T = %.1f N.m\n', P.W+P.F, P.W*P.xW+P.F*P.xF, P.T);
fprintf('Limits: n >= %.1f | total deflection at F <= 1/%g x %.0f mm = %.2f mm (AS 4100 Table B.1)\n', ...
    n_req, 1/Delta_ratio, l_ref*1000, Delta_lim_mm);
if use_twist_limit
    fprintf('        extra twist limit 1 in %g = %.3f deg (drop at F = %.2f mm)\n', tilt_N, theta_lim_deg, tan(deg2rad(theta_lim_deg))*P.L2*1000);
else
    fprintf('        no separate twist limit (twist enters through the total deflection)\n');
end
fprintf('\n');

%% CANDIDATES 
fam_names = {'RHS/SHS','CHS','Solid rect','I-beam'};

cands = [];
for b = 50:10:100, for h = 50:10:100, for t = 3:0.5:12, cands(end+1,:) = [1 b h t 0]; end, end, end
for D = 50:10:100, for t = 3:0.5:12, cands(end+1,:) = [2 D D t 0]; end, end
for b = 50:10:100, for h = 50:10:100, cands(end+1,:) = [3 b h 0 0]; end, end
for b = 50:10:100, for h = 50:10:100, for tf = 3:1:12, for tw = 3:1:12, cands(end+1,:) = [4 b h tf tw]; end, end, end, end

Rall = [];
for i = 1:size(cands,1)
    c = cands(i,:);
    r = eval_section(c(1), c(2)/1000, c(3)/1000, c(4)/1000, c(5)/1000, P);
    if r.valid, Rall(end+1,:) = result_row(c, r); end 
end
Rr = Rall(Rall(:,1)==1, :);              
hdr  = sprintf(' Label      | Family     | b x h (mm) | t (mm)    | Mass (kg) | vM (MPa) | n_static | n_fatigue | Total@F (mm)\n');
hdr2 = sprintf(' Label      | Beam tip defl (mm) | Beam twist (deg) | Bending@F (mm) | Twist drop@F (mm) | Total@F (mm) | Limit (mm)\n');

%% SOLVE
idx_top = pick_top(Rr, n_req, require_fatigue, Delta_lim_mm, theta_lim_deg, n_top);
if isempty(idx_top)
    fprintf('No RHS/SHS candidate meets all constraints. Relax limits or enlarge the search.\n'); return;
end
best_theory = Rr(idx_top(1),:);
Top_theory  = Rr(idx_top,:);
fprintf('TOP %d LIGHTEST FEASIBLE SECTIONS (theoretical, 0.5 mm wall steps)\n%s', numel(idx_top), hdr);
for k = 1:numel(idx_top), print_row(sprintf('#%d', k), Top_theory(k,:), fam_names); end
fprintf('%s', hdr2);
for k = 1:numel(idx_top), print_detail(sprintf('#%d', k), Top_theory(k,:), Delta_lim_mm); end

%% STANDARD CATALOGUE SIZES
std_list = [ 20  20 1.6; 20  20 2.0; ...
             25  25 1.6; 25  25 2.0; 25  25 3.0; ...
             32  32 1.6; 32  32 2.0; 32  32 3.0; ...
             38  38 2.0; 38  38 3.0; ...
             40  40 2.0; 40  40 3.0; 40  40 4.0; ...
             50  50 2.0; 50  50 3.0; 50  50 4.0; 50  50 5.0; ...
             65  65 3.0; 65  65 4.0; 65  65 5.0; ...
             75  75 3.0; 75  75 4.0; 75  75 5.0; ...
             90  90 4.0; 90  90 5.0; ...
            100 100 3.0; 100 100 4.0; 100 100 5.0; 100 100 6.0 ];
Rs = [];
for i = 1:size(std_list,1)
    c = [1 std_list(i,:) 0];
    r = eval_section(1, c(2)/1000, c(3)/1000, c(4)/1000, 0, P);
    if r.valid, Rs(end+1,:) = result_row(c, r); end 
end
idx_s = pick_top(Rs, n_req, require_fatigue, Delta_lim_mm, theta_lim_deg, n_top);
if ~isempty(idx_s)
    best = Rs(idx_s(1),:);
    Top_std = Rs(idx_s,:);
    fprintf('\nTOP %d LIGHTEST FEASIBLE SECTIONS (standard sizes)\n%s', numel(idx_s), hdr);
    for k = 1:numel(idx_s), print_row(sprintf('#%d', k), Top_std(k,:), fam_names); end
    fprintf('%s', hdr2);
    for k = 1:numel(idx_s), print_detail(sprintf('#%d', k), Top_std(k,:), Delta_lim_mm); end
else
    best = best_theory;
    fprintf('\nNo listed standard size meets all constraints; using the theoretical optimum.\n');
end

%% BEST SECTION FROM EACH FAMILY
nfam = numel(fam_names);
Rbest = nan(nfam, size(Rall,2));
infeas = false(nfam,1);
for f = 1:nfam
    Rf = Rall(Rall(:,1)==f, :);
    i = pick_min_mass(Rf, n_req, require_fatigue, Delta_lim_mm, theta_lim_deg);
    if ~isempty(i)
        Rbest(f,:) = Rf(i,:);
    else
        % No section meets every constraint: show the strength-feasible section with the lowest total deflection
        ok2 = Rf(:,7) >= n_req;
        if require_fatigue, ok2 = ok2 & (Rf(:,8) >= n_req); end
        c2 = find(ok2);
        if ~isempty(c2)
            [~, k2] = min(Rf(c2,15));
            Rbest(f,:) = Rf(c2(k2),:);
            infeas(f) = true;
        end
    end
end
m_ref = Rbest(1,5);     
fprintf('\nBEST SECTION PER FAMILY (lightest that meets all constraints)\n');
fprintf(' Family     | Section                  | Mass (kg) | vs RHS | vM (MPa) | n_stat | n_fat | Tip defl (mm) | Twist (deg) | Bend@F (mm) | Drop@F (mm) | Total@F (mm) | EI/m (N.m2/kg) | GJ/m (N.m2/kg)\n');
for f = 1:nfam
    r = Rbest(f,:);
    nm = fam_names{f};  if infeas(f), nm = [nm '*']; end
    if isnan(r(1))
        fprintf(' %-10s | no section passes even the strength checks\n', nm);
    else
        fprintf(' %-10s | %-24s | %8.2f  | %+5.0f%% | %7.1f  | %5.2f  | %5.2f | %13.3f | %11.4f | %11.3f | %11.3f | %12.3f | %14.3g | %14.3g\n', ...
            nm, dim_str(r), r(5), (r(5)/m_ref - 1)*100, r(6), r(7), r(8), r(14), r(10), r(9), r(11), r(15), r(12), r(13));
    end
end
fprintf('  (CHS: OD x wall; solid: width x height; I-beam: width x depth, flange tf, web tw. Thin-wall rule t <= min(b,h)/%g.)\n', P.thin_lim);
if any(infeas)
    fprintf('  * No section in this family meets every constraint. The row shows the strength-feasible section with the\n');
    fprintf('    lowest total deflection at F (compare its Total@F with the %.2f mm limit).\n', Delta_lim_mm);
end
fprintf('  I-beam torsion uses St Venant theory only; warping restraint at Plate 1 is neglected (see notes).\n');

fprintf('\nCONSTRAINT RATIOS FOR SELECTED SECTION (%g x %g x %g mm)\n', best(2), best(3), best(4));
fprintf('  Static strength : n_req/n   = %.2f/%.2f = %.2f\n', n_req, best(7), n_req/best(7));
if require_fatigue
    fprintf('  Fatigue         : n_req/nf  = %.2f/%.2f = %.2f\n', n_req, best(8), n_req/best(8));
end
fprintf('  Total deflection: %.2f / %.2f mm = %.2f   (AS 4100)\n', best(15), Delta_lim_mm, best(15)/Delta_lim_mm);
if use_twist_limit
    fprintf('  Twist angle     : %.3f / %.3f deg = %.2f\n', best(10), theta_lim_deg, best(10)/theta_lim_deg);
end
fprintf('  (ratio closest to 1.00 governs)\n');

fprintf('\nSENSITIVITY TO THE DEFLECTION LIMIT (strength and fatigue constraints unchanged)\n');
fprintf('  Delta/l | l (mm) | limit (mm) | theoretical optimum              | standard-size optimum\n');
for ratio = [1/250 1/125]
    for lref = [P.L, P.xF]
        lim = ratio*lref*1000;
        i1 = pick_min_mass(Rr, n_req, require_fatigue, lim, theta_lim_deg);
        i2 = pick_min_mass(Rs, n_req, require_fatigue, lim, theta_lim_deg);
        fprintf('  1/%-5g | %6.0f | %9.2f  | %-32s | %s\n', 1/ratio, lref*1000, lim, sect_str(Rr, i1), sect_str(Rs, i2));
    end
end
i3 = pick_min_mass(Rr, n_req, false, Delta_lim_mm, theta_lim_deg);
fprintf('  Without the fatigue constraint (current limit): %s\n', sect_str(Rr, i3));
fprintf('  -> if this differs from the optimum above, fatigue (and its assumed Marin factors) is what governs.\n');

%% STIFFNESS PER UNIT MASS BY FAMILY 
fprintf('\nBEST STIFFNESS-TO-MASS BY FAMILY (all sizes, no strength filter)\n');
fprintf('  Family     | max EI/m (N.m^2/kg)  at              | max GJ/m (N.m^2/kg)  at\n');
for f = 1:numel(fam_names)
    Rf = Rall(Rall(:,1)==f,:);
    [e, ie] = max(Rf(:,12));  [j, ij] = max(Rf(:,13));
    fprintf('  %-10s | %10.3g  %-24s | %10.3g  %s\n', fam_names{f}, e, dim_str(Rf(ie,:)), j, dim_str(Rf(ij,:)));
end

%% HAND-CHECK VALUES
r = eval_section(best(1), best(2)/1000, best(3)/1000, best(4)/1000, best(16)/1000, P);
fprintf('\nHAND CHECK for selected section %g x %g x %g mm:\n', best(2), best(3), best(4));
fprintf('  Section: A = %.1f mm^2, Iy = %.4g mm^4, J = %.4g mm^4\n', r.A*1e6, r.I*1e12, r.J*1e12);
fprintf('  Stress : sigma = %.2f MPa, tau_T = %.2f MPa, tau_V = %.2f MPa, vM = %.2f MPa, n = %.2f\n', ...
    r.sig/1e6, r.tauT/1e6, r.tauV/1e6, r.vm/1e6, r.n);
fprintf('  Defl.  : tip = %.3f mm, tip slope = %.4g rad, bending@F = %.3f mm\n', r.d_tip*1e3, r.phi, r.dF*1e3);
fprintf('  Twist  : theta = %.4g rad (%.4f deg), drop@F = %.3f mm, total@F = %.3f mm\n', r.theta_rad, r.theta_deg, r.drop*1e3, r.total*1e3);

%% FOR 100 x 100 SHS
tt = 3:0.5:10;  [m_t, b_t, d_t, tot_t, n_t, nf_t] = deal(nan(size(tt)));
for k = 1:numel(tt)
    r = eval_section(1, 0.100, 0.100, tt(k)/1000, 0, P);
    if r.valid
        m_t(k) = r.mass; b_t(k) = r.dF*1e3; d_t(k) = r.drop*1e3; tot_t(k) = r.total*1e3; n_t(k) = r.n; nf_t(k) = r.nf;
    end
end

font_name = 'Arial';
font_size = 12;

% Figure 1: Mass Sensitivity
fig1 = figure('Name', 'Mass Sensitivity', 'Units', 'centimeters', 'Position', [2 2 12 9]);
plot(tt, m_t, 'o-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'MarkerFaceColor', [0 0.4470 0.7410]);
ylabel('Mass (kg)'); xlabel('Wall Thickness t (mm)'); grid on;
title('100 x 100 mm SHS: Mass Sensitivity', 'FontSize', font_size, 'FontWeight', 'bold');
set(gca, 'FontName', font_name, 'FontSize', font_size);
exportgraphics(fig1, 'PartA_mass_sensitivity.pdf', 'ContentType', 'vector');

% Figure 2: Deflection Trade-off
fig2 = figure('Name', 'Deflection Trade-off', 'Units', 'centimeters', 'Position', [2 2 12 9]);
p1 = plot(tt, b_t, 's--', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'MarkerFaceColor', [0.8500 0.3250 0.0980]); hold on;
p2 = plot(tt, d_t, '^:', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5, 'MarkerFaceColor', [0.9290 0.6940 0.1250]);
p3 = plot(tt, tot_t, 'o-', 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 1.5, 'MarkerFaceColor', [0.4940 0.1840 0.5560]);
p_lim2 = yline(Delta_lim_mm, 'r--', sprintf('AS 4100 Limit (%.2f mm)', Delta_lim_mm), 'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom');
ylabel('Deflection at F (mm)'); xlabel('Wall Thickness t (mm)'); grid on;
ylim([0, max(tot_t)*1.15]);
title('100 x 100 mm SHS: Deflections', 'FontSize', font_size, 'FontWeight', 'bold');

% Correct entries
lgd2 = legend([p1, p2, p3], {'Bending', 'Twist Drop', 'Total Deflection'}, ...
    'FontName', font_name, ...
    'FontSize', 10, ...
    'Location', 'northeast', ...
    'Orientation', 'vertical');
lgd2.Color = [1 1 1 0.85];
lgd2.Box = 'on';
set(gca, 'FontName', font_name, 'FontSize', font_size);
exportgraphics(fig2, 'PartA_deflection_tradeoff.pdf', 'ContentType', 'vector');

% Safety Factors
fig3 = figure('Name', 'Safety Factors', 'Units', 'centimeters', 'Position', [2 2 12 9]);
plot(tt, n_t, 's-', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5, 'MarkerFaceColor', [0.4660 0.6740 0.1880]); hold on;
plot(tt, nf_t, 'o-', 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 1.5, 'MarkerFaceColor', [0.3010 0.7450 0.9330]);
yline(n_req, 'r--', sprintf('n_{req} = %.1f', n_req), 'LineWidth', 1.2);
ylabel('Safety Factor'); xlabel('Wall Thickness t (mm)'); grid on;
title('100 x 100 mm SHS: Safety Factors', 'FontSize', font_size, 'FontWeight', 'bold');

lgd3 = legend({'Static ($n$)', 'Fatigue ($n_f$)'}, ...
    'Interpreter', 'latex', ...
    'FontName', font_name, ...
    'FontSize', font_size, ...
    'Location', 'northwest', ...
    'Orientation', 'vertical');
lgd3.Color = [1 1 1 0.85];
lgd3.Box = 'on';
set(gca, 'FontName', font_name, 'FontSize', font_size);
exportgraphics(fig3, 'PartA_safety_factors.pdf', 'ContentType', 'vector');

%% WALL THICKNESS vs TOTAL DEFLECTION AT F
ok = Rr(:,7) >= n_req;
if require_fatigue, ok = ok & (Rr(:,8) >= n_req); end
fig4 = figure('Name', 'Thickness vs Deflection', 'Units', 'centimeters', 'Position', [2 2 12 9]);
sc = scatter(Rr(ok,4), Rr(ok,15), 35, Rr(ok,10), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5); hold on;

cb = colorbar;
cb.Label.String = 'Twist Angle \theta (deg)';
colormap(turbo); grid on;

p_lim = yline(Delta_lim_mm, 'r--', 'LineWidth', 1.2);
p_opt = plot(best_theory(4), best_theory(15), 'kp', ...
    'MarkerSize', 12, 'MarkerFaceColor', 'y', 'LineWidth', 1.5);
xlabel('Wall Thickness t (mm)'); ylabel('Total Deflection at F (mm)');
title('Feasible RHS/SHS Sections: Deflection vs. Thickness', 'FontSize', font_size, 'FontWeight', 'bold');

lgd4 = legend([p_opt, p_lim], ...
    {'Theoretical Optimum', sprintf('AS 4100 Limit (%.2f mm)', Delta_lim_mm)}, ...
    'FontName', font_name, ...
    'FontSize', 9, ...
    'Location', 'southwest', ...
    'Orientation', 'vertical');
lgd4.Color = [1 1 1 0.85];
lgd4.Box = 'on';

xlim([min(Rr(ok,4))-0.5, max(Rr(ok,4))+0.5]);
ylim([min(Rr(ok,15))-0.5, max(Rr(ok,15))+0.8]);
set(gca, 'FontName', font_name, 'FontSize', font_size);
exportgraphics(fig4, 'PartA_thickness_vs_deflection.pdf', 'ContentType', 'vector');
exportgraphics(fig4, 'PartA_thickness_vs_deflection.png', 'Resolution', 300);

%% SAVE 
writetable(array2table(Rr, 'VariableNames', {'fam','b_mm','h_mm','t_mm','mass_kg','vM_MPa','n_static','n_fatigue','bend_at_F_mm','twist_deg','twist_drop_mm','EI_per_m','GJ_per_m','tip_defl_mm','total_at_F_mm','tw_mm'}), 'PartA_ranked_sections.csv');
save('PartA_selected_section.mat', 'best', 'best_theory', 'Rbest', 'P', 'n_req', 'Delta_lim_mm', 'Delta_ratio', 'l_ref');
fprintf('\nSelected: %s %g x %g x %g mm saved to PartA_selected_section.mat\n', fam_names{best(1)}, best(2), best(3), best(4));

function r = eval_section(fam, b, h, t, tw, P)
    r.valid = false;
    s = section_props(fam, b, h, t, tw, P.T, P.thin_lim);
    if ~s.valid, return; end
    m  = s.A * P.L * P.rho;
    Wb = P.use_wb * m * P.g;                        
    V  = P.W + P.F + Wb;
    M  = P.W*P.xW + P.F*P.xF + Wb*P.L/2;
    sig    = M * s.c / s.I;
    tauV   = V / s.As;
    sig_vm = sqrt(sig^2 + 3*(s.tauT + tauV)^2);      
    EI = P.E * s.I;
    d_tip = P.Ptip*P.L^3/(3*EI) + P.Mtip*P.L^2/(2*EI) + Wb*P.L^3/(8*EI);   
    phi   = P.Ptip*P.L^2/(2*EI) + P.Mtip*P.L/EI      + Wb*P.L^2/(6*EI);   
    theta = P.T * P.L / (P.G * s.J);                                        
    r.d_tip = d_tip;  r.phi = phi;
    r.dF    = d_tip + phi*P.L1;                     
    r.theta_rad = theta;  r.theta_deg = rad2deg(theta);
    r.drop  = theta * P.L2;                          
    r.total = r.dF + r.drop;                         
    r.mass = m;  r.vm = sig_vm;  r.n = P.Sy / sig_vm;
    r.sig = sig; r.tauT = s.tauT; r.tauV = tauV;
    r.nf = fatigue_goodman(fam, b*1000, h*1000, sig_vm, P.Sut);
    r.EIm = EI / m;  r.GJm = P.G * s.J / m;
    r.A = s.A; r.I = s.I; r.J = s.J;
    r.valid = true;
end

function s = section_props(fam, b, h, t, tw, T, thin_lim)
    s.valid = false; s.A = 0; s.I = 1; s.c = 0; s.J = 1; s.tauT = 0; s.As = 1;
    switch fam
        case 1   
            if (b-2*t <= 0) || (h-2*t <= 0) || (t > min(b,h)/thin_lim + 1e-9), return; end
            s.A  = b*h - (b-2*t)*(h-2*t);
            s.I  = (b*h^3 - (b-2*t)*(h-2*t)^3)/12;
            s.c  = h/2;
            Am = (b-t)*(h-t);  pm = 2*((b-t)+(h-t));
            s.J  = 4*Am^2*t/pm;
            s.tauT = T/(2*Am*t);
            s.As = 2*(h-t)*t;                        
        case 2   % circular hollow section
            d = b - 2*t;
            if (d <= 0) || (t > b/thin_lim + 1e-9), return; end
            s.A  = pi/4*(b^2 - d^2);
            s.I  = pi/64*(b^4 - d^4);
            s.c  = b/2;
            s.J  = pi/32*(b^4 - d^4);
            s.tauT = T*(b/2)/s.J;
            s.As = s.A/2;
        case 3   % solid rectangle 
            s.A  = b*h;
            s.I  = b*h^3/12;
            s.c  = h/2;
            a = max(b,h)/2;  bb = min(b,h)/2;
            s.J  = a*bb^3*(16/3 - 3.36*(bb/a)*(1 - bb^4/(12*a^4)));
            s.tauT = T*(3*a + 1.8*bb)/(8*a^2*bb^2);
            s.As = s.A/1.5;
        case 4   
            hw = h - 2*t;
            if (hw <= 0) || (tw <= 0) || (tw >= b) || (max(t,tw) > min(b,h)/thin_lim + 1e-9), return; end
            s.A  = 2*b*t + hw*tw;
            s.I  = (b*h^3 - (b-tw)*hw^3)/12;
            s.c  = h/2;
            s.J  = (2*b*t^3 + hw*tw^3)/3;            
            s.tauT = T*max(t,tw)/s.J;                
            s.As = hw*tw;                            
    end
    s.valid = true;
end

function nf = fatigue_goodman(fam, b_mm, h_mm, sig_vm, Sut)
    Su = Sut/1e6;  Se_p = 0.5*Su;
    ka = 57.7 * Su^(-0.718);
    if fam == 2, de = 0.370*b_mm; else, de = 0.808*sqrt(b_mm*h_mm); end
    if de <= 51, kb = 1.24*de^(-0.107); else, kb = 1.51*de^(-0.157); end
    ke = 0.814;
    Se = ka*kb*ke*Se_p;
    sa = sig_vm/1e6/2;  sm = sa;
    nf = 1/(sa/Se + sm/Su);
end

function row = result_row(c, r)
    row = [c(1) c(2) c(3) c(4) r.mass r.vm/1e6 r.n r.nf r.dF*1e3 r.theta_deg r.drop*1e3 r.EIm r.GJm r.d_tip*1e3 r.total*1e3 c(5)];
end

function idx = pick_top(R, n_req, req_fat, Delta_lim_mm, th_lim_deg, N)
    if isempty(R), idx = []; return; end
    ok = (R(:,7) >= n_req) & (R(:,15) <= Delta_lim_mm) & (R(:,10) <= th_lim_deg);
    if req_fat, ok = ok & (R(:,8) >= n_req); end
    cand = find(ok);
    if isempty(cand), idx = []; return; end
    [~, k] = sortrows([R(cand,5), R(cand,15)]);
    idx = cand(k(1:min(N, numel(k))));
end

function idx = pick_min_mass(R, n_req, req_fat, Delta_lim_mm, th_lim_deg)
    idx = pick_top(R, n_req, req_fat, Delta_lim_mm, th_lim_deg, 1);
end

function s = sect_str(R, i)
    if isempty(i), s = 'none feasible'; else, s = sprintf('%g x %g x %g mm (%.2f kg)', R(i,2), R(i,3), R(i,4), R(i,5)); end
end

function s = dim_str(r)
    if r(1) == 3, s = sprintf('%g x %g solid', r(2), r(3));
    elseif r(1) == 2, s = sprintf('D%g x %g', r(2), r(4));
    elseif r(1) == 4, s = sprintf('%gx%g tf%g tw%g', r(2), r(3), r(4), r(16));
    else, s = sprintf('%g x %g x %g', r(2), r(3), r(4)); end
end

function print_row(label, r, fam_names)
    if r(1) == 3, tstr = 'solid';
    elseif r(1) == 4, tstr = sprintf('%.0f/%.0f', r(4), r(16));   % flange/web
    else, tstr = sprintf('%.1f', r(4)); end
    fprintf(' %-10s | %-10s | %4.0f x %-4.0f | %-9s | %8.2f  | %7.1f  | %7.2f  | %8.2f  | %9.2f\n', ...
        label, fam_names{r(1)}, r(2), r(3), tstr, r(5), r(6), r(7), r(8), r(15));
end

function print_detail(label, r, lim)
    fprintf(' %-10s | %18.3f | %16.4f | %14.3f | %17.3f | %12.3f | %9.2f\n', label, r(14), r(10), r(9), r(11), r(15), lim);
end