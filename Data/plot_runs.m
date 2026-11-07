clc; clear; close all;

runs        = [1, 3];        % log file indices
run_labels  = [1, 2];        % labels shown in the report (log 2 omitted)
valid_range = [2, 10];          % m/s, shaded on the time-series only
log_dir     = "logs/team9";
fig_dir     = "figures";

% --- Style --------------------------------------------------------------
font_name   = "Times New Roman";
font_size   = 8;                % pt
dpi         = 600;
w_double    = 7.16;             % in, IEEE two-column width
w_single    = 3.5;              % in, IEEE one-column width
c_sensor    = hex2rgb("#2a78d6");   % blue
c_gt        = hex2rgb("#eb6834");   % orange
c_band      = [0.88 0.88 0.88];     % valid-range shading
c_zero      = [0.25 0.25 0.25];

% --- Load ---------------------------------------------------------------
R = struct();
for k = 1:numel(runs)
    i = runs(k);
    s  = readtable(fullfile(log_dir, sprintf("Team_9_log%i_student.csv", i)));
    gt = readtable(fullfile(log_dir, sprintf("Team_9_log%i_vdas.csv",    i)));
    al = readtable(fullfile(log_dir, sprintf("Team_9_log%i_aligned.csv", i)));

    t0 = s.timestamp_epoch(1);                  % common time origin
    R(k).run    = run_labels(k);
    R(k).t_s    = s.timestamp_epoch  - t0;
    R(k).v_s    = s.value_raw_mps;
    R(k).t_gt   = gt.timestamp_epoch - t0;
    R(k).v_gt   = gt.value_raw_mps;
    R(k).t_al   = al.timestamp_epoch - t0;
    R(k).gt_al  = al.ground_truth_mps;
    R(k).err    = al.student_mps - al.ground_truth_mps;   % signed, m/s

    in_rng = R(k).gt_al > valid_range(1) & R(k).gt_al < valid_range(2);
    fprintf("Run %i  RMSE: %.3f m/s  (2-10 m/s: %.3f m/s)  bias: %+.3f m/s\n", ...
        R(k).run, rms(R(k).err), rms(R(k).err(in_rng)), mean(R(k).err));
end

% Shared limits so the two runs are directly comparable
v_max   = max(cellfun(@max, [{R.v_s}, {R.v_gt}]));
v_lim   = [0, ceil(v_max)];
e_max   = max(abs(vertcat(R.err)));
e_lim   = [-1 1] * ceil(e_max * 2) / 2;      % symmetric, 0.5 m/s steps

%% Figure 1: time-series + residuals ------------------------------------
fig = newFigure(w_double, 3, font_name, font_size);
tl  = tiledlayout(fig, 3, 2, TileSpacing="compact", Padding="compact");
tags = ["(a)", "(b)"];

for k = 1:numel(R)
    % Time-series (spans two rows)
    ax = nexttile(tl, k, [2 1]);
    hold(ax, "on");
    t_end = max(R(k).t_s(end), R(k).t_gt(end));
    hb = patch(ax, [0 t_end t_end 0], valid_range([1 1 2 2]), c_band, ...
        EdgeColor="none", DisplayName=sprintf("Valid range (%g\x2013%g m/s)", valid_range));
    hs = plot(ax, R(k).t_s,  R(k).v_s,  Color=c_sensor, LineWidth=0.5, ...
        DisplayName="Sensor (200 Hz)");
    hg = plot(ax, R(k).t_gt, R(k).v_gt, Color=c_gt, LineWidth=1.0, LineStyle="--", ...
        DisplayName="VDAS ground truth");
    hold(ax, "off");
    styleAxes(ax, font_name, font_size);
    xlim(ax, [0 t_end]); ylim(ax, v_lim);
    ax.XTickLabel = [];
    if k == 1, ylabel(ax, "Airspeed (m/s)"); end
    title(ax, sprintf("%s Run %i", tags(k), R(k).run), FontWeight="normal", ...
        FontSize=font_size, FontName=font_name);
    ax_top = ax;

    % Residual (one row, directly beneath)
    ax = nexttile(tl, 4 + k);
    hold(ax, "on");
    area(ax, R(k).t_al, R(k).err, 0, FaceColor=c_sensor, FaceAlpha=0.35, ...
        EdgeColor=c_sensor, LineWidth=0.5);
    yline(ax, 0, Color=c_zero, LineWidth=0.5);
    hold(ax, "off");
    styleAxes(ax, font_name, font_size);
    xlim(ax, [0 t_end]); ylim(ax, e_lim);
    xlabel(ax, "Time (s)");
    if k == 1, ylabel(ax, "Error (m/s)"); end
    linkaxes([ax_top, ax], "x");
end

lgd = legend([hs hg hb], Orientation="horizontal", Box="off", ...
    FontName=font_name, FontSize=font_size);
lgd.Layout.Tile = "north";

save_png(fig, fullfile(fig_dir, "run1_run2_timeseries.png"), dpi);

%% Figures 2-3: error vs airspeed scatter --------------------------------
gt_lim = [0, ceil(max(vertcat(R.gt_al)))];

for k = 1:numel(R)
    fig = newFigure(w_single, 2.4, font_name, font_size);
    ax  = axes(fig);
    hold(ax, "on");
    yline(ax, 0, Color=c_zero, LineWidth=0.5);
    scatter(ax, R(k).gt_al, R(k).err, 6, c_sensor, "filled", ...
        MarkerFaceAlpha=0.5, MarkerEdgeColor="none");
    hold(ax, "off");
    styleAxes(ax, font_name, font_size);
    xlim(ax, gt_lim); ylim(ax, e_lim);
    xlabel(ax, "VDAS airspeed (m/s)");
    ylabel(ax, "Error (m/s)");

    save_png(fig, fullfile(fig_dir, sprintf("run%i_error_scatter.png", R(k).run)), dpi);
end

fprintf("Saved figures to %s/\n", fig_dir);

%% Helpers ---------------------------------------------------------------
function fig = newFigure(w_in, h_in, font_name, font_size)
    fig = figure(Visible="off", Color="w", Units="inches", ...
        Position=[1 1 w_in h_in], PaperPositionMode="auto");
    set(fig, DefaultAxesFontName=font_name, DefaultAxesFontSize=font_size, ...
        DefaultTextFontName=font_name, DefaultTextFontSize=font_size);
end

function styleAxes(ax, font_name, font_size)
    set(ax, FontName=font_name, FontSize=font_size, LineWidth=0.5, ...
        Box="on", TickDir="in", XGrid="on", YGrid="on", ...
        GridColor=[0 0 0], GridAlpha=0.12, Layer="top");
    ax.LabelFontSizeMultiplier = 1;
    ax.TitleFontSizeMultiplier = 1;
end

function save_png(fig, path, dpi)
    exportgraphics(fig, path, Resolution=dpi, BackgroundColor="white");
    close(fig);
    fprintf("  wrote %s\n", path);
end

function rgb = hex2rgb(hex)
    hex = char(erase(hex, "#"));
    rgb = double([hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))]) / 255;
end
