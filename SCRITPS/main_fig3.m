clc; clear; close all

%%
config;

set(groot, 'defaultAxesFontName',   'Times New Roman');
set(groot, 'defaultTextFontName',   'Times New Roman');
set(groot, 'defaultAxesFontSize',   10);
set(groot, 'defaultTextFontSize',   18);
set(groot, 'defaultAxesLineWidth',  2);
set(groot, 'defaultLineLineWidth',  2);
set(groot, 'defaultLegendFontSize', 18);
set(groot, 'defaultScatterMarkerFaceAlpha', 0.25);
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');

s = S(2,:);
[s1, s2, s3, s4] = deal(s(1), s(2), s(3), s(4));
[~, ~, bound] = sv_info(s1, s2, s3, s4);
S = [1,                 s1 * exp(1j * s3);
     s1 * exp(1j * s3), s2 * exp(1j * s4)];

pj = [cos(pi / 4); sin(pi / 4) * exp(1j * pi / 2)];

%%
K = 20;
tp = 2;
[pk_2, ~] = gen_pk(K, tp, S);
tp = 5;
[pk_5, ~] = gen_pk(K, tp, S);
tp = 8;
[pk_8, ~] = gen_pk(K, tp, S);

polar_map(pk_2, 2, S);


%%
function polar_map(pk, track, S)
    
    pal = okabe_ito();
    c_tx = pal.blue;
    c_rx = pal.orange;
    halo = [1 1 1];

    if nargin < 3 || isempty(S)
        
        P = jones2xyz(pk);

        f = figure('Color', 'w', ...
            'Name', sprintf('Track %s', track_name(track)));
        set(f, 'NumberTitle', 'off');
        ax = gca; cla(ax); hold(ax,'on');

        draw_poincare(ax, 'glass');
        
        halo = [1 1 1];  sz_halo = 56;  sz = 26;
        
        scatter3(ax, P(:, 1), P(:, 2), P(:, 3), sz_halo, halo, 'filled', ...
            'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor', 'none');
        scatter3(ax, P(:, 1), P(:, 2), P(:, 3), sz, pal.blue, 'filled', ...
            'MarkerFaceAlpha', 0.98, 'MarkerEdgeColor', [0 0 0] + 0.10, ...
            'LineWidth', 0.25);
        
        view(ax, [135 28]); rotate3d(ax, 'on');
        tight_pdf(ax, 600);

    else

        P  = jones2xyz(pk);
        t  = S * pk;  
        t = t ./ max(sqrt(sum(abs(t) .^ 2, 1)), 1e-12);
        T  = jones2xyz(t);

        f = figure('Color', 'w', ...
            'Name', sprintf('Track %d: %s', track, track_name(track)));
        set(f, 'NumberTitle', 'off', 'Position', [400 300 1000 500]);

        ax1 = subplot(1, 2, 1); cla(ax1); hold(ax1, 'on');
        draw_poincare(ax1, 'light');
        sz_halo = 56;  sz = 26;
        scatter3(ax1, P(:, 1), P(:, 2), P(:, 3), sz_halo, [1 1 1], 'filled', ...
                 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
        scatter3(ax1, P(:, 1), P(:, 2), P(:, 3), sz, c_tx, 'filled', ...
                 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor', [0 0 0] + 0.08, ...
                 'LineWidth', 0.25);
             
        view(ax1, [135 28]); rotate3d(ax1, 'on');

        ax2 = subplot(1, 2, 2); cla(ax2); hold(ax2, 'on');
        draw_poincare(ax2, 'light');
        scatter3(ax2, T(:, 1), T(:, 2), T(:, 3), sz_halo, [1 1 1], 'filled', ...
                 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');
        scatter3(ax2, T(:, 1), T(:, 2), T(:, 3), sz, c_rx, 'filled', ...
                 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor', [0 0 0] + 0.08, ...
                 'LineWidth', 0.25);
             
        view(ax2, [135 28]); rotate3d(ax2, 'on');
    
    end
end

function name = track_name(tp)
    switch tp
        case 1,    name = 'RND 随机';
        case 2,    name = 'HV 线极化';
        case 3,    name = 'RCL RC→LC';
        case 4,    name = 'GRD 均匀扫描';
        case 5,    name = 'EQR 赤道环';
        case 6,    name = 'TRG 双环';
        case 7,    name = 'ICO 正二十面体轨道';
        case 8,    name = 'FIB Fibonacci网格轨道 ';
        case 9,    name = 'TAM Tammes数值解';
        case 10,   name = 'COV Covering数值解';
        case 11,   name = 'OPT 最优轨道';
        otherwise, name = sprintf('轨道%d',tp);
    end
end

function tight_pdf(ax, dpi)

    if nargin < 3 || isempty(dpi), dpi = 600; end
    fig = ancestor(ax, 'figure');

    axis(ax, 'off'); grid(ax, 'off'); box(ax, 'off');
    set(ax, 'Units', 'normalized', 'Position', [0 0 1 1]);
    set(ax, 'LooseInset', [0 0 0 0]);
    
    daspect(ax, [1 1 1]); pbaspect(ax, [1 1 1]); axis(ax, 'vis3d');
    lim = 1.02; 
    xlim(ax, [-lim lim]); 
    ylim(ax, [-lim lim]); 
    zlim(ax, [-lim lim]);
    
    set(fig, 'Units', 'centimeters'); 
    pos = get(fig, 'Position');
    
    sz = min(pos(3:4)); 
    set(fig, 'Position', [pos(1:2) sz sz]);
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'auto');
    set(fig, 'PaperSize', [pos(3) pos(4)]);
    set(fig, 'InvertHardcopy', 'off');

end
