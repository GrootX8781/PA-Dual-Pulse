

set(groot, 'defaultAxesFontName',   'Times New Roman');
set(groot, 'defaultTextFontName',   'Times New Roman');
set(groot, 'defaultAxesFontSize',   22);
set(groot, 'defaultTextFontSize',   22);
set(groot, 'defaultAxesLineWidth',  2);
set(groot, 'defaultLineLineWidth',  2);
set(groot, 'defaultLegendFontSize', 18);
set(groot, 'defaultScatterMarkerFaceAlpha', 0.25);
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultColorbarTickLabelInterpreter', 'latex');


p       = [300 200 840 630];

lw      = 2;                
msz     = 6;

cols_ln = {[0.00 0.45 0.74], [0.85 0.33 0.10], [0.49 0.18 0.56], ...
           [0.00 0.50 0.50], [0.70 0.16 0.16], [0.80 0.60 0.00], ...
           [0.30 0.45 0.55], [0.40 0.60 0.20], [0.55 0.27 0.07], ...
           [0.00 0.45 0.63], [0.00 0.45 0.74], [0.85 0.33 0.10]};

cTh     =  [0.35 0.35 0.35];

cols_pt = {[0.80 0.88 0.97], [0.99 0.92 0.82], [0.93 0.86 0.96], ...
           [0.83 0.94 0.93], [0.98 0.86 0.86], [0.99 0.95 0.85], ...
           [0.92 0.94 0.96], [0.92 0.97 0.86], [0.95 0.89 0.85], ...
           [0.85 0.95 0.98], [0.80 0.88 0.97], [0.99 0.92 0.82]};
       
ls      = {'--', '-.', '-', ':', ...
           '--', '-.', '-', ':', ...
           '--', '-.', '-', ':'};
