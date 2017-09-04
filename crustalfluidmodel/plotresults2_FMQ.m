% plot results of calculations from crustalfluidmodel.R
% 1000bar FMQ, graphite satureated

strsplit = @(str,delim) regexp(str,regexptranslate('escape',delim),'split')
% one-line replacement function for strsplit: https://stackoverflow.com/a/35325913

%% Load Model OUtput

fid = fopen('logaeq2e.csv');            % even numbered
    hdr = textscan(fid,'%s',1,'HeaderLines',0)
    fclose(fid);
heads = strsplit(cell2mat(hdr{1}),'","')
heads = heads(2:end)
heads{end} = heads{end}(1:end-1) % get rid of trailing character

sp = heads;
nspecies = length(sp);

conds = csvread('conds2_FMQ.csv', 1,1);       % [T, P, logfO2]
%     temp = csvread('conds2o.csv', 1,1);
%     conds = [conds; temp]
logaeq = csvread('logaeq2_FMQ.csv', 1,1);     % logact [graphite, CO, CO2, ... propane]
%     temp = csvread('logaeq2e.csv', 1,1);     % logact [graphite, CO, CO2, ... propane]
%     logaeq = [logaeq; temp]

steamline = csvread('steamline_1000bar.csv',1,1); %vapor pressure of liq water at 1000bar

TPFMQ =  csvread('conds2_temp.csv', 1,1);

sptoplot = [2 3 4 5 6 8 9];     % omit Cgr and O2

Tint = [120:10:600]';
logaeqint = interp1(conds(:,1), logaeq(:,:), Tint, 'linear');

xl = [00 600];

CH4_CO2 = logaeqint(:,4)-logaeqint(:,3)
CH4_H2O = logaeqint(:,4)-logaeqint(:,6)
CH4_C2 = logaeqint(:,4)-logaeqint(:,8)
C2_C3 = logaeqint(:,8)-logaeqint(:,9)

%% PLot

figure(1); clf
hs = tight_subplot(3,2,[0.03 0.11],[0.11 0.05], [0.12 0.05])


axes(hs(1))
    plot(Tint(Tint<500), logaeqint(Tint<500,sptoplot), 'k-', 'LineWidth', 0.5); hold on;
    plot(Tint(Tint>=500), logaeqint(Tint>=500,sptoplot), 'k--', 'LineWidth', 0.5); hold on;
    
    plot(steamline(:,1), steamline(:,4), '-', 'Color', [0.280000 0.570000 0.810000]); hold on;
    plot([485; 485], [-10 4], '--', 'Color', [0.500000 0.500000 0.500000], 'LineWidth', 0.25); hold on;
    plot([xl; xl]', [3 3; 3 3]', 'k:', 'LineWidth', 0.25); hold off

    for mm = 1:length(sptoplot)
        text(Tint(1), logaeqint(1,sptoplot(mm)), sp{sptoplot(mm)}, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle',...
            'Color', [0.000000 0.000000 0.000000], 'FontSize', 6)
    end
    
    ht = text(550,-3, 'graphite unstable', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom',...
        'Color', [0.500000 0.500000 0.500000], 'FontSize', 8)
    set(ht, 'rotation', 90)
    
    ylabel('log \it{f}')
    xlim(xl)
    
%     legend({sp{sptoplot}}, 'Location', 'best', 'FontSize', 6)
    set(gca(),'XTickLabel',[])
    set(gca(),'TickLength',3*get(gca(),'TickLength'))
    set(gca(),'XMinorTick','on','YMinorTick','on')
    
    axpos = get(gca(),'Position')
    set(gca(),'Position', [axpos(1), axpos(2)-0.29, axpos(3), axpos(4)+0.29])
    
    yl=ylim  %panel label
    text(xl(1)+0.05*diff(xl),diff(yl)*0.97+yl(1),'A', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    
    
axes(hs(3))
    set(gca(), 'Visible', 'off')
    
axes(hs(5))
    plot(TPFMQ(:,1), TPFMQ(:,3), '-', 'Color', [0.000000 0.000000 0.000000])

    ylabel('log {\it{f}}_{O_2}')
    xlabel(['Temperature, ' char(176) 'C'])
    xlim(xl)
    
    set(gca(),'TickLength',3*get(gca(),'TickLength'))
    set(gca(),'XMinorTick','on','YMinorTick','on')
    
    yl=ylim  %panel label
    text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'B', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    
    
axes(hs(2))
    plot(Tint, CH4_CO2, 'k-'); hold on;
    plot(Tint, CH4_H2O, '-', 'Color', [0.500000 0.500000 0.500000]); hold off
    
    text(Tint(1)+10, CH4_CO2(1), 'log {\it{f}}_{CH_4}/{\it{f}}_{CO_2}', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
        'Color', [0.000000 0.000000 0.000000], 'FontSize', 6)
    text(Tint(1)+10, CH4_H2O(1), 'log {\it{f}}_{CH_4}/{\it{f}}_{H_2O}', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top',...
        'Color', [0.500000 0.500000 0.500000], 'FontSize', 6)
    
    ylabel('log {\it{f}}/{\it{f}}')
    xlim(xl)
    
    set(gca(),'XTickLabel',[])
    
    set(gca(),'TickLength',3*get(gca(),'TickLength'))
    set(gca(),'XMinorTick','on','YMinorTick','on')
    
    yl=ylim  %panel label
    text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'C', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    
    
axes(hs(4))
    plot(Tint, CH4_C2, 'k-'); hold on
    plot([xl]', log10(950)*ones(1,2)', 'k:', 'LineWidth', 0.25); hold on % plot C1/C2 range from Charlou et al (2000, 2002), Mcdermott et al., PNAS, and Proskurowosk et al Science
    plot([xl]', log10(4000)*ones(1,2)', 'k:', 'LineWidth', 0.25); hold off
    
    ylabel('log (CH_4/C_2H_6)')
    xlim(xl)
    ylim([2 5])
    
    set(gca(),'XTickLabel',[])
    
    set(gca(),'TickLength',3*get(gca(),'TickLength'))
    set(gca(),'XMinorTick','on','YMinorTick','on')
    
    yl=ylim  %panel label
    text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'D', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    
axes(hs(6))
    plot(Tint, C2_C3, 'k-'); hold on
    plot([xl]', log10(10)*ones(1,2)', 'k:', 'LineWidth', 0.25); hold on % plot C2/C3 range from Charlou et al (2000, 2002), Mcdermott et al., PNAS, and Proskurowosk et al Science
    plot([xl]', log10(1097/48)*ones(1,2)', 'k:', 'LineWidth', 0.25); hold off
    
    ylabel('log (C_2H_6/C_3H_8)')
    xlim(xl)
    ylim([0 4])
    
    set(gca(),'TickLength',3*get(gca(),'TickLength'))
    set(gca(),'XMinorTick','on','YMinorTick','on')
    
    yl=ylim  %panel label
    text(xl(1)+0.05*diff(xl),diff(yl)*0.95+yl(1),'E', 'FontSize', 11, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
    
    
    xlabel(['Temperature, ' char(176) 'C'])
    

set(gcf, 'Position', [621   152   653   511])
print(gcf(), '-depsc2', '-loose', 'crustalfluidmodel_FMQ_1000bar(plotresults2_FMQ).eps');
