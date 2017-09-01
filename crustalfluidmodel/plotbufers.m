% Plot buffered fO2

fid = fopen('buffers.csv')
    hdr = textscan(fid, '%s', 1)
    fclose(fid)
    
heads = strsplit(cell2mat(hdr{1}), '","')
heads = heads(2:end)  % drop the first cell, which is ","
heads = heads(3:end)  % drop the next two coluns, which are T P
heads{end} = heads{end}(1:end-1) % get rid of trailing character
bn = heads;

buffs = csvread('buffers.csv',1,1)
Ts = unique(buffs(:,1));
Ps = unique(buffs(:,2));

%% interpolate at 10 C steps in T
Tint = min(Ts):10:max(Ts);
buffsint = interp1(buffs(:,1), buffs(:,3:end), Tint, 'spline');

[bs, ib] = sortrows(buffsint',1);   % sort by fO2
bs = fliplr(bs');
bn = fliplr(bn(ib));

%% plot logfO2 vs T
figure(1); clf

cm = makelinecmap(length(1:size(bs,2)));
ls = {'-', '--', ':'}; 
ls = repmat(ls,[1, ceil(size(bs,2)/length(ls))]);
ls(1:size(bs,2));

for i = 1:length(1:size(bs,2))
    plot(Tint, bs(:,i),'color',cm(i,:), 'LineStyle', ls{i}); hold on
end
hold off;

legend(bn{1:end}, 'location', 'best')  % first two elements of bn = T P

ylabel('log{\it f}_{O_2}')
xlabel(['Temperature, ' char(176) 'C'])

title('Oxygen fugacity maintained by mineral redox buffers')

xl = xlim;
yl = ylim;
text(xl(1)+0.02*diff(xl), yl(1)+0.98*diff(yl), ['P = ' num2str(Ps) ' bar'],... 
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');%, 'BackgroundColor', 'white', 'EdgeColor', 'k', 'LineWidth', 0.02);


axis square
grid on
% grid minor

colormap(gray)

