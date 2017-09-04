% plot results of calculations from crustalfluidmodel.R

%% Load Model OUtput

fid = fopen('logaeq2e.csv');            % even numbered
    hdr = textscan(fid,'%s',1,'HeaderLines',0)
    fclose(fid);
heads = strsplit(cell2mat(hdr{1}),'","')
heads = heads(2:end)
heads{end} = heads{end}(1:end-1) % get rid of trailing character

sp = heads;
nspecies = length(sp);

conds = csvread('conds2e.csv', 1,1);       % [T, P, logfO2]
    temp = csvread('conds2o.csv', 1,1);
%     conds = [conds; temp]
logaeq = csvread('logaeq2e.csv', 1,1);     % logact [graphite, CO, CO2, ... propane]
    temp = csvread('logaeq2e.csv', 1,1);     % logact [graphite, CO, CO2, ... propane]
%     logaeq = [logaeq; temp]

[A,index] = sortrows(conds,[1,3]);
B = logaeq(index,:);

uniqueT = unique(A(:,1))
% uniqueP = unique(A(:,2))
uniquefO2 = unique(A(:,3))

bb = reshape(B, [length(uniquefO2), length(uniqueT), nspecies])  % [logfO2, T, nspecies]
bb = permute(bb, [2 1 3])   % [T, logfO2, nspecies]

% grid on [T, fO2, logact]
[X, Y] = ndgrid(uniqueT, uniquefO2)

Ti = 300:100:700;
logfO2i = -42:1:-28;
[xi, yi] = ndgrid(Ti, logfO2i)

%% plot gridded results

figure(2); clf;
load('mycolormap.mat')   % contains mycmap

sptoplot = [2 3 4 5 6 8 9];     % omit Cgr and O2

zi = [];
for ii = 1:length(sptoplot)
    zi(:,:,ii) = interpn(X, Y, bb(:,:,sptoplot(ii)), xi, yi)
    mesh(xi,yi,zi(:,:,ii)); hold on;
end

for mm = 1:length(sptoplot)
    xloc = find(xi == max(xi(:,1)),1);
    yloc = find(yi == min(yi(1,:)),1);
    plot3(xi(xloc,yloc), yi(xloc,yloc), zi(xloc, yloc, mm), 'o', 'Color', [0.500000 0.500000 0.500000], 'MarkerSize', 4);
    text(xi(xloc,yloc)+5, yi(xloc,yloc)-1, zi(xloc, yloc, mm), sp{sptoplot(mm)}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle',...
        'Color', [0.500000 0.500000 0.500000], 'FontSize', 8)
end

hold off;

xlabel(['Temperature, ' char(176) 'C'])
ylabel('log{\it f}_{O_2}')
zlabel('log{\it f}')

title('Composition of graphite-saturated C-O-H fluid at 500 bar')

% zlim([-10, 4])
colormap(mycmap)
colorbar

view(-33, 9)

print(gcf(), '-depsc2', '-loose', 'cfm_variablefO2_500bar.eps');


%% make and plot interpolated results at buffered fO2

plotbufers;  % load and plot the buffers model output from T = 300 to 700 C

thebuffer = 'FMQ';
deltabuff = +0;     % Delta logfO2 vs. buffer

buffidx = find(strcmp(bn, 'FMQ'),1);

Tbi = [300:5:700]';
logfO2bi = interp1(Tint', bs(:,buffidx), Tbi, 'spline', NaN) + deltabuff;

logactbi = [];    % [nTbi x nsptoplot]
for ii = 1:length(sptoplot)
    logactbi(:,ii) = interpn(X, Y, bb(:,:,sptoplot(ii)), Tbi, logfO2bi, 'spline', NaN)
%     F = griddedInterpolant(X, Y, bb(:,:,sptoplot(ii)));
%     logactbi2(:,ii) = F(Tbi, logfO2bi);
%     isequaln(logactbi, logactbi2) % TRUE
end


figure(3); clf;

subplot(2,1,1);
plot(Tbi, logactbi)

ylabel('log{\it f}')
xlabel(['Temperature, ' char(176) 'C'])

subplot(2,1,2);
plot(Tbi, logfO2bi)

ylabel('log{\it f}_{O_2}')
xlabel(['Temperature, ' char(176) 'C'])

% zi = [];
% for ii = 1:length(sptoplot)
%     zi(:,:,ii) = interpn(X, Y, bb(:,:,sptoplot(ii)), xi, yi)
%     mesh(xi,yi,zi(:,:,ii)); hold on;
% end

