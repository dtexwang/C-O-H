function cm = makelinecmap(ncolors)
% generate good, perceptually-spaced cmap for lines
% idea from http://stackoverflow.com/a/13039532/1917448
% colors from http://vis.stanford.edu/color-names/analyzer/
% converted to rgb using http://www.zonums.com/online/color_converter/

cmaps = {};

cmaps{1} = [98	30	21
            229	144	118
            18	141	205
            8	60	82
            100	197	242
            97	175	175
            15	115	105
            156	157	161];  % TheEconomist

cmaps{2} = [156	158	222
            115	117	181
            74	85	132
            206	219	156
            181	207	107
            140	162	82
            99	121	57
            231	203	148
            231	186	82
            189	158	57
            140	109	49
            231	150	156
            214	97	107
            173	73	74
            132	60	57
            222	158	214
            206	109	189
            165	81	148
            123	65	115];  % ManyEyes

cmaps{3}  =[31	119	180
            174	199	232
            255	127	14
            255	187	120
            44	160	44
            152	223	138
            214	39	40
            255	152	150
            148	103	189
            197	176	213
            140	86	75
            196	156	148
            227	119	194
            247	182	210
            127	127	127
            199	199	199
            188	189	34
            219	219	141
            23	190	207
            158	218	229];  % Tableau20

sizes = [];
for ii = 1:length(cmaps)
    sizes(ii) = size(cmaps{ii},1);
end

icmap = find(sizes>ncolors, 1);

ic = round([1:ncolors]/length([1:ncolors])*size(cmaps{icmap},1));
cm = cmaps{icmap}(ic,:);
cm = cm./255;
