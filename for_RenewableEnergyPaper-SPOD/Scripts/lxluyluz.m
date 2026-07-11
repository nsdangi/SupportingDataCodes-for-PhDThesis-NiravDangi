clear
clc
load('..\Data\3Dvelarray.mat_3Dup.mat')

%%

for jj=1:length(INP_common.subfolders)
    subfolder=INP_common.subfolders{jj};
    disp(subfolder)
    if jj == 4 || jj == 5 || jj == 6
        dt=0.33;
    else
        dt=0.25;
    end
    X1_u=OPT_3dVel.(subfolder).vel_arr_U;

    for ss=1:size(X1_u,2)

        for aa=1:size(X1_u,2)
            lxhere(ss,aa)=Lx(detrend(squeeze(X1_u(:,ss,aa))),mean(squeeze(X1_u(:,ss,aa)),1),dt);


        end
    end

    LXstabs(:,jj)=mean(lxhere,1);

end
%
for jj=1:length(INP_common.subfolders)
    subfolder=INP_common.subfolders{jj};
    disp(subfolder)
    if jj == 4 || jj == 5 || jj == 6
        dt=0.33;
    else
        dt=0.25;
    end
    X1_u=OPT_3dVel.(subfolder).vel_arr_U;



    for aa=1:size(X1_u,2)
        luyhere(aa)=Lyz(detrend(squeeze(X1_u(:,:,aa))),unique(points.Y));

    end

    LUYstabs(:,jj)=luyhere;

end


%
for jj=1:length(INP_common.subfolders)
    subfolder=INP_common.subfolders{jj};
    disp(subfolder)
    if jj == 4 || jj == 5 || jj == 6
        dt=0.33;
    else
        dt=0.25;
    end
    X1_u=OPT_3dVel.(subfolder).vel_arr_U;



    for aa=1:size(X1_u,2)
        luzhere(aa)=Lyz(detrend(squeeze(X1_u(:,aa,:))),unique(points.Z));

    end

    LUZstabs(:,jj)=luzhere;

end
%%
ieclx=readmatrix('SEEDMEANiecratedlx.txt');
iecluy=readmatrix('SEEDMEANiecratedluy.txt');
iecluz=readmatrix('SEEDMEANiecratedluz.txt');

%% ONLY RATED
close all

abris = [1,2,3];
bris  = [4,5,6];
atris = [7,8,9];

colors = gray(5);
markershere = {'<','*','o','pentagram'};

nMarkers = 5;

%% ===================== FIGURE 1 =====================
figure(1); hold on

xvals = unique(points.Z);
N = length(xvals);
nLines = length(bris) + 1; % + IEC

for jj = 1:length(bris)
    
    y = abs(smooth(filloutliers(LXstabs(:, atris(jj)),'nearest')));
    
    base = linspace(1, N, nMarkers);
    offset = (jj-1) * N/(nMarkers*nLines);
    idxmarker = round(base + offset);
    idxmarker(idxmarker > N) = N;
    
    plot(y, xvals, ...
        'Color', colors(jj,:), ...
        'LineWidth', 2, ...
        'Marker', markershere{jj}, ...
        'MarkerIndices', idxmarker, ...
        'MarkerFaceColor','none', ...
        'MarkerSize',6);
end

% IEC line
y = abs(smooth(filloutliers(ieclx,'nearest')));
base = linspace(1, N, nMarkers);
offset = (nLines-1) * N/(nMarkers*nLines);
idxmarker = round(base + offset);
idxmarker(idxmarker > N) = N;

plot(y, xvals, ...
    'Color', colors(4,:), ...
    'LineWidth', 2, ...
    'Marker', markershere{4}, ...
    'MarkerIndices', idxmarker, ...
    'MarkerFaceColor','none', ...
    'MarkerSize',6);

xlim([0 300])
xlabel('$\mathcal{L}_{u,x} [m]$','Interpreter','latex')
ylabel('$z/D [-]$','Interpreter','latex')

yticks((0:0.25:1.5)*280)
yticklabels({'0','0.25','0.5','0.75','1','1.25'})
yline(170,'k--','LineWidth',2)


%% ===================== FIGURE 2 =====================
figure(2); hold on

for jj = 1:length(bris)
    
    y = abs(smooth(filloutliers(LUYstabs(:, atris(jj)),'nearest')));
    
    base = linspace(1, N, nMarkers);
    offset = (jj-1) * N/(nMarkers*nLines);
    idxmarker = round(base + offset);
    idxmarker(idxmarker > N) = N;
    
    plot(y, xvals, ...
        'Color', colors(jj,:), ...
        'LineWidth', 2, ...
        'Marker', markershere{jj}, ...
        'MarkerIndices', idxmarker, ...
        'MarkerFaceColor','none', ...
        'MarkerSize',6);
end

% IEC line
y = abs(smooth(filloutliers(iecluy,'nearest')));
base = linspace(1, N, nMarkers);
offset = (nLines-1) * N/(nMarkers*nLines);
idxmarker = round(base + offset);
idxmarker(idxmarker > N) = N;

plot(y, xvals, ...
    'Color', colors(4,:), ...
    'LineWidth', 2, ...
    'Marker', markershere{4}, ...
    'MarkerIndices', idxmarker, ...
    'MarkerFaceColor','none', ...
    'MarkerSize',6);

xlim([0 100])
xlabel('$\mathcal{L}_{u,y} [m]$','Interpreter','latex')
ylabel('$z/D [-]$','Interpreter','latex')

yticks((0:0.25:1.5)*280)
yticklabels({'0','0.25','0.5','0.75','1','1.25'})
yline(170,'k--','LineWidth',2)


%% ===================== FIGURE 3 =====================
figure(3); hold on

yvals = unique(points.Y);
N = length(yvals);

for jj = 1:length(bris)
    
    y = abs(smooth(filloutliers(LUZstabs(:, atris(jj)),'nearest')));
    
    base = linspace(1, N, nMarkers);
    offset = (jj-1) * N/(nMarkers*nLines);
    idxmarker = round(base + offset);
    idxmarker(idxmarker > N) = N;
    
    plot(y, yvals, ...
        'Color', colors(jj,:), ...
        'LineWidth', 2, ...
        'Marker', markershere{jj}, ...
        'MarkerIndices', idxmarker, ...
        'MarkerFaceColor','none', ...
        'MarkerSize',6);
end

% IEC line
y = abs(smooth(filloutliers(iecluz,'nearest')));
base = linspace(1, N, nMarkers);
offset = (nLines-1) * N/(nMarkers*nLines);
idxmarker = round(base + offset);
idxmarker(idxmarker > N) = N;

plot(y, yvals, ...
    'Color', colors(4,:), ...
    'LineWidth', 2, ...
    'Marker', markershere{4}, ...
    'MarkerIndices', idxmarker, ...
    'MarkerFaceColor','none', ...
    'MarkerSize',6);

xlim([0 100])
xlabel('$\mathcal{L}_{u,z} [m]$','Interpreter','latex')
ylabel('$y/D [-]$','Interpreter','latex')

yticks(1505-280:0.25*280:1505+280)
yticklabels({'-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1'})
yline(1505,'k--','LineWidth',2)


%% ===================== STYLING =====================
for uu = 1:3
    figure(uu)
    legend('$NABL$','$SABL$','$USABL$','$IEC_{s}$', ...
        'Interpreter','latex','Location','southeast')
    
    ax = gca;
    ax.LineWidth = 3;
    ax.FontSize = 16;
    ax.TickLabelInterpreter = 'latex';
    
    grid off
    box on
end


% ===================== EXPORT =====================
mkdir Images
exportgraphics(figure(1),'Images\luxrene.png')
exportgraphics(figure(2),'Images\luyrene.png')
exportgraphics(figure(3),'Images\luzrene.png')