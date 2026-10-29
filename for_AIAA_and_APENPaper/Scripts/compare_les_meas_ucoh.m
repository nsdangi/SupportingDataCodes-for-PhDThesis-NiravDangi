clear; clc; close all;

% Define separation values and frequency range
separation_values = 20:20:200; % Separation values in meters
addpath("..\Data\")
% Base directory for data files
base_dir = "..\Data\ComparingCoh_shifting_FINAL\";
KaimalFuncYZ = @(r, f, V_hub, L_c) exp(-12 * sqrt((f .* r / V_hub) .^ 2 + (0.12 * r / L_c) .^ 2)); % Kaimal YZ
les_freq = logspace(log10(0.001), log10(0.1), 30)';
les_freq2 = logspace(log10(0.0005), log10(1), 50)';
coherenceFunc_LAT = @(params, r, f, Vhub) ...
    exp(-(params(1) .* ((r .* f ./ Vhub).^params(2)+(0.05.*r./340.2).^params(2))));
paramsN=[35	1.16];
paramsS=[8	0.81];
paramsUS=[39	1.36];
colors = gray(4);
markers={'<','*','o','pentagram'};
L_c = 340.2; 
%%
% Function to plot coherence data
function plot_coherence(nabl_file, les_u_file, les_v_file, row_idx, title_str,subtilte_str)
    % Read NABL dataset
    colors = gray(4);
    markers={'<','*','o','pentagra\,m'};
    data_nabl = readtable(nabl_file, FileType="spreadsheet", VariableNamingRule="preserve", ReadVariableNames=true);
    les_freq = logspace(log10(0.001), log10(0.1), 30)';

    % Read LES data
    data_les_U = readtable(les_u_file);
    data_les_V = readtable(les_v_file);

    % Create figure
    figure();
    semilogx(data_nabl.("Frequency (f_syed)"), data_nabl.Coh_m_0, 'Color',colors(1,:),'Marker',markers{1}, 'LineWidth', 2); hold on;
    
    
    G_fit_U = find_G_fit(data_nabl.Coh_m_0', (table2array(data_les_U(row_idx, :)).^2), 6);
    G_fit_V = find_G_fit(data_nabl.Coh_m_1', (table2array(data_les_V(row_idx, :)).^2), 6);
    
    % Formatting
    xlim([0.0005 1]);
    xticks([1e-3 1e-2 1e-1 1e0])
    ylim([0 1])
    yticks([0:0.2:1])
    grid on;
    xlabel('$f (Hz)$','Interpreter','latex');
    ylabel('$\gamma^{2}_{uu,lat} (-)$','Interpreter','latex');
    title(['$',subtilte_str,'; ',title_str,'$'],'Interpreter','latex');
    ax=gca;
    ax.FontSize=14;
    ax.TickLabelInterpreter='latex';
    ax.LineWidth=2;
    ax.MinorGridLineWidth=1;
    grid on
    box on
   
end
mkdir ../Images/
%% USABL
plot_coherence("Sr_no_17.xlsx", base_dir + "MSCoh_U_Mann_Alt285m_BR_USABL.txt", base_dir + "MSCoh_V_Mann_Alt285m_BR_USABL.txt", 5, '\Delta_{Y}=114\,m','U :: z= 287\,m');

kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(114, les_freq2, 8.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsUS, 114, les_freq2, 8.56);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');

figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)

exportgraphics(gcf, '..\Images\compsr17.png')
close

plot_coherence("Sr_no_20.xlsx", base_dir + "MSCoh_U_Mann_Alt265m_BR_USABL.txt", base_dir + "MSCoh_V_Mann_Alt265m_BR_USABL.txt", 7, '\Delta_{Y}=150\,m','U :: z= 267\,m');
kaimalYZ = (KaimalFuncYZ(150, les_freq2, 8.5, L_c)).^2;
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},'LineWidth',3);
heloocoh=coherenceFunc_LAT(paramsUS, 150, les_freq2, 8.58);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr20.png')
close

plot_coherence("Sr_no_26.xlsx", base_dir + "MSCoh_U_Mann_Alt265m_AR_USABL.txt", base_dir + "MSCoh_V_Mann_Alt265m_AR_USABL.txt", 10, '\Delta_{Y}=200\,m','U :: z= 258\,m');
kaimalYZ = (KaimalFuncYZ(200, les_freq2, 14.5, L_c)).^2;
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},'LineWidth',3);
heloocoh=coherenceFunc_LAT(paramsUS, 200, les_freq2, 13.82);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr26.png')
close

plot_coherence("Sr_no_40.xlsx", base_dir + "MSCoh_U_Mann_Alt265m_AR_USABL.txt", base_dir + "MSCoh_V_Mann_Alt265m_AR_USABL.txt", 9, '\Delta_{Y}=180\,m','U :: z= 270\,m');
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 14.5, L_c)).^2;
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},'LineWidth',3);
heloocoh=coherenceFunc_LAT(paramsUS, 180, les_freq2, 13.82);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr40.png')
close

%% NABL
close all

plot_coherence("Sr_no_11.xlsx", base_dir + "MSCoh_U_Mann_Alt165m_R_NABL.txt", base_dir + "MSCoh_V_Mann_Alt165m_R_NABL.txt", 4, '\Delta_{Y}=114\,m','N :: z= 169\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(114, les_freq2, 11.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 114, les_freq2, 11.43);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)

exportgraphics(gcf, '..\Images\compsr11.png')
close

plot_coherence("Sr_no_16.xlsx", base_dir + "MSCoh_U_Mann_Alt285m_R_NABL.txt", base_dir + "MSCoh_V_Mann_Alt285m_R_NABL.txt", 4, '\Delta_{Y}=114\,m','N :: z= 287\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(114, les_freq2, 11.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 114, les_freq2, 11.99);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr16.png')
close

plot_coherence("Sr_no_31.xlsx", base_dir + "MSCoh_U_Mann_Alt165m_AR_NABL.txt", base_dir + "MSCoh_V_Mann_Alt165m_AR_NABL.txt", 8, '\Delta_{Y}=180\,m','N :: z= 162\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 180, les_freq2, 14.53);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr31.png')
close

plot_coherence("Sr_no_33.xlsx", base_dir + "MSCoh_U_Mann_Alt165m_AR_NABL.txt", base_dir + "MSCoh_V_Mann_Alt165m_AR_NABL.txt", 8, '\Delta_{Y}=180\,m','N :: z= 162\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 180, les_freq2, 14.53);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr33.png')
close

plot_coherence("Sr_no_34.xlsx", base_dir + "MSCoh_U_Mann_Alt165m_R_NABL.txt", base_dir + "MSCoh_V_Mann_Alt165m_R_NABL.txt", 8, '\Delta_{Y}=180\,m','N :: z= 162\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 11.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 180, les_freq2, 11.43);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr34.png')
close

plot_coherence("Sr_no_35.xlsx", base_dir + "MSCoh_U_Mann_Alt165m_AR_NABL.txt", base_dir + "MSCoh_V_Mann_Alt165m_AR_NABL.txt", 8, '\Delta_{Y}=180\,m','N :: z= 162\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 180, les_freq2, 14.53);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr35.png')
close

plot_coherence("Sr_no_39.xlsx", base_dir + "MSCoh_U_Mann_Alt265m_BR_NABL.txt", base_dir + "MSCoh_V_Mann_Alt265m_BR_NABL.txt", 8, '\Delta_{Y}=180\,m','N :: z= 272\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 8.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsN, 180, les_freq2, 8.83);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr39.png')
close

%% SABL
close all
plot_coherence("Sr_no_9.xlsx", base_dir + "COCOHCoh_U_LES_Alt165_AR_SABL.txt", base_dir + "COCOHCoh_V_LES_Alt165_AR_SABL.txt", 4, '\Delta_{Y}=114\,m','S :: z= 169\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(114, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence
semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsS, 114, les_freq2, 14.46);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)

exportgraphics(gcf, '..\Images\compsr9.png')
close

plot_coherence("Sr_no_10.xlsx", base_dir + "COCOHCoh_U_LES_Alt165_AR_SABL.txt", base_dir + "COCOHCoh_V_LES_Alt165_AR_SABL.txt", 4, '\Delta_{Y}=114\,m','S :: z= 169\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(114, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence

semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsS, 114, les_freq2, 14.46);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)
exportgraphics(gcf, '..\Images\compsr10.png')
close

plot_coherence("Sr_no_36.xlsx", base_dir + "COCOHCoh_U_LES_Alt165_AR_SABL.txt", base_dir + "COCOHCoh_V_LES_Alt165_AR_SABL.txt", 8, '\Delta_{Y}=180\,m','S :: z= 162\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(180, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence

semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsS, 180, les_freq2, 14.46);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)

exportgraphics(gcf, '..\Images\compsr36.png')
close

plot_coherence("Sr_no_47.xlsx", base_dir + "COCOHCoh_U_LES_Alt165_AR_SABL.txt", base_dir + "COCOHCoh_V_LES_Alt165_AR_SABL.txt", 10, '\Delta_{Y}=224\,m','S :: z= 159\,m');
kaimalYZ=[];
kaimalYZ = (KaimalFuncYZ(224, les_freq2, 14.5, L_c)).^2; % Compute Kaimal coherence

semilogx(les_freq2,kaimalYZ(:,1),'Color',colors(2,:),'Marker',markers{2},LineWidth=2);
heloocoh=coherenceFunc_LAT(paramsS, 224, les_freq2, 14.46);
semilogx(les_freq2,heloocoh,'Color',colors(3,:),'Marker',markers{3},LineWidth=2)
hold off;
legend({'$Measured$','$IEC$','$Model$'}, 'Location', 'northeast','Interpreter','latex');
figureHandle = gcf;
figureHandle.Units = 'inches'; 

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 21)

exportgraphics(gcf, '..\Images\compsr47.png')
close
%%
function G_fit = find_G_fit(cohs_m, cohs_t, n)
    
    alpha_1s = zeros(1, length(cohs_m)); % Preallocate alpha_1s
    alpha_2s = zeros(1, length(cohs_m)); % Preallocate alpha_2s
    
    for j = 1:length(cohs_t)
        coh = cohs_t(j);
        alpha_1s(j) = 1 - (n-1)/n * (1-coh)^n * hypergeometric2f1(n, n, n+1, coh,100);
        alpha_2s(j) = 1 - (1-coh)^n * (n-1) * (n/(n+1)*hypergeometric2f1(n+1, n, n+2, coh,100) - (n-2)/n*hypergeometric2f1(n, n, n+1, coh,100));
    end
    
    var_theoretical = alpha_2s - alpha_1s.^2;
    var_measured = (alpha_1s - cohs_m).^2;
    
    G_fit = mean(var_measured) / mean(var_theoretical); % Calculate G_fit for i=1
end



function z=hypergeometric2f1(a,b,c,x,n)
% HYPERGEOMETRIC2F1 Computes the hypergeometric function 
% using a series expansion:
%
%    f(a,b;c;x)=
%
%    1 + [ab/1!c]x + [a(a+1)b(b+1)/2!c(c+1)]x^2 +
%    [a(a+1)(a+2)b(b+1)(b+2)/3!c(c+1)(c+2)]x^3 + ...
%
% The series is expanded to n terms
%
% This function solves the Gaussian Hypergeometric Differential Equation:
%
%     x(1-x)y'' + {c-(a+b+1)x}y' - aby = 0
%
% The Hypergeometric function converges only for:
% |x| < 1
% c != 0, -1, -2, -3, ...
%
%
% Comments to:
% Diego Garcia   - d.garcia@ieee.org
% Chuck Mongiovi - mongiovi@fast.net
% June 14, 2002
if nargin ~= 5
    error('Usage: hypergeometric2f1(a,b,c,x,n) --> Wrong number of arguments')
end
if (n <= 0 | n ~= floor(n))
    error('Usage: hypergeometric2f1(a,b,c,x,n) --> n has to be a positive integer')
end
if (abs(x) > 1)
    error('Usage: hypergeometric2f1(a,b,c,x,n) --> |x| has to be less than 1')
end
if (c <= 0 & c == floor(c))
    error('Usage: hypergeometric2f1(a,b,c,x,n) --> c != 0, -1, -2, -3, ...')
end
z = 0;
m = 0;
while (m<n)
    if (m == 0)
        delta = 1;
    else
        delta = delta .* x .* (a + (m - 1)) .* (b + (m-1)) ./ m ./ (c + (m-1));
   end
   z = z + delta;
   m = m + 1;
end
end

