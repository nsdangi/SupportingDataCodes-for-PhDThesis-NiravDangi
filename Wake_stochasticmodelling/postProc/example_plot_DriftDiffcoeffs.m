imgfolder = fullfile(pwd, 'images');
if ~exist(imgfolder, 'dir')
    mkdir(imgfolder);
end
%% DRIFT
atWhatbin = 8; % mid bin; because in FPE settings nBins = 15 was used 
figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        % if bin average
        % func_scatterToContour(coords{i,j}, nanmean(interpD1_all{i,j}{1}(7:10,:),1),50)
        func_scatterToContour(coords{i,j}, interpD1_all{i,j}{1}(atWhatbin,:),50)

        colormap(gray(10))
        clim([-0.1 0.1])
        cb = colorbar;
        cb.Ticks=[-0.25:0.05:0.25];
        cb.Label.String = '$\mathcal{D}_{1}^{u}~[m/s^2]$';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 14;
        cb.TickLabelInterpreter = 'latex';
       
        xticks([-1:0.5:1].*turbDia+turbHub(2))
        xlim([mask.y(1) mask.y(2)])
        yticks([0:0.25:1.25].*turbDia)
        ylim([0 max(mask.z(2))])
        yticklabels({'0','0.25','0.5','0.75','1','1.25'})
        xticklabels({'-1','-0.5','0','0.5','1'})
        ax=gca;
        ax.FontSize=14;
        ax.TickLabelInterpreter="latex";
        xlabel('$y/D[-]$',Interpreter='latex')
        ylabel('$z/D[-]$',Interpreter='latex')
        title(sprintf('$u_{\\infty,hub}=%.1f\\,\\mathrm{m/s}$', uinfhub(i)), ...
            'Interpreter', 'latex')
    end
end
exportgraphics(gcf,'images\D1_contour_midbin_yz5_diffOp.png')
%% DIFFUSION 
atWhatbin = 8;  % mid bin; because in FPE settings nBins = 15 was used
figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        % if bin average
        % func_scatterToContour(coords{i,j}, nanmean(interpD2_all{i,j}{1}(7:10,:),1),50)
        func_scatterToContour(coords{i,j}, interpD2_all{i,j}{1}(atWhatbin,:),50)

        colormap(gray(10))
        clim([-0.1 0.1])
        cb = colorbar;
        cb.Ticks=[-0.25:0.05:0.25];
        cb.Label.String = '$\mathcal{D}_{2}^{u}~[m^2/s^3]$';
        cb.Label.Interpreter = 'latex';
        cb.Label.FontSize = 14;
        cb.TickLabelInterpreter = 'latex';
       
        xticks([-1:0.5:1].*turbDia+turbHub(2))
        xlim([mask.y(1) mask.y(2)])
        yticks([0:0.25:1.25].*turbDia)
        ylim([0 max(mask.z(2))])
        yticklabels({'0','0.25','0.5','0.75','1','1.25'})
        xticklabels({'-1','-0.5','0','0.5','1'})
        ax=gca;
        ax.FontSize=14;
        ax.TickLabelInterpreter="latex";
        xlabel('$y/D[-]$',Interpreter='latex')
        ylabel('$z/D[-]$',Interpreter='latex')
        title(sprintf('$u_{\\infty,hub}=%.1f\\,\\mathrm{m/s}$', uinfhub(i)), ...
            'Interpreter', 'latex')
    end
end
exportgraphics(gcf,'images\D2_contour_midbin_yz5_diffOp.png')