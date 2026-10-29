imgfolder = fullfile(pwd, 'images');
if ~exist(imgfolder, 'dir')
    mkdir(imgfolder);
end
%% MEAN LES

figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        
        func_scatterToContour(coords{i,j}, nanmean(windfield{i,j}{1},1)./uinfhub(i),10)

        colormap(gray(10))
        clim([0.3 1])
        cb = colorbar;
        cb.Ticks=[0:0.2:1];
        cb.Label.String = '$u/u_{\infty,hub} [-]$';
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

exportgraphics(gcf,'images\mean_contour_yz5_diffOp_LES_oneseed.png')
%% MEAN Stoch

figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        stoch_u = nanmean(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,:),1)./uinfhub(i);
        
        func_scatterToContour(coords{i,j}, ...
            stoch_u,10)

        colormap(gray(10))
        clim([0.3 1])
        cb = colorbar;
        cb.Ticks=[0:0.2:1];
        cb.Label.String = '$u/u_{\infty,hub} [-]$';
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

exportgraphics(gcf,'images\mean_contour_yz5_diffOp_stoch_oneseed.png')

%% MEAN Diff

figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        stoch_u_abs = nanmean(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,:),1);
        les_u_abs   = nanmean(windfield{i,j}{1},1);
        func_scatterToContour(coords{i,j}, ...
            stoch_u_abs-les_u_abs,10)

        colormap(gray(10))
        clim([-2 2])
        cb = colorbar;
        cb.Ticks=[-2:0.5:2];
        cb.Label.String = '$u_{stoch}-u_{LES} [m/s]$';
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

exportgraphics(gcf,'images\meanDiff_contour_yz5_diffOp_stoch_oneseed.png')
