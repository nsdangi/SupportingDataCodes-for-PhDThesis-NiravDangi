imgfolder = fullfile(pwd, 'images');
if ~exist(imgfolder, 'dir')
    mkdir(imgfolder);
end
%% TI LES

figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        
        func_scatterToContour(coords{i,j}, 100.*nanstd(windfield{i,j}{1},1)./nanmean(windfield{i,j}{1},1),10)

        colormap(gray(10))
        clim([3 20])
        cb = colorbar;
        cb.Ticks=[0:5:20];
        cb.Label.String = '$TI_{u} [\%]$';
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

exportgraphics(gcf,'images\TI_contour_yz5_diffOp_LES_oneseed.png')
%% TI Stoch

figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        stoch_TIu = 100.*nanstd(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,:),1)./nanmean(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,:),1);
        
        func_scatterToContour(coords{i,j}, ...
            stoch_TIu,10)

        colormap(gray(10))
        clim([3 20])
        cb = colorbar;
        cb.Ticks=[0:5:20];
        cb.Label.String = '$TI_{u} [\%]$';
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

exportgraphics(gcf,'images\TI_contour_yz5_diffOp_stoch_oneseed.png')

%% TI Diff

figure(Position=[241.6667  129.6667  948.0000  498.0000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        les_TIu = 100.*nanstd(windfield{i,j}{1},1)./nanmean(windfield{i,j}{1},1);
        stoch_TIu = 100.*nanstd(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,:),1)./nanmean(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,:),1);
        
        func_scatterToContour(coords{i,j}, ...
            stoch_TIu-les_TIu,10)

        colormap(gray(10))
        clim([-5 5])
        cb = colorbar;
        cb.Ticks=[-5:2.5:5];
        cb.Label.String = '$TI_{u,stoch} -TI_{u,LES} [\%]$';
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

exportgraphics(gcf,'images\TIdiff_contour_yz5_diffOp_stoch_oneseed.png')
