imgfolder = fullfile(pwd, 'images');
if ~exist(imgfolder, 'dir')
    mkdir(imgfolder);
end
colorchoices=gray(3);

%%
figure(Position=[51	171.666666666667	1223.33333333333	472.000000000000])
for i = 1:numel(cases)
        nexttile
    for j = 1%1:numel(sliceNames)
        
        timearr = outtime{i,j}-min(outtime{i,j});
        dtscal  = median(diff(timearr));
        langArr = timearr(1) + (0:3*numel(timearr)-1)*dtscal;

        % indexes at 50% spanwise location and hub height +- 50m

        [~,idxHub] = min((coords{i,j}(:,2)-(turbHub(2)+0.5*turbDia)).^2 + ...
              (coords{i,j}(:,3)-turbHub(3)).^2);

        [~,idxHub_altp50] = min((coords{i,j}(:,2)-(turbHub(2)+0.5*turbDia)).^2 + ...
              (coords{i,j}(:,3)-(turbHub(3)+50)).^2);

        [~,idxHub_altm50] = min((coords{i,j}(:,2)-(turbHub(2)+0.5*turbDia)).^2 + ...
              (coords{i,j}(:,3)-(turbHub(3)-50)).^2);
        
        plot(timearr,windfield{i,j}{1}(:,idxHub),LineWidth=2,Color=colorchoices(1,:))
        hold on
        plot(langArr,u_syn_ALL{i,j}{1}(:,idxHub),LineWidth=2,Color=colorchoices(2,:))
        ax=gca;
        ax.FontSize=14;
        ax.TickLabelInterpreter="latex";
        xticks(0:300:max(langArr))
        yticks([0:3:30])
        xlabel('$Time [s]$',Interpreter="latex");
        ylabel('$u [m/s]$',Interpreter="latex");

        title(sprintf('$u_{\\infty,hub}=%.1f\\,\\mathrm{m/s}$', uinfhub(i)), ...
            'Interpreter', 'latex')
        
        
        
        

    end

end
lg = legend('$LES$','$Stochastic$','Interpreter', 'latex',Location='northoutside',Orientation='horizontal');
lg.Position= [0.4109    0.4876    0.1769    0.0512]; 
exportgraphics(gcf,'images/les_vs_stoch_ucomp_hub.png')