imgfolder = fullfile(pwd, 'images');
if ~exist(imgfolder, 'dir')
    mkdir(imgfolder);
end
%%
colorchoices = gray(3);
figure(Position=[51	171.666666666667	1223.33333333333	472.000000000000])
for i = 1:numel(cases)
    nexttile
    for j = 1%1:numel(sliceNames)

        timearr = outtime{i,j}-min(outtime{i,j});
        dtscal  = median(diff(timearr));
        langArr = timearr(1) + (0:3*numel(timearr)-1)*dtscal;
        Fs = 1/dtscal;
        Nblock = 6;
        Ncoh = round(numel(timearr) / Nblock);

        % indexes at 50% spanwise location and hub height +- 50m

        [~,idxHub] = min((coords{i,j}(:,2)-(turbHub(2)+0.5*turbDia)).^2 + ...
            (coords{i,j}(:,3)-turbHub(3)).^2);

        [~,idxHub_altp50] = min((coords{i,j}(:,2)-(turbHub(2)+0.5*turbDia)).^2 + ...
            (coords{i,j}(:,3)-(turbHub(3)+50)).^2);

        [~,idxHub_altm50] = min((coords{i,j}(:,2)-(turbHub(2)+0.5*turbDia)).^2 + ...
            (coords{i,j}(:,3)-(turbHub(3)-50)).^2);


        y_les   = detrend(windfield{i,j}{1}(:,idxHub));
        y_stoch = detrend(u_syn_ALL{i,j}{1}(end-round(size(u_syn_ALL{i,j}{1},1)/n_transients):end,idxHub));

        [freq_les, ft_les] = func_spectrum(y_les, Ncoh, Fs);
        [freq_stoch, ft_stoch] = func_spectrum(y_stoch, Ncoh, Fs);



        loglog(freq_les, ft_les,LineWidth=2,Color=colorchoices(1,:))
        hold on
        loglog(freq_stoch, ft_stoch,LineWidth=2,Color=colorchoices(2,:))
        loglog(freq_les,freq_les.^(-5/3),Color=colorchoices(1,:),LineWidth=2,LineStyle="--")
        ax=gca;
        ax.FontSize=14;
        ax.TickLabelInterpreter="latex";
       

        xlabel('$f [Hz]$',Interpreter='latex')
        ylabel('$S_{xx}(f)\;[\mathrm{(m/s)^2/Hz}]$', 'Interpreter', 'latex');
        grid on
        ax.XMinorGrid='off';
        ax.YMinorGrid="off";
        ax.LineWidth=3;

        box on
        axis tight
        xlim([1e-2 1e0])
        xticks([1e-3 1e-2 1e-1 1e0])
        yticks([1e-6 1e-4 1e-2 1e0 1e1 1e2 1e3])
        title(sprintf('$u_{\\infty,hub}=%.1f\\,\\mathrm{m/s}$', uinfhub(i)), ...
            'Interpreter', 'latex')





    end

end
lg = legend('$LES$','$Stochastic$','$f^{5/3}$','Interpreter', 'latex',Location='northoutside',Orientation='horizontal');
lg.Position= [0.3796    0.4905    0.2419    0.0529];
exportgraphics(gcf,'images/les_vs_stoch_ucomp_hub_spec.png')