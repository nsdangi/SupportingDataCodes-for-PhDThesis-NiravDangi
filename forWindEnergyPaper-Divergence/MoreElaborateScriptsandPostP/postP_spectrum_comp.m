%% Spectrum comparison
close all
% comparing with "NOT" constrained one.. you can notice the high
% frequency loss.. this is allevaited by the constrained algoritm as shown 
% in the article
for b = 1:length(basePaths)

    figure('Visible','on')
    hold on

    legend_entries = {};

    for c = 1:2

        puu_seed_all = [];

        for sss = 1:numSeeds

            dataog   = windData.(cases{1}).(basePaths{b}){sss};
            datahere = windData.(cases{c}).(basePaths{b}){sss};

            Fs = 1/dt(b);
            Nblock = 4;
            Ncoh = round(size(dataog,1)/Nblock);

            Nz = size(datahere,3);   % spanwise dimension

            puu_span_all = [];

            % ==========================================
            % Loop over spanwise locations
            % ==========================================
            for zz = 18%1:35

                signal = detrend(squeeze(datahere(:,1,zz,18)));

                [Puu,Fu] = pwelch(signal, ...
                                  Ncoh, ...
                                  round(Ncoh/2), ...
                                  Ncoh, ...
                                  Fs);

                puu_span_all = [puu_span_all, Puu];

            end

            % ==========================================
            % Average across spanwise direction FIRST
            % ==========================================
            puu_span_avg = mean(puu_span_all,2);

            % Store for seed statistics
            puu_seed_all = [puu_seed_all, puu_span_avg];

        end

        % ==========================================
        % Compute min–max across seeds
        % ==========================================
        puu_min  = min(puu_seed_all,[],2);
        puu_max  = max(puu_seed_all,[],2);
        puu_mean = mean(puu_seed_all,2);  % optional mean line

        freq = Fu(2:end);

        ylow  = puu_min(2:end) .* freq;
        yhigh = puu_max(2:end) .* freq;
        ymean = puu_mean(2:end) .* freq;

        % Avoid negative values for log scale
        ylow(ylow <= 0) = min(ymean(ymean>0))*1e-6;

        % ==========================================
        % Plot min–max band
        % ==========================================
        fill([freq; flipud(freq)], ...
             [ylow; flipud(yhigh)], ...
             colors(c,:), ...
             'FaceAlpha',0.25, ...
             'EdgeColor','none');

        % Mean line
        loglog(freq, ymean, ...
               'LineWidth',3, ...
               'Color',colors(c,:));

    end

    % ==========================================
    % Axis formatting
    % ==========================================
    set(gca,'XScale','log','YScale','log')
    xlim([1e-3 1e0])
    xticks([1e-3 1e-2 1e-1 1e0])
    yticks([10e-5 10e-4 10e-3 10e-2 10e-1 10e0 10e1])
    ylim([10e-6 10e-1])
    ylabel('$<S(f)*f>_{hub} (m^{2}/s^{2})$',Interpreter='latex')
    xlabel('$ f (Hz)$',Interpreter='latex')

    box on
    grid on
    % grid minor

    ax = gca;
    ax.FontSize = 16;
    ax.TickLabelInterpreter = "latex";
    ax.LineWidth = 3;
    ax.MinorGridLineWidth=1;
    
    legend({'','$U_{o}$','','$U_{div-free}$'},'Location','northeast','Interpreter','latex')

    % ==========================================
    % Export
    % ==========================================
    % keyboard
    exportgraphics(gcf, ...
        ['Images\USpec_', ...
        cases{1},'_',basePaths{b},'.png'])

end
