clear
clc
load('..\Data\3Dvelarray.mat_3Dup.mat')
addpath("spod_matlab-1.1.2")
%%
close all

f_lowpass       = [0.03,0.06,0.25];
f_highpass      = [0.01,0.01,0.06];
legendlabels={'$Original$','$0.01\, Hz < f < 0.03\, Hz$','$0.01\, Hz < f < 0.06\, Hz$','$0.06\, Hz < f < 0.25\, Hz$'};
legendlabelsSS={'$Original$','$0.01\, Hz < f < 0.03\, Hz$','$0.01\, Hz < f < 0.06\, Hz$','$0.06\, Hz < f < 0.25\, Hz$',''};
fnamesave={'f0p01to0p03','f0p01to0p06','f0p06tof0p25'};

Nblock = 4;
markers= {'pentagram','<','*','o'};
colors=(gray(5));

for jj=1
    subfolder=INP_common.subfolders{jj};
    disp(subfolder)

    X1_u=OPT_3dVel.(subfolder).vel_arr_U(:,:)';
    X1_v=OPT_3dVel.(subfolder).vel_arr_V(:,:)';
    X1_w=OPT_3dVel.(subfolder).vel_arr_W(:,:)';
    zarr=unique(points.Z);
    if jj == 4 || jj == 5 || jj == 6
        X1_u=X1_u(:,1:2100);
        X1_v=X1_v(:,1:2100,:);
        X1_w=X1_w(:,1:2100,:);
    end
    X1_u_avg = mean (X1_u ,2) ;

    X1_v_avg = mean (X1_v ,2) ;

    X1_w_avg = mean (X1_w ,2) ;
    if jj == 4 || jj == 5 || jj == 6
        dt=0.33;
        whatsizetimedim=2100;
        WINDOW=300;
        NOVERLAP=round(WINDOW/2);
    else
        dt=0.25;
        whatsizetimedim=2800;
        WINDOW=400;
        NOVERLAP=round(WINDOW/2);
    end

    [L_u,P_u,F_u,~,A_u] = spod(transpose(X1_u-X1_u_avg),WINDOW,[],NOVERLAP,dt);
    [L_v,P_v,F_v,~,A_v] = spod(transpose(X1_v-X1_v_avg),WINDOW,[],NOVERLAP,dt);
    [L_w,P_w,F_w,~,A_w] = spod(transpose(X1_w-X1_w_avg),WINDOW,[],NOVERLAP,dt);

    for mm=1:length(f_lowpass)
        A_u_here = A_u;
        A_v_here = A_v;
        A_w_here = A_w;

        A_u_here(F_u>=f_lowpass(mm) | F_u<=f_highpass(mm),:,:) = 0;
        A_v_here(F_v>=f_lowpass(mm) | F_v<=f_highpass(mm),:,:) = 0;
        A_w_here(F_w>=f_lowpass(mm) | F_w<=f_highpass(mm),:,:) = 0;

        Lu_here = L_u;
        Lu_here(F_u>=f_lowpass(mm) | F_u<=f_highpass(mm),:,:) = 0;

        p_rec_u = invspod(P_u,A_u_here,rectwin(WINDOW),NOVERLAP);
        p_rec_v = invspod(P_v,A_v_here,rectwin(WINDOW),NOVERLAP);
        p_rec_w = invspod(P_w,A_w_here,rectwin(WINDOW),NOVERLAP);

        totalenergy=sum(L_u(:));
        energyhere=sum(Lu_here(:))/totalenergy;
        energyused=energyhere;
        disp([num2str(sum(energyused)*100),'%'])

        savinenergy(mm)=sum(energyused)*100;
       
        reconstructed_field_U = p_rec_u' + X1_u_avg;
        reconstructed_field_V = p_rec_v' + X1_v_avg;
        reconstructed_field_W = p_rec_w' + X1_w_avg;

        numY=height(unique(points.Y));
        numZ=height(unique(points.Z));
        recons_Vel_U=reshape(reconstructed_field_U', whatsizetimedim,numY, numZ);
        recons_Vel_V=reshape(reconstructed_field_V', whatsizetimedim,numY, numZ);
        recons_Vel_W=reshape(reconstructed_field_W', whatsizetimedim,numY, numZ);
        OG_Vel_U=reshape(X1_u', whatsizetimedim,numY, numZ);

        
        if mm==1
            figure(1)
            plot(timesteps-min(timesteps),X1_u(613,:),'Color',colors(1,:),'LineWidth',2,'Marker',markers{1}, 'MarkerIndices', 1:round(length(timesteps)/10):length(timesteps))
            hold on

            figure(7)
            plot(timesteps-min(timesteps),X1_u(613,:),'Color',colors(1,:),'LineWidth',2,'Marker',markers{1}, 'MarkerIndices', 1:30:length(timesteps))
            hold on

            figure(2)
            uarruseOG=X1_u(:,1400);
            uarr_plotOG=reshape(uarruseOG,35,[]);
            plot(uarr_plotOG(17,:),zarr,'Color',colors(1,:),'LineWidth',2,Marker=markers{1},MarkerIndices= 1:4:length(zarr))
            hold on

            figure(3)
            if mm==1

                Puu=[];
                Fu=[];

                for hh=1:width(OG_Vel_U)
                       Ncoh = round(size(recons_Vel_U(:,18,18),1)/Nblock);
                       Fs = 1/dt;
                    [Puu(:,hh),Fu(:,hh)] = pwelch(detrend(OG_Vel_U(:,hh,18)), ...
                        Ncoh, ...
                        round(Ncoh/2), ...
                        Ncoh, ...
                        Fs);

                end

                puuavg = mean(Puu,2);

                funew = logspace(log10(0.005),log10(1),50);
                puuavgnew = interp1(Fu(:,1),puuavg,funew);

                semilogx(funew,funew.*puuavgnew, ...
                    'Color',colors(1,:), ...
                    'LineWidth',2, ...
                    'Marker',markers{1},'MarkerIndices', 1:5:length(funew))

                hold on

            end

        end

        
            figure(1)
            plot(timesteps-min(timesteps),reconstructed_field_U(613,:),'Color',colors(mm+1,:),'LineWidth',2,'Marker',markers{mm},'MarkerIndices', 1:round(length(timesteps)/10):length(timesteps))
            hold on

            figure(7)
            plot(timesteps-min(timesteps),reconstructed_field_U(613,:),'Color',colors(mm+1,:),'LineWidth',2,'Marker',markers{mm},'MarkerIndices', 1:30:length(timesteps))
            hold on


        
            figure(2)
            zarr=unique(points.Z);
            uarruse=reconstructed_field_U(:,1400);
            uarr_plot=reshape(uarruse,35,[]);
            plot(uarr_plot(17,:),zarr,'Color',colors(mm+1,:),'LineWidth',2,Marker=markers{mm+1},MarkerIndices= 1:4:length(zarr))
            hold on
        

        figure(3)
        Puu=[];
        Fu=[];
        Fs = 1/dt;
        Nblock = 4;
        Ncoh = round(size(recons_Vel_U(:,18,18),1)/Nblock);

        for hh=1:width(recons_Vel_U(:,:,18))
            [Puu(:,hh),Fu(:,hh)] = pwelch(detrend(recons_Vel_U(:,hh,18)), ...
                Ncoh, ...
                round(Ncoh/2), ...
                Ncoh, ...
                Fs);
        end
        puuavg=mean(Puu,2);


        funew=logspace(log10(0.005),log10(1),50);
        puuavgnew=interp1(Fu(:,1),puuavg,funew);
   
        semilogx(funew,funew.*puuavgnew, ...
            'Color',colors(mm+1,:), ...
            'LineWidth',2, ...
            'Marker',markers{mm+1},'MarkerIndices', 1:5:length(funew))
        hold on

    end
    figure(1)
    legend(legendlabels,Interpreter="latex",Location="southeast")
    xlim([0 700])
    ylim([12 16])
    ax=gca;
    grid off
    box on
    ax.LineWidth=3;
    ax.FontSize=16;
    ax.TickLabelInterpreter="latex";
    ylabel('$u\, [m/s]$',Interpreter="latex")
    xlabel('$t\,[s]$',Interpreter="latex")
    %exportgraphics(gca,'showspodinstsignal.png')


    figure(7)
    legend(legendlabels,Interpreter="latex",Location="southeast")
    xlim([0 100])
    ylim([12 16])
    ax=gca;
    grid off
    box on
    ax.LineWidth=3;
    ax.FontSize=16;
    ax.TickLabelInterpreter="latex";
    ylabel('$u\, [m/s]$',Interpreter="latex")
    xlabel('$t\,[s]$',Interpreter="latex")
    %exportgraphics(gca,'showspodinstsignal_zoomed.png')
    figure(2)

    ax=gca;
    grid off
    box on
    ax.LineWidth=3;
    ax.FontSize=16;
    ax.TickLabelInterpreter="latex";
    xlabel('$u\, [m/s]$',Interpreter="latex")
    ylabel('$z/D\, [-] $',Interpreter="latex")
    ylim([0 350])
    yticks([0:0.25:1.5].*280);
    yticklabels({'0','0.25','0.5','0.75','1','1.25','1.5'})
    yline(170,'k--','LineWidth',2)
    legend(legendlabelsSS,Interpreter="latex",Location="northwest")
    %exportgraphics(gca,'showspodshearsignal.png')


    figure(3)
    legend(legendlabels,Interpreter="latex",Location="northeast")
    ax=gca;
    grid off
    box on
    ax.LineWidth=3;
    ax.FontSize=16;
    ax.TickLabelInterpreter="latex";
    set(gca,'XScale','log','YScale','linear')
    yticks([0:0.02:0.1])
    ylim([0 0.1])
    xlim([6*1e-3 1e0])
    grid on
    grid minor
    ax.MinorGridLineWidth=1;
    xticks([1e-3 1e-2 1e-1 1e0])
    ylabel('$<S(f)*f>_{hub} [m^{2}/s^{2}]$',Interpreter='latex')
    xlabel('$ f [Hz]$',Interpreter='latex')
    %exportgraphics(gca,'showspodspecsignal.png')
    
    
    clearvars -except INP_common OPT_3dVel points sorted_indices timesteps
end
