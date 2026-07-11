clear
clc

INP_common.parentDir     = '..\Data\CompDiffs';
INP_common.folders       = {'f0p01to0p03','f0p01to0p06','f0p06to0p25','Original'};
INP_common.folderstr     = {'0 < f < 0.03 Hz','0 < f < 0.06 Hz','0.06 < f < 0.25 Hz','Original'};
INP_common.subfolders    = {'Below','Rated','Above'};
INP_common.files         = {'NABL.out','SABL.out','USABL.out','iea22mw_dtalm.0.out'};
INP_common.filenames     = strrep(INP_common.files,'.out','');
INP_common.filenames{4}  ='IEC_NTM_Scaled_TI';
addpath("..\Data\iea22mw_3p5_p1");

INP_common.bladefile     =readtable("..\Data\iea22mw_3p5_p1\IEA-22-280-RWT_AeroDyn15_blade.dat");
INP_common.span          =INP_common.bladefile.BlSpn;

INP_common.binEdges      = 0:5:360;

for mj=1:length(INP_common.folders)
    folder=INP_common.folders{mj};
    for kj=1:length(INP_common.subfolders)
        subfolder=INP_common.subfolders{kj};
        for lj=1:length(INP_common.files)
            file=INP_common.files{lj};
            filestr=INP_common.filenames{lj};
            folderPath = fullfile(INP_common.parentDir, folder,subfolder);
            OPT.(folder).(subfolder).(filestr).fastout = postpfast(file,folderPath);
        end
    end
end

fieldNames = {'AB1N0', 'AB1N0','AB1N0'};
fieldNamesALSO = {'AB1N00', 'AB1N00', 'AB1N00'};
fieldSuffix = {'Fx', 'Fy', 'Alpha'};
fieldToUpdate = {'Bin_AVG_Fx', 'Bin_AVG_Fy', 'Bin_AVG_AOA'};
fieldToUpdateSTD = {'Bin_AVG_FxSTD', 'Bin_AVG_FySTD', 'Bin_AVG_AOASTD'};
fieldBinavgstr={'F_{x}/F_{x} (IEC NTM Scaled TI)',...
                'F_{y}/F_{y} (IEC NTM Scaled TI)',...
                'AOA/AOA (IEC NTM Scaled TI)'};
%%
for mj=length(INP_common.folders)
    folder=INP_common.folders{mj};
    for kj = 1:length(INP_common.subfolders)
        subfolder = INP_common.subfolders{kj};
        for lj=1:length(INP_common.files)
            filestr=INP_common.filenames{lj};
            for fieldIdx = 1:numel(fieldNames)
                outputY = [];
                outputYstd=[];
                for jj = 1:length(INP_common.span)
                    inputX = OPT.(folder).(subfolder).(filestr).fastout.Azimuth;
                    if jj>=10
                        inputY = OPT.(folder).(subfolder).(filestr).fastout.([fieldNames{fieldIdx} num2str(jj) fieldSuffix{fieldIdx}]);
                    else
                        inputY = OPT.(folder).(subfolder).(filestr).fastout.([fieldNamesALSO{fieldIdx} num2str(jj) fieldSuffix{fieldIdx}]);
                    end

                    [outputYcomp, newtheta,outputYcompstd]=binAveraging(inputY,inputX,'binWidth',2);
                    outputYcompstd=outputYcompstd';
                    outputY=[outputY,outputYcomp];
                    outputYstd=[outputYstd,outputYcompstd];
                end
                OPT_BIN.(folder).(subfolder).(filestr).(fieldToUpdate{fieldIdx}) = outputY;
                OPT_BIN.(folder).(subfolder).(filestr).theta=newtheta;
                OPT_BIN.(folder).(subfolder).(filestr).(fieldToUpdateSTD{fieldIdx}) = outputYstd;
            end
            
        end
    end
end
%%
clearvars -except INP_common OPT
mean_podS=struct();
max_podS=struct();
std_podS=struct();
whatFields={'RtAeroPwr','RotThrust','RootMFlp1','TwrBsMyt','TwrBsMxt'};
for uj=1:length(whatFields)
    fieldHeremsm=whatFields{uj};
    for mj=1:length(INP_common.folders)
        folder=INP_common.folders{mj};
        for kj=1:length(INP_common.subfolders)
            subfolder=INP_common.subfolders{kj};
            for lj=1:length(INP_common.files)
                file=INP_common.files{lj};
                filestr=INP_common.filenames{lj};
                folderPath = fullfile(INP_common.parentDir, folder,subfolder);
                mean_podS.(folder).(subfolder).(filestr).(fieldHeremsm)=mean(OPT.(folder).(subfolder).(filestr).fastout.(fieldHeremsm));
                max_podS.(folder).(subfolder).(filestr).(fieldHeremsm)=max(OPT.(folder).(subfolder).(filestr).fastout.(fieldHeremsm));
                std_podS.(folder).(subfolder).(filestr).(fieldHeremsm)=std(OPT.(folder).(subfolder).(filestr).fastout.(fieldHeremsm));
            end
        end
    end
end
clearvars -except INP_common OPT *podS
%%
fields={'f0p01to0p03','f0p01to0p06','f0p06to0p25','Original'};
whatoperation={'Below','Rated','Above'};
whatFields={'RtAeroPwr','RotThrust','RootMFlp1','TwrBsMyt','TwrBsMxt'};
for mj=1:length(whatoperation)
    subfolder=whatoperation{mj};
    NABL_data = [];
    SABL_data = [];
    USABL_data = [];
    IEC_data = [];
    for i = 1:numel(fields)
        field_data = mean_podS.(fields{i}).(subfolder);
        NABL_data = [struct2array(field_data.NABL)'];
        SABL_data = [struct2array(field_data.SABL)'];
        USABL_data = [struct2array(field_data.USABL)'];
        IEC_data = [struct2array(field_data.IEC_NTM_Scaled_TI)'];
        T_mean.(subfolder).(fields{i}) = table(NABL_data, SABL_data, USABL_data, IEC_data,'RowNames',whatFields);

        field_data = max_podS.(fields{i}).(subfolder);
        NABL_data = [struct2array(field_data.NABL)'];
        SABL_data = [struct2array(field_data.SABL)'];
        USABL_data = [struct2array(field_data.USABL)'];
        IEC_data = [struct2array(field_data.IEC_NTM_Scaled_TI)'];
        T_max.(subfolder).(fields{i}) = table(NABL_data, SABL_data, USABL_data, IEC_data,'RowNames',whatFields);

        field_data = std_podS.(fields{i}).(subfolder);
        NABL_data = [struct2array(field_data.NABL)'];
        SABL_data = [struct2array(field_data.SABL)'];
        USABL_data = [struct2array(field_data.USABL)'];
        IEC_data = [struct2array(field_data.IEC_NTM_Scaled_TI)'];
        T_std.(subfolder).(fields{i}) = table(NABL_data, SABL_data, USABL_data, IEC_data,'RowNames',whatFields);
    end
    
end
clearvars -except INP_common OPT *podS T_*
%% NORM MEAN
fields = {'f0p01to0p03','f0p01to0p06','f0p06to0p25','Original'};
legendnames={'0 < f < 0.03 Hz','0 < f < 0.06 Hz','0.06 < f < 0.25 Hz','Original'};
colors=  {'k','m','b','r'};
markers= {'o','<','>','*'};
whatoperation={'Below','Rated','Above'};
whatste=strrep(strrep({'NABL.out','SABL.out','USABL.out','IEC_NTM_Scaled_TI.out'}, '.out', ''), '_', ' ');
whatvar={'RtAeroPwr','RotThrust','RootMFlp1','TwrBsMyt','TwrBsMxt'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','Q/ Q(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','TDx_{flap}/ TDx_{flap} (IEC NTM Scaled TI)'};
% ylabel_names={'\sigma_{P}/\sigma_{P} (IEC NTM Scaled TI)','\sigma_{T}/\sigma_{T} (IEC NTM Scaled TI)','\sigma_{Q}/\sigma_{Q} (IEC NTM Scaled TI)','\sigma_{M_{flap}}/\sigma_{M_{flap}} (IEC NTM Scaled TI)','\sigma_{TDx_{flap}}/\sigma_{TDx_{flap}} (IEC NTM Scaled TI)'};

ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI)','M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI)'};


% ylabel_names={'\sigma_{P}/\sigma_{P} (IEC NTM Scaled TI)','\sigma_{T}/\sigma_{T} (IEC NTM Scaled TI)','\sigma_{M_{flap}}/\sigma_{M_{flap}} (IEC NTM Scaled TI)','\sigma_{M_{T-fa}}/ \sigma_{M_{T-fa}} (IEC NTM Scaled TI)','\sigma_{M_{T-ss}}/ \sigma_{M_{T-ss}} (IEC NTM Scaled TI)'};

% ylabel_names={'Max(P/P(IEC NTM Scaled TI))','Max(T/ T(IEC NTM Scaled TI))','Max(M_{flap}/ M_{flap} (IEC NTM Scaled TI))','Max(M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI))','Max(M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI))'};

for jj=1:length(whatoperation)
    for lj=1:length(whatvar)
        figure()
        for mj=1:length(fields)
        max_Y=[];

        max_Y=[max_Y;table2array(T_mean.(whatoperation{jj}).(fields{mj})(lj,:))];
        max_Y_original=table2array(T_mean.(whatoperation{jj}).(fields{4})(lj,4));
        norm_max_Y= max_Y./ max_Y_original;
        plot(norm_max_Y,Marker=markers{mj},Color=colors{mj},LineWidth=2)
        hold on
        
        end
        ylabel(ylabel_names{lj})
        xticks([1,2,3,4])
        xlim([1 4])
        xticklabels(whatste)
        ax=gca;
        ax.FontSize=12;
        box on
        grid on
        ylim([0.8 1.2])
        yticks([0.8:0.05:1.2])
        % ylim([0 2]) %std
        % yticks([0:0.25:2]) %std
        legend(legendnames,Location="northwest")
        if exist(['.\images_Mean\'], 'dir') ~= 7
                    mkdir(['.\images_Mean\']);
        end
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        saveas(gcf,['.\images_Mean\',whatvar{lj},whatoperation{jj},'.pdf']);
        saveas(gcf,['.\images_Mean\',whatvar{lj},whatoperation{jj},'.png']);
        close all
        
    end

end

%% NORM STD

fields = {'f0p01to0p03','f0p01to0p06','f0p06to0p25','Original'};

legendnames={'$0\, Hz < f < 0.03\, Hz$','$0\, Hz < f < 0.06\, Hz$','$0.06\, Hz < f < 0.25\, Hz$','$Original$'};

colors=  gray(5);
markers= {'<','*','o','pentagram'};
whatoperation={'Below','Rated','Above'};
% whatste=strrep(strrep({'NABL.out','SABL.out','USABL.out','IEC_NTM_Scaled_TI.out'}, '.out', ''), '_', ' ');
whatste={'$NABL$','$SABL$','$USABL$','$IEC_{s}$'};
whatvar={'RtAeroPwr','RotThrust','RootMFlp1','TwrBsMyt','TwrBsMxt'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','Q/ Q(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','TDx_{flap}/ TDx_{flap} (IEC NTM Scaled TI)'};
% ylabel_names={'\sigma_{P}/\sigma_{P} (IEC NTM Scaled TI)','\sigma_{T}/\sigma_{T} (IEC NTM Scaled TI)','\sigma_{Q}/\sigma_{Q} (IEC NTM Scaled TI)','\sigma_{M_{flap}}/\sigma_{M_{flap}} (IEC NTM Scaled TI)','\sigma_{TDx_{flap}}/\sigma_{TDx_{flap}} (IEC NTM Scaled TI)'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI)','M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI)'};


ylabel_names={'$\sigma_{P}/\sigma_{P}\, (IEC_{s})$','$\sigma_{T}/\sigma_{T}\, (IEC_{s})$','$\sigma_{M_{flap}}/\sigma_{M_{flap}}\, (IEC_{s})$','$\sigma_{M_{T-fa}}/ \sigma_{M_{T-fa}}\, (IEC_{s})$','$\sigma_{M_{T-ss}}/ \sigma_{M_{T-ss}}\, (IEC_{s})$'};

% ylabel_names={'Max(P/P(IEC NTM Scaled TI))','Max(T/ T(IEC NTM Scaled TI))','Max(M_{flap}/ M_{flap} (IEC NTM Scaled TI))','Max(M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI))','Max(M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI))'};

for jj=1:length(whatoperation)
    for lj=1:length(whatvar)
        figure()
        for mj=1:length(fields)
        max_Y=[];

        max_Y=[max_Y;table2array(T_std.(whatoperation{jj}).(fields{mj})(lj,:))];
        max_Y_original=table2array(T_std.(whatoperation{jj}).(fields{4})(lj,4));
        norm_max_Y= max_Y./ max_Y_original;
        plot(norm_max_Y,Marker=markers{mj},Color=colors(mj,:),LineWidth=2)
        hold on
        
        end
        ylabel(ylabel_names{lj},Interpreter="latex")
        xticks([1,2,3,4])
        xlim([1 4])
        xticklabels(whatste)
        ax=gca;
        ax.FontSize=14;
        ax.LineWidth=3;
        ax.MinorGridLineWidth=1;
        ax.TickLabelInterpreter="latex";
        box on
        
        % ylim([0.8 1.2])
        % yticks([0.8:0.05:1.2])
        ylim([0 2]) %std
        yticks([0:0.25:2]) %std
        legend(legendnames,Location="northwest",Interpreter="latex")
        if exist(['.\images_std\'], 'dir') ~= 7
                    mkdir(['.\images_std\']);
        end

        exportgraphics(gcf,['.\images_std\',whatvar{lj},whatoperation{jj},'.png']);
        close all
        
    end

end
 
%% NORM MAX


fields = {'f0p01to0p03','f0p01to0p06','f0p06to0p25','Original'};

legendnames={'$0\, Hz < f < 0.03\, Hz$','$0\, Hz < f < 0.06\, Hz$','$0.06\, Hz < f < 0.25\, Hz$','$Original$'};

colors=  gray(5);
markers= {'<','*','o','pentagram'};
whatoperation={'Below','Rated','Above'};
% whatste=strrep(strrep({'NABL.out','SABL.out','USABL.out','IEC_NTM_Scaled_TI.out'}, '.out', ''), '_', ' ');
whatste={'$NABL$','$SABL$','$USABL$','$IEC_{s}$'};
whatvar={'RtAeroPwr','RotThrust','RootMFlp1','TwrBsMyt','TwrBsMxt'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','Q/ Q(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','TDx_{flap}/ TDx_{flap} (IEC NTM Scaled TI)'};
% ylabel_names={'\sigma_{P}/\sigma_{P} (IEC NTM Scaled TI)','\sigma_{T}/\sigma_{T} (IEC NTM Scaled TI)','\sigma_{Q}/\sigma_{Q} (IEC NTM Scaled TI)','\sigma_{M_{flap}}/\sigma_{M_{flap}} (IEC NTM Scaled TI)','\sigma_{TDx_{flap}}/\sigma_{TDx_{flap}} (IEC NTM Scaled TI)'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI)','M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI)'};

ylabel_names={'$max({P})/max({P})\, (IEC_{s})$','$max({T})/max({T})\, (IEC_{s})$','$max({M_{flap}})/max({M_{flap}})\, (IEC_{s})$','$max({M_{T-fa}})/ max({M_{T-fa}})\, (IEC_{s})$','$max({M_{T-ss}})/ max({M_{T-ss}})\, (IEC_{s})$'};

% ylabel_names={'$\sigma_{P}/\sigma_{P}\, (IEC_{s})$','$\sigma_{T}/\sigma_{T}\, (IEC_{s})$','$\sigma_{M_{flap}}/\sigma_{M_{flap}}\, (IEC_{s})$','$\sigma_{M_{T-fa}}/ \sigma_{M_{T-fa}}\, (IEC_{s})$','$\sigma_{M_{T-ss}}/ \sigma_{M_{T-ss}}\, (IEC_{s})$'};

% ylabel_names={'Max(P/P(IEC NTM Scaled TI))','Max(T/ T(IEC NTM Scaled TI))','Max(M_{flap}/ M_{flap} (IEC NTM Scaled TI))','Max(M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI))','Max(M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI))'};

for jj=1:length(whatoperation)
    for lj=1:length(whatvar)
        figure()
        for mj=1:length(fields)
        max_Y=[];

        max_Y=[max_Y;table2array(T_max.(whatoperation{jj}).(fields{mj})(lj,:))];
        max_Y_original=table2array(T_max.(whatoperation{jj}).(fields{4})(lj,4));
        norm_max_Y= max_Y./ max_Y_original;
        plot(norm_max_Y,Marker=markers{mj},Color=colors(mj,:),LineWidth=2)
        hold on
        
        end
        ylabel(ylabel_names{lj},Interpreter="latex")
        xticks([1,2,3,4])
        xlim([1 4])
        xticklabels(whatste)
        ax=gca;
        ax.FontSize=14;
        ax.LineWidth=3;
        ax.MinorGridLineWidth=1;
        ax.TickLabelInterpreter="latex";
        box on
        
        ylim([0.8 1.3])
        yticks([0.8:0.05:1.3])
        % ylim([0 2]) %std
        % yticks([0:0.25:2]) %std
        legend(legendnames,Location="northwest",Interpreter="latex")
        if exist(['.\images_max\'], 'dir') ~= 7
                    mkdir(['.\images_max\']);
        end

        exportgraphics(gcf,['.\images_max\',whatvar{lj},whatoperation{jj},'.png']);
        close all
        
    end

end
%% NOT NORMALISED
fields = {'f0p01to0p03','f0p01to0p06','f0p06to0p25','Original'};
legendnames={'0 < f < 0.03 Hz','0 < f < 0.06 Hz','0.06 < f < 0.25 Hz','Original'};
colors=  {'k','m','b','r'};
markers= {'o','<','>','*'};
whatoperation={'Below','Rated','Above'};
whatste=strrep(strrep({'NABL.out','SABL.out','USABL.out','IEC_NTM_Scaled_TI.out'}, '.out', ''), '_', ' ');
whatvar={'RtAeroPwr','RotThrust','RootMFlp1','TwrBsMyt','TwrBsMxt'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','Q/ Q(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','TDx_{flap}/ TDx_{flap} (IEC NTM Scaled TI)'};
% ylabel_names={'\sigma_{P}/\sigma_{P} (IEC NTM Scaled TI)','\sigma_{T}/\sigma_{T} (IEC NTM Scaled TI)','\sigma_{Q}/\sigma_{Q} (IEC NTM Scaled TI)','\sigma_{M_{flap}}/\sigma_{M_{flap}} (IEC NTM Scaled TI)','\sigma_{TDx_{flap}}/\sigma_{TDx_{flap}} (IEC NTM Scaled TI)'};

% ylabel_names={'P/P(IEC NTM Scaled TI)','T/ T(IEC NTM Scaled TI)','M_{flap}/ M_{flap} (IEC NTM Scaled TI)','M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI)','M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI)'};


ylabel_names={'\sigma_{P}/\sigma_{P} (IEC NTM Scaled TI)','\sigma_{T}/\sigma_{T} (IEC NTM Scaled TI)','\sigma_{M_{flap}}/\sigma_{M_{flap}} (IEC NTM Scaled TI)','\sigma_{M_{T-fa}}/ \sigma_{M_{T-fa}} (IEC NTM Scaled TI)','\sigma_{M_{T-ss}}/ \sigma_{M_{T-ss}} (IEC NTM Scaled TI)'};

% ylabel_names={'Max(P/P(IEC NTM Scaled TI))','Max(T/ T(IEC NTM Scaled TI))','Max(M_{flap}/ M_{flap} (IEC NTM Scaled TI))','Max(M_{T-fa}/ M_{T-fa} (IEC NTM Scaled TI))','Max(M_{T-ss}/ M_{T-ss} (IEC NTM Scaled TI))'};

for jj=1:length(whatoperation)
    for lj=1:length(whatvar)
        figure()
        for mj=1:length(fields)
        max_Y=[];

        max_Y=[max_Y;table2array(T_std.(whatoperation{jj}).(fields{mj})(lj,:))];
        max_Y_original=table2array(T_std.(whatoperation{jj}).(fields{4})(lj,4));
        norm_max_Y= max_Y./ 1;
        plot(norm_max_Y,Marker=markers{mj},Color=colors{mj},LineWidth=2)
        hold on
        
        end
        ylabel(ylabel_names{lj})
        xticks([1,2,3,4])
        xlim([1 4])
        xticklabels(whatste)
        ax=gca;
        ax.FontSize=12;
        box on
        grid on

        legend(legendnames,Location="northwest")
        if exist(['.\NOTNORMimages_STD\'], 'dir') ~= 7
                    mkdir(['.\NOTNORMimages_STD\']);
        end
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        saveas(gcf,['.\NOTNORMimages_STD\',whatvar{lj},whatoperation{jj},'.pdf']);
        saveas(gcf,['.\NOTNORMimages_STD\',whatvar{lj},whatoperation{jj},'.png']);
        close all
        
    end

end
%% bin avg ones
fieldToUpdate = {'Bin_AVG_Fx', 'Bin_AVG_Fy', 'Bin_AVG_AOA'};
fieldBinavgstr={'$F_{x}/F_{x} (IEC_{s})$',...
                '$F_{y}/F_{y} (IEC_{s})$',...
                '$AOA/AOA (IEC_{s})$'};
addpath("C:\Users\nsdangi\OneDrive - Delft University of Technology\syncingtheoldfiles\fromDesktop\somescripts");
maxcolors=[];
custom = customcolormap_preset('red-white-blue'); % choose color map scheme
close all
for uu=1%:length(fieldToUpdate)
    fieldHere=fieldToUpdate{uu};
    fieldstrhere=fieldBinavgstr{uu};
    
    for mj=length(INP_common.folders)
        folder=INP_common.folders{mj};
        folderstr=INP_common.folderstr{mj};
        for kj = 1:length(INP_common.subfolders)
            subfolder = INP_common.subfolders{kj};
            for lj=1:length(INP_common.files)
                filestr=INP_common.filenames{lj};
                thetahere=OPT_BIN.(folder).(subfolder).(filestr).theta;
                thetahere=linspace(0,360,length(thetahere))';
                figure(Visible="off")

                color_ax_numnorm=OPT_BIN.(folder).(subfolder).(filestr).(fieldHere);%Bin_AVG_Fx(:,i);
                color_ax_denomnorm=OPT_BIN.(INP_common.folders{4}).(subfolder).(INP_common.filenames{4}).(fieldHere);%Bin_AVG_Fx(:,i);
                try
                color_ax=color_ax_numnorm./color_ax_denomnorm;
                [polcontH,polcontC]=polarPcolor(ones(1,59).*transpose(INP_common.span)./142, ...
                    thetahere', ...
                    color_ax,'Nspokes',7,'Ncircles',5,'RtickLabel',{'$0$','$0.25$','$0.5$','$0.75$','$1$'});
                
                catch error
                    disp(['Hello I have encountered an issue for',folder,' ', subfolder, ' ', filestr])
                    color_ax=color_ax_numnorm./color_ax_denomnorm(1:height(color_ax_numnorm),:);
                    [polcontH,polcontC]=polarPcolor(ones(1,59).*transpose(INP_common.span)./142, ...
                    ones(1,length(INP_common.binEdges)-1).*INP_common.binEdges(1:end-1), ...
                    color_ax,'Nspokes',7,'Ncircles',5,'RtickLabel',{'$0$','$0.25$','$0.5$','$0.75$','$1$'});
                end

                clim([0.85 1.15])
                colorbar
                % title(strrep(subfolder,'_',' '))
                box on
                ax = gca;
                ax.LineWidth=3;
                ax.FontSize=14;
                ax.TickLabelInterpreter="latex";
                colormap(gray(7))
                cb = colorbar;
               

                
                cb=colorbar('Ticks', linspace(0.85, 1.15, 7));  
  
                cb.FontSize=14;
                % title([subfolder,' rated operation :: ',strrep(filestr,'_',' ')])
                % subtitle(folderstr)
                if exist([folder, '\images_BINAVGmean'], 'dir') ~= 7
                    mkdir([folder, '\images_BINAVGmean']);
                end
                 
                cb.Label.String = fieldstrhere;
                cb.Label.Interpreter = "latex";
                cb.TickLabelInterpreter = "latex";
                
                % keyboard

                exportgraphics(gcf,[folder,'\images_BINAVGmean\',subfolder,filestr,fieldHere,'.png']);
                close all
            end
        end
    end
    disp(uu)
end


%%
function [datafastout] = postpfast(filename,folderPath)
fullFilePath = fullfile(folderPath, filename);
datafastout=readtable(fullFilePath,"FileType","text","NumHeaderLines",6);
idxkeep=datafastout.Time>=min(datafastout.Time)+60 & datafastout.Time<=660;
datafastout=datafastout(idxkeep,:);
meantable=mean(datafastout);
meantable=removevars(meantable, "Time");
stdtable=std(datafastout);
stdtable=removevars(stdtable, "Time");


end