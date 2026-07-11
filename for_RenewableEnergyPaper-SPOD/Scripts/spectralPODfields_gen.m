clear
clc
load('..\Data\3Dvelarray.mat_3Dup.mat')
addpath("spod_matlab-1.1.2")
%%
close all
writeOPT=1;

f_lowpass       = [0.03,0.06,0.25];
f_highpass      = [0.01,0.01,0.06];
fnamesave={'f0p01to0p03','f0p01to0p06','f0p06tof0p25'};

Nblock = 4;
for jj=1%:size(INP_common.subfolders,2)
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
        OG_Vel_U=reshape(X1_u', whatsizetimedim,numY, numZ);



        numY=height(unique(points.Y));
        numZ=height(unique(points.Z));


        recons_Vel_U=reshape(reconstructed_field_U', whatsizetimedim,numY, numZ);
        recons_Vel_V=reshape(reconstructed_field_V', whatsizetimedim,numY, numZ);
        recons_Vel_W=reshape(reconstructed_field_W', whatsizetimedim,numY, numZ);

        velocity = zeros(whatsizetimedim, 3, numY, numZ);
        velocity(:,1,:,:)=recons_Vel_U;
        velocity(:,2,:,:)=recons_Vel_V;
        velocity(:,3,:,:)=recons_Vel_W;



        dy=max(unique(diff(points.Y)));
        dz=max(unique(diff(points.Z)));
        zHub=170;
        z1=min(points.Z);

        gridHt=max(points.Z)-min(points.Z);

        zOffset=z1+gridHt/2;
        z0=0.0001;
        % summary file {zHub, Clockwise, UBAR, TI_u, TI_v, TI_w}
        Clockwise=1;
        UBAR=mean(recons_Vel_U(:,18,18));
        TI_u=100*(std(recons_Vel_U(:,18,18))/mean(recons_Vel_U(:,18,18)));
        TI_v=100*(std(recons_Vel_V(:,18,18))/mean(recons_Vel_U(:,18,18)));
        TI_w=100*(std(recons_Vel_W(:,18,18))/mean(recons_Vel_U(:,18,18)));
        SummVars=[zHub; Clockwise; UBAR; TI_u; TI_v; TI_w];

        pathsaving = ['./', char(subfolder), '/', fnamesave{mm}];

        filename = fullfile(pathsaving, 'forum_code_wnd_v2');

        if ~exist(pathsaving, 'dir')
            mkdir(pathsaving);
        end

        if writeOPT==1
            writeBLgrid(filename, velocity, dy, dz, dt, zOffset, z0, SummVars)
            write_wind_file(dt, UBAR, TI_u, TI_v, TI_w,pathsaving)
            mkdir([pathsaving,'/',strrep(num2str(sum(energyused)*100),'.','_'),'perc'])
        end

        

    end
    clearvars -except INP_common OPT_3dVel points sorted_indices timesteps
end

%%
function writeBLgrid(FileName, velocity, dy, dz, dt, zOffset, z0, SummVars)


% writeBLgrid('testWind', velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars)
%
% writes wind velocity data to binary .wnd files
%

% inputs - FileName: Name of output .wnd-file (extension will be added)
%        - velocity: 4D-array: (time, 3D-windcomp, y, z)
%        - dz
%        - dy
%        - dt
%        - SumVars: 6 variables from the summary file {zHub, Clockwise, UBAR, TI_u, TI_v, TI_w}
%                                                        90          1    12    12   9.6     6
%        - zOffset: (=90) Reference height (m) = Z(1) + GridHeight / 2.0
%        - z0 = .03;      Roughness length (m)

%-----------------------------------------
% INITIALIZE VARIABLES
%-----------------------------------------

fileFmt  = 'int16';

fc = 4;                 % should be 4 to allow turbulence intensity to be stored in the header
lat = 0;                % latitude (deg)

zHub = SummVars(1);     % hub height [m]
MFFWS = SummVars(3);    % mean full-field wind speed
TI_U = SummVars(4);     % Turbulence Intensity of u component (%)
TI_V = SummVars(5);     % Turbulence Intensity of v component (%)
TI_W = SummVars(6);     % Turbulence Intensity of w component (%)

[nt, nffc, ny, nz] = size(velocity); % determin dimensions of windfield, e.g. [1286,3,23,23]

z1 =5;  %this is the bottom of the grid

dx = dt*MFFWS;          % delta x in m
nt_header = floor(nt/2);% half the number of time steps


%-----------------------------------------
% OPEN FILE
%-----------------------------------------
len    = length(FileName);
ending = FileName(len-3:len);

if strcmpi( ending, '.wnd' )
    FileName = FileName(1:len-4);
end
fid_wnd   = fopen( [ FileName '.wnd' ], 'w' ); % Open file, or create new file, for writing; discard existing contents, if any.
if ( fid_wnd <= 0 )
    error( 'Wind file could not be opened.' );
    return;
end


%-----------------------------------------
%WRITE THE HEADER OF THE BINARY FILE
%-----------------------------------------
% THE NEWER-STYLE AERODYN WIND FILE

fwrite( fid_wnd, -99, 'int16' );             % number of components
fwrite( fid_wnd, fc, 'int16' );              % should be 4 to allow turbulence intensity to be stored in the header
fwrite( fid_wnd, nffc, 'int32' );            % number of components (should be 3)
fwrite( fid_wnd, lat, 'float32' );           % latitude (deg)
fwrite( fid_wnd, z0, 'float32' );            % Roughness length (m)
fwrite( fid_wnd, zOffset, 'float32' );       % Reference height (m) = Z(1) + GridHeight / 2.0
fwrite( fid_wnd, TI_U, 'float32' );          % Turbulence Intensity of u component (%)
fwrite( fid_wnd, TI_V, 'float32' );          % Turbulence Intensity of v component (%)
fwrite( fid_wnd, TI_W, 'float32' );          % Turbulence Intensity of w component (%)

fwrite( fid_wnd, dz, 'float32' );            % delta z in m
fwrite( fid_wnd, dy, 'float32' );            % delta y in m
fwrite( fid_wnd, dx, 'float32' );            % delta x in m
fwrite( fid_wnd, nt_header, 'int32' );       % half the number of time steps
fwrite( fid_wnd, MFFWS, 'float32');          % mean full-field wind speed

fwrite( fid_wnd, [0 0 0], 'float32' );    % unused variables (for BLADED): write zeros
fwrite( fid_wnd, [0 0], 'int32' );        % unused variables (for BLADED): write zeros
fwrite( fid_wnd, nz, 'int32' );              % number of points in vertical direction
fwrite( fid_wnd, ny, 'int32' );              % number of points in horizontal direction
fwrite( fid_wnd, zeros(3*(nffc-1),1), 'int32' );     % unused variables (for BLADED): write zeros
%-----------------------------------------
% WRITE GRID DATA
%-----------------------------------------
Scale    = 0.00001*SummVars(3)*SummVars(4:6);
Offset   = [SummVars(3) 0 0];


if SummVars(2) > 0 %clockwise rotation
    %flip the y direction....
    %let's change the dimension of velocity so that it's 4-d instead of 3-d
    y_ix = ny:-1:1;
else
    y_ix = 1:ny;
end

v = zeros(nz*ny*nffc,1);

for it = 1:nt
    cnt = 1;
    for iz = 1:nz
        for iy = y_ix
            for k=1:nffc
                v(cnt) = (velocity(it,k,iy,iz) - Offset(k))/Scale(k);
                cnt = cnt + 1;
                %velocity(it,k,iy,iz) = v(cnt2)*Scale(k) + Offset(k);
            end %for k
        end %iy
    end % iz

    fwrite( fid_wnd, v, fileFmt );

end %it

%-----------------------------------------
% CLOSE .WND FILE
%-----------------------------------------
fclose(fid_wnd);



% y    = [0:ny-1]*dy - dy*(ny-1)/2;
% zHub = SummVars(1);
% z1   = zHub - ZGoffset - dz*(nz-1)/2;  %this is the bottom of the grid
% z    = [0:nz-1]*dz + z1;
end
%%
function write_wind_file(dt, UBar, TI_u, TI_v, TI_w, pathsaving)
% Open file for writing

fileID = fopen(fullfile(pathsaving, 'forum_code_wnd_v2.sum'), 'w');

% Write header
fprintf(fileID, 'T  Clockwise\n');
fprintf(fileID, '0 -  NONE  IEC turbulence models scaled to exact specified standard deviation\n\n');
fprintf(fileID, 'Turbine/Model Specifications:\n');
fprintf(fileID, '         35  Vertical grid-point matrix dimension\n');
fprintf(fileID, '         35  Horizontal grid-point matrix dimension\n');
fprintf(fileID, '         %.2f  Time step [seconds]\n', dt);
fprintf(fileID, '         700  Analysis time [seconds]\n');
fprintf(fileID, '         700  Usable output time [seconds]\n\n');

% Write hub height
fprintf(fileID, '170.000  Hub height [m]\n\n');

% Write UBar and TI values
fprintf(fileID, 'UBar  =  %.4f m/s\n', UBar);
fprintf(fileID, 'TI(u) =  %.4f %%\n', TI_u);
fprintf(fileID, 'TI(v) =  %.4f %%\n', TI_v);
fprintf(fileID, 'TI(w) =  %.4f %%\n', TI_w);

% Close file
fclose(fileID);
end

