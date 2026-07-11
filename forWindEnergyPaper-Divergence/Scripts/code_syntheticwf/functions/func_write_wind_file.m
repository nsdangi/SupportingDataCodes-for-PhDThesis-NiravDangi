function func_write_wind_file(dt, UBar, TI_u, TI_v, TI_w,filename,nyy,nzz,timemax,hubheight)
    % Open file for writing
    % fileID = fopen('forum_code_wnd_v2.sum', 'w');
    
    % pathsaving = ['./', char(subfolder)];

    
    fileID = fopen([filename,'.sum'], 'w');


    % Write header
    fprintf(fileID, 'T  Clockwise\n');
    fprintf(fileID, '0 -  NONE  IEC turbulence models scaled to exact specified standard deviation\n\n');
    fprintf(fileID, 'Turbine/Model Specifications:\n');
    fprintf(fileID, '         %.0f  Vertical grid-point matrix dimension\n',nyy);
    fprintf(fileID, '         %.0f  Horizontal grid-point matrix dimension\n',nzz);
    fprintf(fileID, '         %.2f  Time step [seconds]\n', dt);
    fprintf(fileID, '         %.0f  Analysis time [seconds]\n',timemax);
    fprintf(fileID, '         %.0f  Usable output time [seconds]\n\n',timemax);

    % Write hub height
    fprintf(fileID, '%2f  Hub height [m]\n\n',hubheight);

    % Write UBar and TI values
    fprintf(fileID, 'UBar  =  %.4f m/s\n', UBar);
    fprintf(fileID, 'TI(u) =  %.4f %%\n', TI_u);
    fprintf(fileID, 'TI(v) =  %.4f %%\n', TI_v);
    fprintf(fileID, 'TI(w) =  %.4f %%\n', TI_w);

    % Close file
    fclose(fileID);
end

