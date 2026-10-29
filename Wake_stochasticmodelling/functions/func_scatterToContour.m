function func_scatterToContour(coords, values,levels)

    x = coords(:,2);
    y = coords(:,3);
    F = scatteredInterpolant(x,y,values(:),'natural','none');
    [X,Y] = meshgrid(linspace(min(x),max(x),100), ...
                     linspace(min(y),max(y),100));
    contourf(X,Y,F(X,Y),levels,'LineColor','none');
    colorbar; axis equal;
end
