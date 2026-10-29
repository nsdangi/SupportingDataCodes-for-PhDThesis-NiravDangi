function [D1_all, D2_all] = func_interpolateDriftDiffusion( ...
    D1_all, D2_all, centersALL)

% Supplementing function to interpolate missing first and second 
% Kramers-Moyal coefficients 
% Author : N.S.Dangi ; TU Delft

for c = 1:numel(D1_all)

    D1 = D1_all{c};
    D2 = D2_all{c};
    centers = centersALL{c};

    for j = 1:size(D1,2)

        valid = isfinite(D1(:,j)) & isfinite(D2(:,j));
        missing = ~valid;

        if sum(valid) >= 2

            d1 = D1(:,j);
            d2 = D2(:,j);

            d1(~valid) = NaN;
            d2(~valid) = NaN;

            
            D1(:,j) = fillmissing(d1, 'linear', ...
                'SamplePoints', centers, ...
                'EndValues', 'nearest');

            D2(:,j) = fillmissing(d2, 'linear', ...
                'SamplePoints', centers, ...
                'EndValues', 'nearest');

        elseif sum(valid) == 1

            
            D1(missing,j) = D1(valid,j);
            D2(missing,j) = D2(valid,j);

        end
        
    end



    D1_all{c} = D1;
    D2_all{c} = D2;
end

end