function [D1_all, D2_all, centersALL] = ...
    func_estimateDriftDiffusion(timearr, dataArr, nbins, nComp)

% Function to estimate the first and second Kramers-Moyal coefficients 
% using conditional moments
% Author : N.S.Dangi ; TU Delft

if nargin < 3 || isempty(nbins)
    nbins = 10;
end
if nargin < 4 || isempty(nComp)
    nComp = numel(dataArr);
end

nComp = min(nComp, numel(dataArr));

dt = mean(diff(timearr));
Npts = size(dataArr{1}, 2);

D1_all = cell(1, nComp);
D2_all = cell(1, nComp);
centersALL = cell(1, nComp);

for c = 1:nComp

    fprintf('Estimating drift/diffusion for component %d...\n', c);

    data = dataArr{c};

    edges = linspace(min(data(:)), max(data(:)), nbins+1);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));

    centersALL{c} = centers;

    D1_comp = nan(nbins, Npts);
    D2_comp = nan(nbins, Npts);
    n_samples = nan(nbins, Npts);

    for j = 1:Npts

        u = data(:,j);

        if range(u) == 0
            D1_comp(:,j) = 0;
            D2_comp(:,j) = 0;
            continue;
        end

        du = diff(u);
        u_t = u(1:end-1);

        for i = 1:nbins

            if i < nbins
                idx = u_t >= edges(i) & u_t < edges(i+1);
            else
                idx = u_t >= edges(i) & u_t <= edges(i+1);
            end

            if sum(idx) > 30
                D1_comp(i,j) = mean(du(idx)) / dt;
                D2_comp(i,j) = var(du(idx)) / (2*dt);
                % D2_comp(i,j) = mean(du(idx).^2)/(2*dt);
                n_samples(i,j)=sum(idx);
            end
        end
        
    end

    % Optional smoothing of estimated coefficients
    for j = 1:Npts
        D1_comp(:,j) = smooth(D1_comp(:,j), 3, 'sgolay', 2);
        D2_comp(:,j) = smooth(D2_comp(:,j), 3, 'sgolay', 2);
    end

    D1_all{c} = D1_comp;
    D2_all{c} = D2_comp;

end

end
