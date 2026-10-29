function u_syn_ALL = func_sim_Langevin( ...
    timearr, dataArr, D1, D2, centers, nComp, n_transients)

% Function to siulative Langevin trajectories of underlying PDF;
% uses the first and second Kramers-Moyal coefficients 

% example n_transients = 3; means 3 times longer than original series is generated
% This transient enables mean to be found from the randn start point

% Author : N.S.Dangi ; TU Delft

if nargin < 6 || isempty(nComp)
    nComp = numel(dataArr);
end

nComp = min([nComp, numel(dataArr), ...
    numel(D1), numel(D2), numel(centers)]);

Nt = length(timearr);
dt = mean(diff(timearr));

u_syn_ALL = cell(1, nComp);

for c = 1:nComp

    fprintf('Generating Langevin trajectories for component %d...\n', c);

    data = dataArr{c};
    D1_comp = D1{c};
    D2_comp = D2{c};
    centers_comp = centers{c}(:);

    Npts = size(data, 2);

    dataMin = min(data(:), [], 'omitnan');
    dataMax = max(data(:), [], 'omitnan');

    % Generate Langevin trajectories
    u_syn = zeros(n_transients * Nt, Npts);

    for j = 1:Npts

        valid = isfinite(D1_comp(:,j)) & ...
            isfinite(D2_comp(:,j)) & ...
            isfinite(centers_comp);

        validCenters = centers_comp(valid);
        d1_valid = D1_comp(valid,j);
        d2_valid = D2_comp(valid,j);

        % Initial condition
        u_syn(1,j) = randn;

        % -------------------------------------------------
        % Euler-Maruyama integration
        % -------------------------------------------------
        for t = 1:n_transients * Nt - 1

            xcurr = u_syn(t,j);

            % Estimate drift and diffusion at current state
            if numel(validCenters) >= 2

                if xcurr < validCenters(1) || xcurr > validCenters(end)

                    % Use nearest valid coefficient bin outside the range
                    [~, idx_near] = min(abs(validCenters - xcurr));

                    drift = d1_valid(idx_near);
                    diffu = d2_valid(idx_near);

                else

                    % Linear interpolation inside the valid range
                    drift = interp1( ...
                        validCenters, d1_valid, xcurr, 'linear');

                    diffu = interp1( ...
                        validCenters, d2_valid, xcurr, 'linear');

                end

            elseif numel(validCenters) == 1

                drift = d1_valid(1);
                diffu = d2_valid(1);

            else

                drift = 0;
                diffu = 0;

            end

            % Ensure nonnegative diffusion
            diffu = max(diffu, 0);

            % Euler-Maruyama update
            u_next = xcurr + drift * dt ...
                + sqrt(2 * diffu * dt) * randn;

            % Optional bound by observed data range
            u_next = min(max(u_next, dataMin), dataMax);

            u_syn(t+1,j) = u_next;

        end

        % Progress update every 5%
        progress = 100 * j / Npts;

        if j == 1 || j == Npts || ...
                floor(progress/5) > floor(100*(j-1)/Npts/5)

            fprintf('Component %d/%d: %.1f%% complete\n', ...
                c, nComp, progress);
        end

    end

    u_syn_ALL{c} = u_syn;

end

end
