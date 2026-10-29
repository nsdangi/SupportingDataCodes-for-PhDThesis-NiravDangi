function [freqs, ft] = func_spectrum(ll, n, sampl_freq, spec_only)
% Author : N.S.Dangi ; TU Delft
if nargin < 4
    spec_only = false;
end

% Ensure the number of samples is an integer multiple of n
num_samples = floor(length(ll) / n) * n;
ll_trimmed = ll(1:num_samples); % Trim the signal
% disp(['You trimmed: ', num2str(length(ll) - length(ll_trimmed))]);

% Reshape into n rows
lreal = reshape(ll_trimmed, n, []);

% Compute power spectrum
ft = sum(abs(fft(lreal, [], 1)).^2, 2) / (size(lreal, 1) * sampl_freq);

% Output results
if spec_only
    freqs = [];
    ft = ft(1:floor(n/2)+1);
else
    freqs = (0:floor(n/2))' * (sampl_freq / n); % Ensure integer indices
    ft = ft(1:floor(n/2)+1);
end
end