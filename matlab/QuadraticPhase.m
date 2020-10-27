% * * * * * * * * * * * * * * * * * * * * *
%
% Quadratic phase distribution
%
% * * * * * * * * * * * * * * * * * * * * *

function [phi_wrapped]=QuadraticPhase(n, kappa, snr, debug, output)

% Dimension
N = 2^n;

% amplitude default = 0.001

if kappa <= 0
    kappa1 = 0.001;
else 
    kappa1 = kappa;
end


% Compute continuous phase

[posy, posx] = meshgrid(1:N, 1:N);
phi = kappa1 * ((posy - N/2).^2 + (posx - N/2).^2);


% Plot continuous phase for debug

if debug==1
    phi_plot = normalize(phi);
    figure, imagesc(phi_plot)
    title('Quadratic phase distribution')
    colormap(gray)
    set(gca,'FontSize',13)
end


% Complex signal
phi = exp(1i*phi);

% Add noise if it is required

phi = awgn(phi, snr, 'measured'); 

% Wrapped phase

phi_wrapped = atan2(imag(phi), real(phi));

% Plot wrapped phase for debug
if debug==1
    phi_plot = normalize(phi_wrapped);
    figure, imagesc(phi_plot)
    title('Wrapped Quadratic phase')
    colormap(gray)
    set(gca,'FontSize',13)
end


% Save wrapped phase distribution into a binary file

fprintf('Saving into file: %s\n',output);

fid = fopen(output,'w');

if fid==-1
    fprintf('Cannot create output file %s\n',output);
else
    fwrite(fid, phi_wrapped,'float');
    fclose(fid);
end
