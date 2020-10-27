% * * * * * * * * * * * * * * * * * * * * *
%
% Peaks phase distribution
%
% * * * * * * * * * * * * * * * * * * * * *

function [phi_wrapped] = PeaksPhase(n, kappa, snr, debug, output)

% Dimension
N = 2^n;

fprintf('Dimension %d\n',N)

% amplitude default = 5

if kappa <= 0
    kappa1 = 5;
else 
    kappa1 = kappa;
end


% Compute continuous phase

phi = kappa1 * peaks(N);


% Plot continuous phase for debug

if debug==1
    phi_plot = normalize(phi);
    figure, imagesc(phi_plot)
    title('Peaks phase distribution')
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
    title('Wrapped Peaks phase')
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
