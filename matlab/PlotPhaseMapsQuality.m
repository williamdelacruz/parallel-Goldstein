close all
clc


% * * * * * * * * * * * * * * * * * *
%  Phase, correlation, mask,
%  unwrapped, residues, branch-cuts
% * * * * * * * * * * * * * * * * * *


filename = 'ifsar.512x512';
ysize = 512;
xsize = 512;
type = 'float';
type1 = 'char';
 

% - - - - - - - - - - - - - -
%        Phase map
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.phase'));
if fid ~= -1
    data = fread(fid, type);
    fclose(fid);

    data = reshape(data, ysize, xsize);
    data1 = data;
    data = data/max(max(data));
    figure(1), 
     subplot(2,2,1), 
    imagesc(data), colormap(gray)
    title('Phase map')
    %set(gca,'FontSize',15)
else
    fprintf('File %s does not exists.\n', strcat(filename,'.phase'))
end


% - - - - - - - - - - - - - -
%     Residues and cuts
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.brc'));
if fid ~= -1
    data = fread(fid, 'char');
    fclose(fid);

    data =reshape(data, ysize, xsize);
    data = data/max(max(data));
    figure(1), subplot(2,2,2), imagesc(data), colormap(gray)
    title('Residues and cuts')
else
    fprintf('File %s does not exists.\n', strcat(filename,'.qual'))
end


% - - - - - - - - - - - - - -
%        Path integration
% - - - - - - - - - - - - - -

fid = fopen(strcat(filename,'.path'));
if fid ~= -1
    data = fread(fid,'int');
    fclose(fid);

    data = reshape(data, ysize, xsize);
    data = data/max(max(data));
    figure(1), subplot(2,2,3), imagesc(data), colormap(gray)
    title('Path integration')
else
    fprintf('File %s does not exists.\n', strcat(filename,'.path'))
end


% - - - - - - - - - - - - - -
%    Unwrapped phase map
% - - - - - - - - - - - - - -

fid = fopen(filename);
if fid ~= -1
    data = fread(fid,'float');
    fclose(fid);

    %data = data - min(min(data));
    data = reshape(data, ysize, xsize);
    data = data/max(max(data));
    figure(1), subplot(2,2,4), imagesc(data), colormap(gray)
    title('Unwrapped phase map')
else
    fprintf('File %s does not exists.\n', filename)
end
