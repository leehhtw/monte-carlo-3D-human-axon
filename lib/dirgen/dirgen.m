function gdir = dirgen(ndir)
%DIRGEN create the gradient direction in ndir directions using mrtrix.

filePath = matlab.desktop.editor.getActiveFilename;
root = fileparts(filePath);

system(sprintf('dirgen %u %s -cartesian -force',...
    ndir,...
    fullfile(root,'tmp.txt') ));
fid = fopen(fullfile(root,'tmp.txt'),'r');
fgetl(fid);
gdir = zeros(ndir,3);
for i = 1:ndir
    tline = fgetl(fid);
    gdir(i,:) = str2num(tline);
end
fclose(fid);

delete(fullfile(root,'tmp.txt'));
end