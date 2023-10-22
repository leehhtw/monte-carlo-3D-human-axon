function btab = btabgen(gdir,bval,t,root)
%BTABGEN create a b-table for simulation.
% btab = btabgen(gdir, bval, t, root) create the b-table file for diffusion
% simulation. The b-table is composed of the gradient directions gdir,
% b-value bval, and t = [inter-pulse duration, pulse width, echo time].
%
% A b=0 is added in the beginning of the file.

Nt = size(t,1);
Nb = numel(bval);
Ng = size(gdir,1);
btab = [];
for i = 1:Nt
    ti = t(i,:);
    btab = cat(1,btab,[0 0 0 0 ti]);
    for j = 1:Nb
        bvalj = bval(j);
        for k = 1:Ng
            gdirk = gdir(k,:);
            btab = cat(1,btab,[gdirk bvalj ti]);
        end
    end
end
if nargin>3
fid = fopen(fullfile(root),'w');
for i = 1:size(btab,1)
    fprintf(fid,'% .16f % .16f % .16f %.6f %.6f %.6f %.6f \n',btab(i,:));
end
fclose(fid);
end
end