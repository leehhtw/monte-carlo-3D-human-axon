% This example demonstrates the diffusion simualtions in synthetic fibers
% of scaled undulations and caliber variations.
%
% Author: Hong-Hsi Lee, October, 2023 (orcid.org/0000-0002-3663-6559)

clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root = fileparts(filePath);
addpath(genpath(fullfile(root,'lib')));

%% Create artificial cylinder with scaled undulation and caliber variation
root_ias = fullfile(root,'data','ias_align');
root_input = fullfile(root,'data','cyl_input'); mkdir(root_input);
files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% scaling for undulation amplitude
un = [0 0.25 0.50 0.75 1];

% scaling for caliber variation, coefficient of vairation of radius
cv = [0 0.25 0.50 0.75 1];

tic;
for i = 1:nfiles
    for j = 1:numel(un)
        for k = 1:numel(cv)
            uni = un(j);
            cvi = cv(k);
            
            % load axonal shape
            filename = sprintf('%06u.mat',i);
            filei = load(fullfile(root_ias,filename));
            ias = filei.ias_align;
            
            % scale the undulation and caliber variation of the real axon, 
            % and translate the scaled shape into a fiber of circular cross
            % sections
            RMS = rmsobj();
            ias_cyl = RMS.bwaxon2cyl(ias, uni, cvi);
            
            % Length of each pixel, nm
            pix_l = 64;
            
            filename = sprintf('%06u_un%u_cv%u.bin', i, uni*100, cvi*100);
            RMS.saveSubstrate(fullfile(root_input, filename), ias_cyl, pix_l);
        end
    end
end
toc;


%% Visualize 3D cylinder
root_ias = fullfile(root,'data','ias_align');
root_input = fullfile(root,'data','cyl_input');
files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% scaling for undulation amplitude
un = [0 0.25 0.50 0.75 1];

% scaling for caliber variation, coefficient of vairation of radius
cv = [0 0.25 0.50 0.75 1];

figure('unit','inch','position',[0 0 10 10]);
tic;
RMS = rmsobj();
for i = 1:nfiles
    for j = 1:numel(un)
        for k = 1:numel(cv)
            uni = un(j);
            cvi = cv(k);
            
            % load cylindrical shape
            filename = sprintf('%06u_un%u_cv%u.bin', i, uni*100, cvi*100);
            cyl = RMS.readSubstrate(fullfile(root_input,filename));
            
            % crop and down sample to accelerate the plotting
            cyl = cyl(:,:,floor(end/3):floor(end*2/3));
            cyl = cyl(1:2:end,1:2:end,1:2:end);
            
            % Plot 3D triangulation
            subplot(5, 5, (j-1)*numel(cv) + k)
            h = RMS.plotaxon(cyl);
            set(h,'edgealpha', 0, 'facecolor', [1 1 1]*0.6);
            axis equal off;
            camlight;
            material dull;
            title(sprintf('p_w %u%%, p_r %u%%',uni*100, cvi*100),...
                'FontWeight', 'normal');
        end
    end
end
toc;

%% Design protocol
root_input   = fullfile(root,'data','cyl_input');
ndir = 60;
gdir = dirgen(ndir);
bval = [1, 2, 3, 5, 7, 12, 17, 26];
DEL  = 20;
del  = 10;
TE   = DEL + del;

btab = btabgen(gdir, bval, [DEL, del, TE], fullfile(root_input, 'btable.txt'));

%% Compile CUDA C++ code
fid = fopen(fullfile(root,'lib','rms','compile_rms.sh'),'w');
fprintf(fid,sprintf('nvcc %s -o %s',...
    fullfile(root,'lib','rms','main_rms.cu'),...
    fullfile(root,'lib','rms','myrms') ));
fclose(fid);

% Please open the folder ./lib/rms and run the shell script in terminal: sh ./compile_rms.sh

%% Perform simulation
root_ias    = fullfile(root,'data','ias_align');
root_input  = fullfile(root,'data','cyl_input');
root_output = fullfile(root,'data','cyl_output');
mkdir(root_output);

files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% scaling for undulation amplitude
un = [0 0.25 0.50 0.75 1];

% scaling for caliber variation, coefficient of vairation of radius
cv = [0 0.25 0.50 0.75 1];

% Total time to perform simulation
Tmax = max(TE);

% Number of random walkers
Npar = 1e5;

% Create shell scrip to run myrms
fileID = fopen(fullfile(root_output,'job.sh'),'w');
fprintf(fileID,'#!/bin/bash\n');
for i = 1:nfiles
    for j = 1:numel(un)
        for k = 1:numel(cv)
            uni = un(j);
            cvi = cv(k);
            
            filename = sprintf('%06u_un%u_cv%u.bin', i, uni*100, cvi*100);
            filename_pgse = sprintf('%06u_un%u_cv%u_pgse.bin', i, uni*100, cvi*100);
            fprintf(fileID,'%s/myrms %s %s %s -pgse %s -time %u -particle %u -space 1 \n',...
                fullfile(root,'lib','rms'),...
                fullfile(root_input,filename),...
                fullfile(root_input,'btable.txt'),...
                fullfile(root_output,filename),...
                fullfile(root_output,filename_pgse),...
                Tmax,...
                Npar );
        end
    end
end
fclose(fileID);

% Please open the folder ./data/ias_output and run the shell script in terminal: sh ./myrms

%% Read simulation result
root_ias    = fullfile(root,'data','ias_align');
root_input  = fullfile(root,'data','cyl_input');
root_output = fullfile(root,'data','cyl_output');

files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% scaling for undulation amplitude
un = [0 0.25 0.50 0.75 1];

% scaling for caliber variation, coefficient of vairation of radius
cv = [0 0.25 0.50 0.75 1];

% number of b-values
Nb  = 8;

% Spherical mean signal
S_sm = zeros(Nb, nfiles, numel(un), numel(cv));

RMS = rmsobj();
for i = 1:nfiles
    for j = 1:numel(un)
        for k = 1:numel(cv)
            uni = un(j);
            cvi = cv(k);
            filename      = sprintf('%06u_un%u_cv%u.bin', i, uni*100, cvi*100);
            filename_pgse = sprintf('%06u_un%u_cv%u_pgse.bin', i, uni*100, cvi*100);
            filename_btab = 'btable.txt';
            RMS.readPGSE(fullfile(root_output,filename),...
                fullfile(root_input, filename_btab),...
                fullfile(root_output,filename_pgse));

            t = [RMS.DEL RMS.del RMS.TE];
            [C,IA,IC] = unique(t,'rows');
            Si = [];
            for ii = 1:numel(IA)
                ti = C(ii,:);
                TDi = ti(1); Tdi = ti(2);
                Td(ii) = Tdi;
                TD(ii) = TDi;
                listi = IC==ii;
                bvali = RMS.bval(listi);
                bveci = RMS.bvec(listi,:);
                sigi  = RMS.pgse(listi);

                sig0  = sigi(1);
                sigi  = sigi(2:end)/sig0;
                bvali = bvali(2:end);
                bveci = bveci(2:end,:);
                [c,ia,ic] = unique(bvali,'rows');
                Ng    = nnz(ic==1);
                Sj = zeros(numel(ia),Ng);
                for jj = 1:numel(ia)
                    listj = ic==jj;
                    Sj(jj,:) = sigi(listj);
                end
                Si = cat(2, Si, Sj);
            end
            
            % Spherical mean signal
            S_sm(:, i, j, k) = mean(Si,2);
        end
    end
end

%% Plot figure

% Model fitting
b = c;
delta = 10;
Delta = 20;
D0 = 2;
smt = AxCaliberSMT1(b, delta, Delta, D0);
model = 'Neuman';

r    = zeros(nfiles, numel(un), numel(cv));
beta = zeros(nfiles, numel(un), numel(cv));
for i = 1:nfiles
    for j = 1:numel(un)
        for k = 1:numel(cv)
            [r(i, j, k), beta(i, j, k)] = smt.getAxonRadius(S_sm(:, i, j, k), model);
        end
    end
end

i = 1;
j = 5;
k = 5;
% Fitting curve
b_fit = 1./linspace(0.01,1,100).^2;
smt_fit = AxCaliberSMT1(b_fit, delta, Delta, D0);
S_fit = smt_fit.AxonDiameterFWD([r(i,j,k), beta(i,j,k)], model);

% Plot signal vs 1/sqrt(b)
figure;
hold on;
plot(1./sqrt(b), S_sm(:,i,j,k), 'o');
plot(1./sqrt(b_fit), S_fit, '-');
xlabel('$1/\sqrt{b}$','interpreter','latex','fontsize',20);
ylabel('$S/S_0$','interpreter','latex','fontsize',20);
xlim([0 1]);
ylim([0 1]);
box on; grid on; pbaspect([1 1 1]);
text(0.1, 0.9, sprintf('$r_{\\rm MR}$=%.2f $\\mu$m', r(i,j,k)),'interpreter','latex',...
    'fontsize',20);

%% Calculate fiber radius
root_input  = fullfile(root,'data','cyl_input');
files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% scaling for undulation amplitude
un = [0 0.25 0.50 0.75 1];

% scaling for caliber variation, coefficient of vairation of radius
cv = [0 0.25 0.50 0.75 1];

r_eff  = zeros(nfiles, numel(un), numel(cv));   % effective radius, um
CV_r   = zeros(nfiles, numel(un), numel(cv));   % coefficient of variationa of radius
CV_A   = zeros(nfiles, numel(un), numel(cv));   % coefficient of variation of cross-sectional area
w0     = zeros(nfiles, numel(un), numel(cv));   % undulation amplitude, um
lambda = zeros(nfiles, numel(un), numel(cv));   % undulation wavelength, um
for i = 1%:nfiles
    for j = 1:numel(un)
        for k = 1:numel(cv)
            uni = un(j);
            cvi = cv(k);
            filename = sprintf('%06u_un%u_cv%u.bin', i, uni*100, cvi*100);
            
            % load axonal shape
            cyl = RMS.readSubstrate(fullfile(root_input,filename));

            % length of each pixel, nm
            pix_l = 64e-3;

            % cross-sectional area
            Az = squeeze(sum(sum(cyl,1),2)) * pix_l^2;
            
            % coefficient of variation of cross-sectional area
            CV_A(i, j, k)  = std(Az)/mean(Az);

            % equivalent circle radius
            rz = sqrt(Az/pi);
            
            % coefficient of variation of radius
            CV_r(i, j, k)  = std(rz)/mean(rz);

            % effective radius
            r_eff(i, j, k) = ( sum(rz.^6)/sum(rz.^2) )^(1/4);
            
            % axonal skeleton
            cm = RMS.centermass(cyl, pix_l);
            
            % undulation amplitude and wavelength
            [w0(i, j, k), lambda(i, j, k)] = RMS.undulation(cm);
        end
    end
end

j = 5; k = 5;
figure;
plot(r_eff(:,j,k), r(:,j,k), 'o');
xlabel('effective radius','interpreter','latex','fontsize',20);
ylabel('fitted radius','interpreter','latex','fontsize',20);
xlim([0 2]);
ylim([0 2]);
hr = refline(1,0); set(hr,'color','k');
box on; grid on; pbaspect([1 1 1]);


