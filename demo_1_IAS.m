% This example demonstrates the diffusion simualtions in realistic axons
% from a human brain EM sample.
%
% Author: Hong-Hsi Lee, October, 2023 (orcid.org/0000-0002-3663-6559)

clear
restoredefaultpath
filePath = matlab.desktop.editor.getActiveFilename;
root = fileparts(filePath);
addpath(genpath(fullfile(root,'lib')));

%% Visualize 3D intra-axonal space (IAS)
root_ias = fullfile(root,'data','ias_align');
files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

figure('unit','inch','position',[0 0 10 5]);
tic;
j = 0;
for i = sort(randsample(nfiles,10).')
    j = j+1;
    
    % load axonal shape
    filei = load(fullfile(root_ias,files(i).name));
    ias   = filei.ias_align;
    
    % Create 3D triangulation
    FV = isosurface(ias,0);
    TR = triangulation(FV.faces, FV.vertices);
    
    % Plot 3D triangulation
    subplot(2,5,j)
    trisurf(TR,'edgealpha',0);
    axis equal off
    material dull
    camlight
    title(sprintf('# %u',i),'FontWeight','normal');
end
toc;

%% Save the axonal shape into the bin file for simulation
root_ias   = fullfile(root,'data','ias_align');
root_input = fullfile(root,'data','ias_input');
files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

RMS = rmsobj();
for i = 1:nfiles
    % load axonal shape
    filei = load(fullfile(root_ias,files(i).name));
    ias   = filei.ias_align;
    
    % Length of each pixel, nm
    pix_l = 64;
    
    % Save in bin file
    filename = sprintf('%06u.bin',i);
    RMS.saveSubstrate(fullfile(root_input,filename), ias, pix_l);
end

%% Design protocol
root_input   = fullfile(root,'data','ias_input');
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
root_input  = fullfile(root,'data','ias_input');
root_output = fullfile(root,'data','ias_output');
mkdir(root_output);

files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% Total time to perform simulation
Tmax = max(TE);

% Number of random walkers
Npar = 1e5;

% Create shell scrip to run myrms
fileID = fopen(fullfile(root_output,'job.sh'),'w');
fprintf(fileID,'#!/bin/bash\n');
for i = 1:nfiles
    filename = sprintf('%06u.bin',i);
    filename_pgse = sprintf('%06u_pgse.bin',i);
    fprintf(fileID,'%s/myrms %s %s %s -pgse %s -time %u -particle %u -space 1 \n',...
        fullfile(root,'lib','rms'),...
        fullfile(root_input,filename),...
        fullfile(root_input,'btable.txt'),...
        fullfile(root_output,filename),...
        fullfile(root_output,filename_pgse),...
        Tmax,...
        Npar );
end
fclose(fileID);

% Please open the folder ./data/ias_output and run the shell script in terminal: sh ./myrms

%% Read simulation result
root_ias    = fullfile(root,'data','ias_align');
root_input  = fullfile(root,'data','ias_input');
root_output = fullfile(root,'data','ias_output');

files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

% number of b-values
Nb  = 8;

% Spherical mean signal
S_sm = zeros(Nb, nfiles);

RMS = rmsobj();
for i = 1:nfiles
    filename      = sprintf('%06u.bin',i);
    filename_pgse = sprintf('%06u_pgse.bin',i);
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
    S_sm(:, i) = mean(Si, 2);
end

%% Plot figure

% Model fitting
b = c;
delta = 10;
Delta = 20;
D0 = 2;
smt = AxCaliberSMT1(b, delta, Delta, D0);
model = 'Neuman';

r = zeros(nfiles,1);
beta = zeros(nfiles,1);
for i = 1:nfiles
    [r(i), beta(i)] = smt.getAxonRadius(S_sm(:,i), model);
end

i = 1;
% Fitting curve
b_fit = 1./linspace(0.01,1,100).^2;
smt_fit = AxCaliberSMT1(b_fit, delta, Delta, D0);
S_fit = smt_fit.AxonDiameterFWD([r(i), beta(i)], model);

% Plot signal vs 1/sqrt(b)
figure;
hold on;
plot(1./sqrt(b), S_sm(:,i), 'o');
plot(1./sqrt(b_fit), S_fit, '-');
xlabel('$1/\sqrt{b}$','interpreter','latex','fontsize',20);
ylabel('$S/S_0$','interpreter','latex','fontsize',20);
xlim([0 1]);
ylim([0 1]);
box on; grid on; pbaspect([1 1 1]);
text(0.1, 0.9, sprintf('$r_{\\rm MR}$=%.2f $\\mu$m', r(i)),'interpreter','latex',...
    'fontsize',20);

%% Calculate axon radius
root_ias = fullfile(root,'data','ias_align');
files    = dir(fullfile(root_ias,'*.mat'));
nfiles   = numel(files);

r_eff  = zeros(nfiles, 1);   % effective radius, um
CV_r   = zeros(nfiles, 1);   % coefficient of variationa of radius
CV_A   = zeros(nfiles, 1);   % coefficient of variation of cross-sectional area
w0     = zeros(nfiles, 1);   % undulation amplitude, um
lambda = zeros(nfiles, 1);   % undulation wavelength, um
RMS = rmsobj();
for i = 1:nfiles
    % load axonal shape
    filei = load(fullfile(root_ias,files(i).name));
    ias   = filei.ias_align;
    
    % length of each pixel, nm
    pix_l = 64e-3;
    
    % cross-sectional area
    Az = squeeze(sum(sum(ias,1),2)) * pix_l^2;
    
    % coefficient of variation of cross-sectional area
    CV_A(i)  = std(Az)/mean(Az);
    
    % equivalent circle radius
    rz = sqrt(Az/pi);
    
    % coefficient of variation of radius
    CV_r(i)  = std(rz)/mean(rz);
    
    % effective radius
    r_eff(i) = ( sum(rz.^6)/sum(rz.^2) )^(1/4);
    
    % axonal skeleton
    cm = RMS.centermass(ias, pix_l);

    % undulation amplitude and wavelength
    [w0(i), lambda(i)] = RMS.undulation(cm);
    
end

figure;
plot(r_eff, r, 'o');
xlabel('effective radius','interpreter','latex','fontsize',20);
ylabel('fitted radius','interpreter','latex','fontsize',20);
xlim([0 2]);
ylim([0 2]);
hr = refline(1,0); set(hr,'color','k');
box on; grid on; pbaspect([1 1 1]);


