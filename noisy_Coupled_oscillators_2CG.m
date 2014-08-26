%% Coupled oscillators with Granule dendrites

% Wutup!

clear all

noisy_InitNetwork_2CG

%% Simulate LFP activity using V
save_path = '/Users/boleslawosinski/Documents/MATLAB/Kay/OB_PC_model/pics/2014/Jan 18';
app = 'PyrGrapt1';

numtrials = 1;

% ALWAYS CHECK tsim AND tfinal !!!!!!!
% VLFP
% [OSN Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] ...
%     = noisy_VLFP_2CG(numtp, numtrials, input_file, Delay, Wfrac, Hd);
% MitLFPs = MitLFPs(:,1);
% ILFP

[Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] ...
    = noisy_ILFP_2CG(numtp, numtrials, input_file, Delay, Wfrac, Hd);
MitLFPs = MitLFPs.GradistMit(:,1);

% Raster plots
mitH = PlotRasterplot(Mitral,param);
graproxH = PlotRasterplot(GraProximal,param);
pyrH = PlotRasterplot(Pyramidal,param);
fbaH = PlotRasterplot(Feedback,param);

mode = 1;
conn = 1;
comment = 'near critical point';
fname = '2CGfull_vary-dGCgmax_4';
writeparams(input_file,mode,comment,conn,fname)


% spike count profile
scount_mit = zeros(nmit,1);
for ii = 1:nmit
    scount_mit(ii) = sum(Mitral{ii}.S);
end   
bar(scount_mit,'k')
xlim([0 nmit])
ylim([0 150])
set(gca,'fontsize',14)
xlabel('Cell');ylabel('Counts')



% save rasters
cd(save_path)
saveas(mitH,['Raster_mit_',app],'tiff')
saveas(graproxH,['Raster_graprox_',app],'tiff')
saveas(pyrH,['Raster_pyr_',app],'tiff') 
saveas(fbaH,['Raster_fba_',app],'tiff') 


% Spike Coherence Plots
binsize = 2;
[Mit_Synch Mit_PSTH times] = SpikeSynch(Mitral,param,binsize);
[Pyr_Synch Pyr_PSTH times] = SpikeSynch(Pyramidal,param,binsize);
[Gra_Synch Gra_PSTH times] = SpikeSynch(GraProximal,param,binsize);
[Fba_Synch Fba_PSTH times] = SpikeSynch(Feedback,param,binsize);

figure(1)
subplot(2,1,1)
bar(times,Mit_PSTH,'facecolor',[0 0 0.5])
set(gca,'fontsize',14)
title(['Mit Coherence = ',num2str(roundn(Mit_Synch,-3)),',  Total N _{spikes} = ', num2str(sum(Mit_PSTH)),' '])
xlabel('Time (ms) ');ylabel('Spike Count ')
xlim([0 tsim])
subplot(2,1,2)
bar(times,Pyr_PSTH,'facecolor',[0 0.5 0])
set(gca,'fontsize',14)
title(['Pyr Coherence = ',num2str(roundn(Pyr_Synch,-3)),',  Total N _{spikes} = ', num2str(sum(Pyr_PSTH)),' '])
xlabel('Time (ms) ');ylabel('Spike Count ')
xlim([0 tsim])

figure(2)
subplot(2,1,1)
bar(times,Gra_PSTH,'facecolor',[0 0 0.9])
set(gca,'fontsize',14)
title(['Gra Coherence = ',num2str(roundn(Gra_Synch,-3)),',  Total N _{spikes} = ', num2str(sum(Gra_PSTH)),' '])
xlabel('Time (ms) ');ylabel('Spike Count ')
xlim([0 tsim])
subplot(2,1,2)
bar(times,Fba_PSTH,'facecolor',[0 0.8 0])
set(gca,'fontsize',14)
title(['Fba Coherence = ',num2str(roundn(Fba_Synch,-3)),',  Total N _{spikes} = ', num2str(sum(Fba_PSTH)),' '])
xlabel('Time (ms) ');ylabel('Spike Count ')
xlim([0 tsim])




% plot ISI and PSTH
twin = 10; % time window for binning (ms)
ISImit = [];
PSTHmit = zeros(1,dt*numtp/twin);
for ii = 1:nmit
    spiketimes = dt*find(Mitral{ii}.S == 1);
    ISImit = [ISImit diff(spiketimes)];
    PSTHmit = PSTHmit + histc(spiketimes,1:twin:tsim);
end   
figure(121)
% subplot(1,2,1)
% hist(ISImit,200)
% set(gca,'fontsize',14)
% xlabel('ISI (ms)');ylabel('Counts')
% title('Mitral')
% subplot(1,2,2)
% bar(1:twin:tsim,PSTHmit,'k')
% xlim([0 tsim])
% set(gca,'fontsize',14)
% xlabel('time (ms)');ylabel('Counts')
bar(1:twin:tsim,PSTHmit,'k')
xlim([-5 tsim]);ylim([0 nmit])
set(gca,'fontsize',14)
xlabel('time (ms)');ylabel('Counts')


ISIgraprox = [];
PSTHgraprox = zeros(1,dt*numtp/twin);
for ii = 1:ngraprox
    spiketimes = dt*find(GraProximal{ii}.S == 1);
    ISIgraprox = [ISIgraprox diff(spiketimes)];
    PSTHgraprox = PSTHgraprox + histc(spiketimes,1:twin:tsim);
end
 figure(212)
% subplot(1,2,1)
% hist(ISIgraprox,200)
% set(gca,'fontsize',14)
% xlabel('ISI (ms)');ylabel('Counts')
% title('Proximal Granule')
% subplot(1,2,2)
% bar(1:twin:tsim,PSTHgraprox,'k')
% xlim([0 tsim])
% set(gca,'fontsize',14)
% xlabel('time (ms)');ylabel('Counts')
bar(1:twin:tsim,PSTHgraprox,'k')
xlim([-5 tsim]);ylim([0 ngraprox])
set(gca,'fontsize',14)
xlabel('time (ms)');ylabel('Counts')


ISIpyr = [];
PSTHpyr = zeros(1,dt*numtp/twin);
for ii = 1:npyr
    spiketimes = dt*find(Pyramidal{ii}.S == 1);
    ISIpyr = [ISIpyr diff(spiketimes)];
    PSTHpyr = PSTHpyr + histc(spiketimes,1:twin:tsim);
end
figure(3)
subplot(1,2,1)
hist(ISIpyr,200)
set(gca,'fontsize',14)
xlabel('ISI (ms)');ylabel('Counts')
title('Pyramidal')
subplot(1,2,2)
bar(1:twin:tsim,PSTHpyr,'k')
xlim([0 tsim])
xlabel('time (ms)');ylabel('Counts')


% take a look at membrane voltage
figure(333)
choice = 2;
plot(timevec,GraProximal{choice}.V)
set(gca,'fontsize',14)
% title(['GraDistal',num2str(choice)])
% ylim([-0.08 0.01])
xlabel('time (ms)');ylabel('V')



% LFP plots with superimposed firing rates
figure(11)
hold on
mitH = plot(timevec,MitLFPs);
% bar(1:twin:tsim,0.005*PSTHmit,'k')
hold off
set(gca,'fontsize',14)
xlim([0 tsim])
% ylim([-0.15 0])
% title('Mitral')
xlabel('time (ms)');ylabel('sLFP')

figure(2)
hold on
graproxH = plot(timevec,GraProxLFPs(:,1));
% bar(1:twin:tsim,0.005*PSTHgraprox,'k')
hold off
set(gca,'fontsize',14)
xlim([0 tsim])
title('Proximal Granule')
xlabel('time (ms)');ylabel('LFP')

figure(3)
% gradistH = plot(timevec,GraDistLFPs(:,1));
gradistH = plot(timevec,GraDistLFPs.MitGradist(:,1));
set(gca,'fontsize',14)
xlim([0 tsim])
title('Distal Granule')
xlabel('time (ms)');ylabel('LFP')

% plot(timevec,(MitLFPs(:,1) + (0.1*GraDistLFPs(:,1)))/2)

figure(4)
hold on
pyrH = plot(timevec,PyrLFPs(:,1)+6);
% bar(1:twin:tsim,0.005*PSTHpyr,'k')
hold off
set(gca,'fontsize',14)
xlim([0 tsim])
title('Pyramidal')
xlabel('time (ms)');ylabel('LFP')

% save LFP plots
cd(save_path)
saveas(mitH,['LFP_mit_',app],'tiff')
saveas(graproxH,['LFP_graprox_',app],'tiff')
saveas(pyrH,['LFP_pyr_',app],'tiff')


%% Power Spectrum
trim = 5000; % trim beginning and end to avoid edge effects
% trim:end-50
L = length(timevec(trim:end-100));  % Length of simulation
NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
f = sampf/2*linspace(0,1,NFFT/2+1);

mitFFT = fft(detrend(MitLFPs(trim:end-100,1),'constant'),NFFT)/L;

% graproxFFT = fft(detrend(GraProxLFPs(trim:end-100,1),'constant'),NFFT)/L;
% gradistFFT = fft(detrend(GraDistLFPs(trim:end-100,1),'constant'),NFFT)/L;
% pyrFFT = fft(detrend(PyrLFPs(trim:end-100,1),'constant'),NFFT)/L;


figure(111)
mitH = plot(f,2*abs(mitFFT(1:NFFT/2+1)));
xlim([1 100])
ylim([0 0.025])
set(gca,'fontsize',13)
title('Mitral FFT')
xlabel('Frequency (Hz) ')
ylabel('Power')
figure(2)
graproxH = plot(f,2*abs(graproxFFT(1:NFFT/2+1)));
xlim([1 100])
set(gca,'fontsize',14)
title('Proximal Granule FFT')
xlabel('Frequency (Hz)')
ylabel('Power')
figure(3)
gradistH = plot(f,2*abs(gradistFFT(1:NFFT/2+1)));
xlim([1 100])
set(gca,'fontsize',14)
title('Distal Granule FFT')
xlabel('Frequency (Hz)')
ylabel('Power')
figure(4)
pyrH = plot(f,2*abs(pyrFFT(1:NFFT/2+1)));
xlim([1 100])
set(gca,'fontsize',14)
title('Pyramidal FFT')
xlabel('Frequency (Hz)')
ylabel('Power')

% save LFP power plots
cd(save_path)
saveas(mitH,['LFPpwr_mit_',app],'tiff') 
saveas(graproxH,['LFPpwr_graprox_',app],'tiff')
saveas(pyrH,['LFPpwr_pyr_',app],'tiff')

% 
% % mtm power
% Hs=spectrum.mtm(3);
% psd(Hs,GraProxLFPs(:,1),'Fs',sampf)
% xlim([0 0.2])

%% Input Currents with 2 compartment granule cell (prox and dist)
mitchoice = 1;
distchoice = 2;
proxchoice = 8;
pyrchoice = 3;

glomitH = figure(1)
plot(timevec,InputCurrent.Iglomit(mitchoice,:))
set(gca,'fontsize',14)
title(['Glo -> Mit',num2str(mitchoice),' input'])
xlabel('time (ms)');ylabel('current')
mitdistH = figure(2)
plot(timevec,InputCurrent.ImitgradistAMPA(distchoice,:))
set(gca,'fontsize',14)
title(['Mit -> Gradist AMPA',num2str(distchoice),' input'])
xlabel('time (ms)');ylabel('current')
mitdistH = figure(3)
plot(timevec,InputCurrent.ImitgradistNMDA(distchoice,:))
set(gca,'fontsize',14)
title(['Mit -> Gradist NMDA',num2str(distchoice),' input'])
xlabel('time (ms)');ylabel('current')
ylim([0 6e-3])
distmitH = figure(4)
plot(timevec,InputCurrent.Igradistmit(mitchoice,:))
set(gca,'fontsize',14)
title(['Gradist -> Mit',num2str(mitchoice),' input'])
xlabel('time (ms)');ylabel('current')
distproxH = figure(5)
plot(timevec,InputCurrent.Idistprox(proxchoice,:))
set(gca,'fontsize',14)
title(['Gradist -> Gra soma',num2str(proxchoice),' input'])
xlabel('time (ms)');ylabel('current')
ylim([0 8e-3])
proxdistH = figure(6)
plot(timevec,InputCurrent.Iproxdist(distchoice,:))
set(gca,'fontsize',14)
title(['Gra soma -> Gradist',num2str(distchoice),' input'])
xlabel('time (ms)');ylabel('current')
mitpyrH = figure(7)
plot(timevec,InputCurrent.Imitpyr(pyrchoice,:))
set(gca,'fontsize',14)
title(['Mit -> Pyr',num2str(pyrchoice),' input'])
xlabel('time (ms)');ylabel('current')
pyrproxH = figure(8)
plot(timevec,InputCurrent.Ipyrgraprox(proxchoice,:))
set(gca,'fontsize',14)
title(['Pyr -> Gra soma',num2str(proxchoice),' input'])
xlabel('time (ms)');ylabel('current')
pyrfbaH = figure(9)
plot(timevec,InputCurrent.Ipyrfba(pyrchoice,:))
set(gca,'fontsize',14)
title(['Pyr -> Fba',num2str(pyrchoice),' input'])
xlabel('time (ms)');ylabel('current')

% ImitgradistNMDA
imagesc(InputCurrent.ImitgradistNMDA);colorbar;set(gca,'fontsize',14)

% Pfiregradist
plot(0.0001:0.0001:1,InputCurrent.Pfiregradist(distchoice,:));ylabel('Pfire');xlabel('time (s)')

% Pfiregraprox
plot(0.0001:0.0001:1,InputCurrent.Pfiregraprox(proxchoice,:));ylabel('Pfire');xlabel('time (s)')


% save ICs
cd(save_path)
saveas(mitdistH,['Input_mit-dist_',app],'tiff') 
saveas(distmitH,['Input_dist-mit_',app],'tiff')
saveas(distproxH,['Input_dist-prox_',app],'tiff')
saveas(proxdistH,['Input_dist-prox_',app],'tiff')
saveas(mitpyrH,['Input_mit-pyr_',app],'tiff') 
saveas(pyrproxH,['Input_pyr-prox_',app],'tiff')


%% Calculate Spectrogram and Coherogram

t_scale = 0.3; % s
movingwin = [0.3 t_scale]; % [window length, step size] in (s)

dmat1 = PyrLFPs;
label1 = 'Pyramidal';
dmat2 = FfoLFPs;
label2 = 'Feedforward';

dmat1 = MitLFPs;
label1 = 'Mitral';
dmat2 = GraProxLFPs;
label2 = 'Granule';

dmat1 = MitLFPs;
label1 = 'Mitral';
dmat2 = PyrLFPs;
label2 = 'Pyramidal';


[Coh phi Cross Pwr_spec1 Pwr_spec2 t f] = get_spectra_generic(dmat1,dmat2,sampf,movingwin);
% Pwr_spec has dimensions #stim X #fq X #trials

nums = length(t); % number of time bins from cohgramc
numfq = length(f); % number of fqs from cohgramc

% visualize spectra
Pwr_spec1ave = zeros(nums,numfq);
Pwr_spec2ave = zeros(nums,numfq);
for ii = 1:(numtrials)
    Pwr_spec1ave = Pwr_spec1ave + squeeze(Pwr_spec1(:,:,ii));
    Pwr_spec2ave = Pwr_spec2ave + squeeze(Pwr_spec2(:,:,ii));
end
    Pwr_spec1ave = Pwr_spec1ave/numtrials;
    Pwr_spec2ave = Pwr_spec2ave/numtrials;

min_pwr = min([log(Pwr_spec1ave(:))+20;log(Pwr_spec2ave(:))+20]); % arbitrarily adding 20 to make (+) for ease of normalization
max_pwr = max([log(Pwr_spec1ave(:))+20;log(Pwr_spec2ave(:))+20]);
    
yheight = [2 f(end)];
figure(1)
subplot(2,1,1)
hold on
pcolor(t,f,(20+log(Pwr_spec1ave))'/max_pwr)
plot([0 0],yheight,'-','color',[0 0.5 0],'linewidth',3)
hold off
shading interp
lighting phong
set(gca,'fontsize',14)
% xlim([-1 4])
ylim(yheight)
title([ label1 ' Spectrum ']);
ylabel('fq (Hz)');
subplot(2,1,2)
hold on
pcolor(t,f,(20+log(Pwr_spec2ave))'/max_pwr)
plot([0 0],yheight,'-','color',[0 0.5 0],'linewidth',3)
hold off
shading interp
lighting phong
B=colorbar;
set(B, 'Position', [0.914 .1 0.034 .8150])
set(gca,'fontsize',14)
% xlim([-1 4])
ylim(yheight)
title([ label2 ' Spectrum ']);
ylabel('fq (Hz)');xlabel('time (s)');

COHave = zeros(nums,numfq);
for ii = 1:(numtrials)
    COHave = COHave + squeeze(Coh(:,:,ii));
end
    COHave = COHave/numtrials;

pcolor(t,f,COHave')
shading interp
lighting phong
set(gca,'fontsize',14)
ylim(yheight)
title([ label1 ' - ' label2 ' Coherence ']);
ylabel('fq (Hz)');xlabel('time (s)');
colorbar


% Pwr spec at fixed time
hold on
plot(f,log(Pwr_spec1ave(1,:)))
plot(f,log(Pwr_spec2ave(1,:)),'color','r')
hold off
set(gca,'fontsize',14)
xlabel('fq (Hz)');ylabel('log(PWR)');
legend('Mitral','Granule')
legend boxoff


% 
% Spikesyn_lesion = zeros(4,numtp);
% Gradesyn_lesion = zeros(4,numtp);

Spikesyn_lesion(4,:) = MitLFPs(:,1);
Gradesyn_lesion(4,:) = MitLFPs(:,1);

subplot(4,2,1)
plot(timevec,Spikesyn_lesion(1,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
title('Spiking Granule')
ylabel('LFP')
legend('E_{GABA} = -15mV')
subplot(4,2,3)
plot(timevec,Spikesyn_lesion(2,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
ylabel('LFP')
legend('E_{GABA} = -35mV')
subplot(4,2,5)
plot(timevec,Spikesyn_lesion(3,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
ylabel('LFP')
legend('E_{GABA} = -55mV')
subplot(4,2,7)
plot(timevec,Spikesyn_lesion(4,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
xlabel('time (ms)');ylabel('LFP')
legend('E_{GABA} = -75mV')

subplot(4,2,2)
plot(timevec,Gradesyn_lesion(1,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
title('Graded Granule')
ylabel('LFP')
legend('E_{GABA} = -15mV')
subplot(4,2,4)
plot(timevec,Gradesyn_lesion(2,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
ylabel('LFP')
legend('E_{GABA} = -35mV')
subplot(4,2,6)
plot(timevec,Gradesyn_lesion(3,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
ylabel('LFP')
legend('E_{GABA} = -55mV')
subplot(4,2,8)
plot(timevec,Gradesyn_lesion(4,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-0.8 0.1])
xlabel('time (ms)');ylabel('LFP')
legend('E_{GABA} = -75mV')







Spikesyn_intact = zeros(4,numtp);
Gradesyn_intact = zeros(4,numtp);

Spikesyn_intact(4,:) = MitLFPs(:,1);
Gradesyn_intact(4,:) = MitLFPs(:,1);

subplot(4,2,1)
plot(timevec,Spikesyn_intact(1,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
title('Spiking Granule')
ylabel('LFP')
legend('E_{GABA} = -15mV')
subplot(4,2,3)
plot(timevec,Spikesyn_intact(2,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
ylabel('LFP')
legend('E_{GABA} = -35mV')
subplot(4,2,5)
plot(timevec,Spikesyn_intact(3,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
ylabel('LFP')
legend('E_{GABA} = -55mV')
subplot(4,2,7)
plot(timevec,Spikesyn_intact(4,:),'b')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
xlabel('time (ms)');ylabel('LFP')
legend('E_{GABA} = -75mV')

subplot(4,2,2)
plot(timevec,Gradesyn_intact(1,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
title('Graded Granule')
ylabel('LFP')
legend('E_{GABA} = -15mV')
subplot(4,2,4)
plot(timevec,Gradesyn_intact(2,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
ylabel('LFP')
legend('E_{GABA} = -35mV')
subplot(4,2,6)
plot(timevec,Gradesyn_intact(3,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
ylabel('LFP')
legend('E_{GABA} = -55mV')
subplot(4,2,8)
plot(timevec,Gradesyn_intact(4,:),'r')
set(gca,'fontsize',14)
xlim([0 tsim])
ylim([-1.5 0.1])
xlabel('time (ms)');ylabel('LFP')
legend('E_{GABA} = -75mV')




y = (GraProximal{55}.V ) ./ (0.013);
J = y <= 0;
y(J) = 0;
J = y > 1;
y(J) = 1;
y = y.^3;
I = rand(length(y),1)' <= y * 1;





