

clear all
Delay = [10 10 0.1 0.1 0.1 0.1 0.1 0.1];
% Delay = [LOTdelay OPdelay MitGABAdelay GraAMPAdelay GraGABAdelay
% PyrGABAdelay PyrAMPAdelay IntAMPAdelay]
input_file = 'noisy_OB_PC_params_2CG.txt';
[crap, nums] = textread(input_file,'%s %s',10);
dt = str2num(nums{2});
tsim = str2num(nums{3});
numtp = length(1:round(tsim/dt));
nmit = str2num(nums{7});
ngra = str2num(nums{8});
npyr = str2num(nums{9});
nfb = str2num(nums{10});



% set weight vector
% Wfrac = [linspace(0,1,40000) ones(1,60000)];
Wfrac = 1./(1+exp(-(linspace(-500,500,tsim/dt))/10)); % sigmoidal weight change over time
% timevec = dt*(1:numtp);
% plot(timevec,Wfrac,'k','linewidth',3)
% ylabel('W_{P-G} ');xlabel('Time (ms) ')

% design filter
sampf = 1/(dt*1e-3);
Fp = 150; % pass fq
Fst = 850; % stop fq
d=fdesign.lowpass('Fp,Fst,Ap,Ast',pi*Fp/sampf,pi*Fst/sampf,1,60);
Hd = design(d,'equiripple');
% fvtool(Hd);
timevec = dt*(1:numtp);






%% Mitral max Pwr and max FQ as a function of wMitGra and wGraMit


% NOTE!!!! change_line1 must match wstring1 below
change_line1 = 162; % wMitGraAMPA is line 162 in noisy_OB_PC_params_2CG.txt
% change_line1 = 163; % wMitGraNMDA is line 163 in noisy_OB_PC_params_2CG.txt
change_line2 = 126; % wGraMit is line 126 in noisy_OB_PC_params_2CG.txt


weights = logspace(-3,0,30);

maxpwrmit = zeros(length(weights));
fmaxmit = zeros(length(weights));
maxpwrgra = zeros(length(weights));
fmaxgra = zeros(length(weights));
RvsWmit = zeros(length(weights));
RvsWgra = zeros(length(weights));

% FFT params
trim = 100; % trim beginning and end to avoid edge effects
L = length(timevec(trim:end-100));  % Length of simulation
NFFT = 2^nextpow2(L); % Next power of 2 from length of simulation
f = sampf/2*linspace(0,1,NFFT/2+1);

for w = 1:length(weights)
    for v = 1:length(weights)

    TextCell = regexp( fileread('noisy_OB_PC_params.txt'), '\n', 'split');
    wstring1 = ['wAMPAMI ',num2str(weights(w))]; % use for wMitGraAMPA
    % wstring1 = ['wNMDAMI ',num2str(weights(w))]; % use for wMitGraNMDA
    wstring2 = ['wGABAGR ',num2str(weights(v))]; % use for wGraMit

    
    TextCell{change_line1} = sprintf('%s',wstring1);
    TextCell{change_line2} = sprintf('%s',wstring2);
    fid = fopen('noisy_OB_PC_params.txt', 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);


    numtrials = 1;

    % ALWAYS CHECK tsim AND tfinal !!!!!!!
    [OSN Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] ...
        = noisy_VLFP(numtp, numtrials, input_file, Delay, Wfrac, Hd);


    % perform FFT, find max fq
    %mit
    mitFFT = fft(detrend(MitLFPs(trim:end-100,1),'constant'),NFFT)/L;
    absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
    maxpwrmit(w,v) = max(absmitFFT);
    maxind = find(absmitFFT == maxpwrmit(w,v));
    fmaxmit(w,v) = f(maxind);
    %gra
    graFFT = fft(detrend(GraProxLFPs(trim:end-100,1),'constant'),NFFT)/L;
    absgraFFT = 2*abs(graFFT(1:NFFT/2+1));
    maxpwrgra(w,v) = max(absgraFFT);
    maxind = find(absgraFFT == maxpwrgra(w,v));
    fmaxgra(w,v) = f(maxind);
    
    
    % calculate firing rate (s)
    RvsWmit(w,v) = sum(Mitral{25}.S)/(tsim/1000);
    RvsWgra(w,v) = sum(GraProximal{25}.S)/(tsim/1000);
    
    end
end



figure(1)
surf(weights,weights,log(maxpwrmit),'EdgeColor','none');
set(gca,'fontsize',14,'xscale','log','yscale','log')
xlabel('wGraMit')
ylabel('wMitGra')
title('Mitral max LFP pwr')
colorbar
% caxis([log(min(maxpwrmit(:))) log(max(maxpwrmit(:)))])

figure(2)
surf(weights,weights,fmaxmit,'EdgeColor','none');
% surf(weights,weights,fmaxmit,'FaceColor', 'interp', 'edgecolor', 'none', 'FaceLighting', 'phong');
% camlight right;
set(gca,'fontsize',14,'xscale','log','yscale','log')
xlabel('wGraMit')
ylabel('wMitGra')
title('Mitral max LFP fq')
colorbar
% caxis([6 22])

figure(3)
surf(weights,weights,log(maxpwrgra),'EdgeColor','none');
set(gca,'fontsize',14,'xscale','log','yscale','log')
xlabel('wGraMit')
ylabel('wMitGra')
title('Granule max LFP pwr')
colorbar
% caxis([log(min(maxpwrmit(:))) log(max(maxpwrmit(:)))])

figure(4)
surf(weights,weights,fmaxgra,'EdgeColor','none');
% surf(weights,weights,fmaxgra,'FaceColor', 'interp', 'edgecolor', 'none', 'FaceLighting', 'phong');
% camlight left;
set(gca,'fontsize',14,'xscale','log','yscale','log')
xlabel('wGraMit')
ylabel('wMitGra')
title('Granule max LFP fq')
colorbar
% caxis([6 22])

figure(5)
surf(weights,weights,RvsWmit,'EdgeColor','none');
set(gca,'fontsize',14,'xscale','log','yscale','log')
xlabel('wGraMit')
ylabel('wMitGra')
title('Mitral firing rate')
colorbar

figure(6)
surf(weights,weights,RvsWgra,'EdgeColor','none');
set(gca,'fontsize',14,'xscale','log','yscale','log')
xlabel('wGraMit')
ylabel('wMitGra')
title('Granule firing rate')
colorbar






