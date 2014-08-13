%% LFP fq as a function of GC Vthresh


clear all
noisy_InitNetwork_2CG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Vthresh = 1e-3*[-65 -60 -55 -50 -45];
% Vthresh = 1e-3*[-68 -45];
change_line = 164; % granule cell Vthresh is on line 164
pstring = 'FThresh ';

% pstring2 = 'wNMDAMI ';
pstring2 = 'Iext ';
% pstring2 = 'tauNMDA1 ';
% pstring2 = 'tauNMDA2 ';

numval = 0.0068:0.0002:0.0104;
% change_line2 = 166; % Gradistal wNMDA
change_line2 = 35; % Iext
% change_line2 = 157; % tauNMDA1
% change_line2 = 158; % tauNMDA2

FmaxMAT = zeros(length(Vthresh),length(numval));
maxpwrmit = zeros(length(Vthresh),length(numval));

for w = 1:length(numval)
    w
    TextCell = regexp( fileread('noisy_OB_PC_params_2CG.txt'), '\n', 'split');
    wstring = [pstring2,num2str(numval(w))];
    
    TextCell{change_line2} = sprintf('%s',wstring);
    fid = fopen('noisy_OB_PC_params_2CG.txt', 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);
    
    [FmaxMAT(:,w) maxpwrmit(:,w)] = Fmax_vs_Param_2CG(Vthresh,change_line,pstring,timevec,input_file, Delay, Wfrac, Hd);
end

% save('FmaxMAT.mat','FmaxMAT')


% 
% plot(Vthresh,Fmax,'-o')
% set(gca,'fontsize',14)
% xlabel('V_{thresh}');ylabel('f (Hz)')
% xlim([-70e-3 -54e-3])

% Load variables from Jul 22 folder

Vthresh = 1e-3*(-69:0.5:-55);
load Fmax1.mat
Fmax1 = Fmax;
load Fmax2.mat
Fmax2 = Fmax;
load FmaxBo1.mat
FmaxBo1 = Fmax;
load FmaxBo2.mat
FmaxBo2 = FmaxMAT;
load FmaxBo3.mat
FmaxBo3 = FmaxMAT;
load FmaxServ1.mat
FmaxServ1 = Fmax;
load FmaxServ2.mat
FmaxServ2 = Fmax;
load FmaxAMPA1.mat
FmaxAMPA1 = Fmax;
load FmaxNMDA1.mat
FmaxNMDA1 = Fmax;
load FmaxNMDA2.mat
FmaxNMDA2 = Fmax;

VthreshExt = 1e-3*(-70:1:0);
load FmaxNMDAext.mat
FmaxNMDAext = Fmax;
load FmaxAMPAext.mat
FmaxAMPAext = Fmax;
load FmaxBoExt.mat
FmaxBoExt = Fmax;


% plot(Vthresh,Fmax1,'-o',Vthresh,Fmax2,'-o')
% set(gca,'fontsize',14)
% xlabel('V_{thresh}');ylabel('f (Hz)')
% xlim([-70e-3 -55e-3])

% Combined model with errorbars
plot(Vthresh,FmaxBo1,'-o',Vthresh,FmaxBo2,'-o',Vthresh,FmaxBo3,'-o',Vthresh,FmaxServ1,'-o',Vthresh,FmaxServ2,'-o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f (Hz)')
xlim([-70e-3 -54e-3])

FmaxMAT = [FmaxBo1 FmaxBo2 FmaxBo3 FmaxServ1 FmaxServ2];
stdFmaxMAT = std(FmaxMAT');
mFmaxMAT = mean(FmaxMAT');

errorbar(Vthresh,mFmaxMAT,stdFmaxMAT,'o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('fq (Hz)')
xlim([-70e-3 -54e-3])



% looking at both NMDA runs
plot(Vthresh,FmaxBo1,'-o',Vthresh,FmaxAMPA1,'-o',Vthresh,FmaxNMDA1,'-o',Vthresh,FmaxNMDA2,'-o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f (Hz)')
legend('Composite','AMPA','NMDA')
legend boxoff
xlim([-70e-3 -54e-3])


% The extended Vthresh range (showing NMDA non-linearity)
plot(VthreshExt,FmaxBoExt,'-o',VthreshExt,FmaxAMPAext,'-o',VthreshExt,FmaxNMDAext,'-o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f (Hz)')
legend('Composite','AMPA','NMDA','location','best')
legend boxoff
xlim([-70e-3 -10e-3])


plot(VthreshExt,FmaxNMDAext,'-o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
xlim([-70e-3 0])



%% Fmax vs Vthresh vs wNMDA

% shorter range
Vthresh = 1e-3*(-65:0.5:-40);
wNMDA = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1];
load FmaxMAT
FmaxMAT_NMDA = FmaxMAT;
clear FmaxMAT

plot(Vthresh,FmaxMAT_NMDA,'-o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
legend(num2str(wNMDA'))
legend boxoff
xlim([-70e-3 0])

% full range (and more weights)
Vthresh_full = 1e-3*(-70:1:0);
wNMDA_full = logspace(-3,-1.2,15);
load FmaxMAT_Vt_wNMDA
FmaxMAT_Vt_wNMDA = FmaxMAT;
clear FmaxMAT

plot(Vthresh_full,FmaxMAT_Vt_wNMDA,'-o')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
legend(num2str(wNMDA_full(1:2:end)'),'location','eastoutside')
legend boxoff
xlim([-70e-3 0])

wNMDA_samp = wNMDA_full(1:end-3);

Vcrit_empirical = [-0.0665 -0.0652 -0.0645 -0.062 -0.0587 -0.056 -0.0524 -0.048 -0.04 -0.0335 -0.025 -0.014];% from eyeballing the values on the x axis
plot(wNMDA_samp,Vcrit_empirical,'.','markersize',13)
set(gca,'fontsize',14)
xlabel('W_{M-G}');ylabel('V_{CRIT} (mV)')


%% Fmax vs Vthresh vs Iext

Vthresh = 1e-3*(-70:1:0);
Iext = 1e-3 * [6.8 8 9.2 10.4];
load FmaxMAT_Vt_Iext
FmaxMAT_Vt_Iext = FmaxMAT;
clear FmaxMAT

plot(Vthresh,FmaxMAT_Vt_Iext,'o-')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
legend(num2str(Iext'),'location','eastoutside')
legend boxoff
xlim([-70e-3 0])


Vthresh = 1e-3*[-65 -60 -55 -50 -45];
Iext = 0.0068:0.0002:0.0104;
load FmaxMAT_Iext_Vt
FmaxMAT_Iext_Vt = FmaxMAT;
clear FmaxMAT
maxpwrmit_Iext_Vt = maxpwrmit;
clear maxpwrmit

figure(1)
plot(Iext,FmaxMAT_Iext_Vt','o-')
set(gca,'fontsize',14)
xlabel('I_{ext}');ylabel('f_{LFP} (Hz)')
legend(num2str(Vthresh'),'location','eastoutside')
legend boxoff
xlim([0.0067 0.0105])
figure(2)
plot(Iext,maxpwrmit_Iext_Vt','o-')
set(gca,'fontsize',14)
xlabel('I_{ext}');ylabel('Power')
legend(num2str(Vthresh'),'location','eastoutside')
legend boxoff
xlim([0.0067 0.0105])





%% Fmax vs Vthresh vs tauNMDA

% tauNMDA2
Vthresh = 1e-3*(-70:1:0);
tauNMDA2 = 5:15:200;
load FmaxMAT_Vt_tauNMDA2
FmaxMAT_Vt_tauNMDA2 = FmaxMAT;
clear FmaxMAT

% LFP fq vs Vthresh
figure(1)
hold on
for t = 1:length(tauNMDA2)
plot(Vthresh,FmaxMAT_Vt_tauNMDA2(:,t),'.-','color',[0,t/length(tauNMDA2),1-(t/length(tauNMDA2))])
end
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
legend(num2str(tauNMDA2'),'location','eastoutside')
legend boxoff
xlim([-70e-3 0]);ylim([0 95])

% LFP fq vs tauNMDA2
figure(2)
indVthresh = 6:4:50;
hold on
for ii = 1:length(indVthresh)
    plot(tauNMDA2,FmaxMAT_Vt_tauNMDA2(indVthresh(ii),:),'.-','color',[1-(ii/length(indVthresh)),0,ii/length(indVthresh)])
end
set(gca,'fontsize',14)
xlabel('\tau_{decay}');ylabel('f_{LFP} (Hz)')
legend(num2str(Vthresh(indVthresh)'),'location','eastoutside')
legend boxoff


%tauNMDA1
Vthresh = 1e-3*(-70:1:0);
tauNMDA1 = 1:1:10;
load FmaxMAT_Vt_tauNMDA1
FmaxMAT_Vt_tauNMDA1 = FmaxMAT;
clear FmaxMAT

% LFP fq vs Vthresh
figure(1)
hold on
for t = 1:length(tauNMDA1)
plot(Vthresh,FmaxMAT_Vt_tauNMDA1(:,t),'.-','color',[0,t/length(tauNMDA1),1-(t/length(tauNMDA1))])
end
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
legend(num2str(tauNMDA1'),'location','eastoutside')
legend boxoff
xlim([-70e-3 0])

% LFP fq vs tauNMDA1
figure(22)
indVthresh = 6:4:50;
hold on
for ii = 1:length(indVthresh)
    plot(tauNMDA1,FmaxMAT_Vt_tauNMDA1(indVthresh(ii),:),'.-','color',[1-(ii/length(indVthresh)),0,ii/length(indVthresh)])
end
hold off
set(gca,'fontsize',14)
xlabel('\tau_{rise}');ylabel('f_{LFP} (Hz)')
legend(num2str(Vthresh(indVthresh)'),'location','eastoutside')
legend boxoff




%% Fmax vs Vthresh Vs input uniformity fraction


clear all
noisy_InitNetwork_2CG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Vthresh = 1e-3*[-70:1:-30];
% Vthresh = 1e-3*[-68 -45];
change_line = 161; % granule cell Vthresh is on line 161
pstring = 'FThresh ';

numval = [0 0.25 0.5 0.75 1];

FmaxMAT = zeros(length(Vthresh),length(numval));
maxpwrmit = zeros(length(Vthresh),length(numval));



numtp = length(timevec);
dt = timevec(2)-timevec(1);
sampf = 1/(dt*1e-3);
numtrials = 1;

trim = 500; % MAKE SURE TRIM IS SET CORRECTLY!!!!

L = length(timevec(trim:end-100));  % Length of data that will be analyzed
NFFT = 2^nextpow2(L); % Next power of 2 from L
f = sampf/2*linspace(0,1,NFFT/2+1);
ROI = ceil(8/(f(2)-f(1))):ceil(140/(f(2)-f(1)));
% 



for n = 1:length(numval)
for v = 1:length(Vthresh)
    
    TextCell = regexp( fileread('noisy_OB_PC_params_2CG.txt'), '\n', 'split');
    wstring = [pstring,num2str(Vthresh(v))];
    
    TextCell{change_line} = sprintf('%s',wstring);
    fid = fopen('noisy_OB_PC_params_2CG.txt', 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);

    % ALWAYS CHECK tsim AND tfinal !!!!!!!
    [OSN Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] ...
        = noisy_VLFP_2CG(numtp, numtrials, input_file, Delay, Wfrac, numval(n), Hd);
    
    mitFFT = fft(detrend(MitLFPs(trim:end-100,1),'constant'),NFFT)/L;
    
    absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
    maxpwrmit(v,n) = max(absmitFFT(ROI));
    
    maxind = find(absmitFFT == maxpwrmit(v,n));
    FmaxMAT(v,n) = f(maxind);

    disp([num2str(length(Vthresh) - v),' more iterations'])
    
end
end


% tauNMDA2
Vthresh = 1e-3*(-70:1:-30);
Unif = numval;
load FmaxMAT_Vt_Uniform
FmaxMAT_Vt_Uniform = stuff;
clear stuff
load maxpwrmit_Vt_Uniform
maxpwrmit_Vt_Uniform = pwrstuff;
clear pwrstuff

yo = [' % Uniform';' % Uniform';' % Uniform';' % Uniform';' % Uniform'];

% LFP fq vs Vthresh
figure(1)
plot(Vthresh,FmaxMAT_Vt_Uniform,'.-')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
% legend([num2str(100*Unif') yo],'location','southeast')
% legend boxoff
xlim([-70e-3 -40e-3]);


% LFP PWR vs Vthresh
figure(2)
plot(Vthresh,maxpwrmit_Vt_Uniform,'.-')
set(gca,'fontsize',14)
xlabel('V_{thresh}');ylabel('Power')
legend([num2str(100*Unif') yo],'location','eastoutside')
legend boxoff
xlim([-70e-3 -40e-3]);




