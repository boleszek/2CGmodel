%% LFP fq as a function of GC numval1


clear all
noisy_InitNetwork_2CG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Vthresh = 1e-3*[-70:1:-40];
% Vthresh = 1e-3*[-68 -45];

uniformity = 0:0.05:1;
% uniformity = [0.5 0.75];

change_line1 = 98; % Fthresh
pstring1 = 'FThresh ';
numval1 = Vthresh;


change_line2 = 10; % randfrac
pstring2 = 'randfrac ';
numval2 = uniformity;

FmaxMAT = zeros(length(numval1),length(numval2));
maxpwrMAT = zeros(length(numval1),length(numval2));


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


for n1 = 1:length(numval1)
    TextCell = regexp( fileread('noisy_OB_PC_params_2CG.txt'), '\n', 'split');
    wstring = [pstring1,num2str(numval1(n1))];
    
    TextCell{change_line1} = sprintf('%s',wstring);
    fid = fopen('noisy_OB_PC_params_2CG.txt', 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);
    for n2 = 1:length(numval2)
        
        TextCell = regexp( fileread('noisy_OB_PC_params_2CG.txt'), '\n', 'split');
        wstring = [pstring2,num2str(numval2(n2))];
        
        TextCell{change_line2} = sprintf('%s',wstring);
        fid = fopen('noisy_OB_PC_params_2CG.txt', 'w');
        fprintf(fid, '%s\n', TextCell{:});
        fclose(fid);
    
        % ALWAYS CHECK tsim AND tfinal !!!!!!!
        [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] ...
            = noisy_VLFP_2CG(numtp, numtrials, input_file, Delay, Wfrac, Hd);
        
        mitFFT = fft(detrend(MitLFPs(trim:end-100,1),'constant'),NFFT)/L;
        
        absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
        maxpwrMAT(n1,n2) = max(absmitFFT(ROI));
        
        maxind = find(absmitFFT == maxpwrMAT(n1,n2));
        FmaxMAT(n1,n2) = f(maxind);
    
        disp([num2str(length(numval1)*length(numval2) - ((n1-1)*length(numval2) + n2)),' more iterations'])
        
    end
end


save('FmaxMAT.mat','FmaxMAT')
save('maxpwrMAT.mat','maxpwrMAT')


% % tauNMDA2
% numval1 = 1e-3*(-70:1:-30);
% Unif = numval2;
% load FmaxMAT_Vt_Uniform
% FmaxMAT_Vt_Uniform = stuff;
% clear stuff
% load maxpwrMAT_Vt_Uniform
% maxpwrMAT_Vt_Uniform = pwrstuff;
% clear pwrstuff
% 
% yo = [' % Uniform';' % Uniform';' % Uniform';' % Uniform';' % Uniform'];
% 
% % LFP fq vs numval1
% figure(1)
% plot(numval1,FmaxMAT_Vt_Uniform,'.-')
% set(gca,'fontsize',14)
% xlabel('V_{thresh}');ylabel('f_{LFP} (Hz)')
% % legend([num2str(100*Unif') yo],'location','southeast')
% % legend boxoff
% xlim([-70e-3 -40e-3]);
% 
% 
% % LFP PWR vs numval1
% figure(2)
% plot(numval1,maxpwrMAT_Vt_Uniform,'.-')
% set(gca,'fontsize',14)
% xlabel('V_{thresh}');ylabel('Power')
% legend([num2str(100*Unif') yo],'location','eastoutside')
% legend boxoff
% xlim([-70e-3 -40e-3]);
% 



