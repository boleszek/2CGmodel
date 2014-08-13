function [Fmax maxpwrmit] = Fmax_vs_Param_2CG(ParamVec,change_line,pstring,timevec,input_file, Delay, Wfrac, Hd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function can only be ran after noisy_InitNetwork_2CG.m is called
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


Fmax = zeros(length(ParamVec),1);
maxpwrmit = zeros(length(ParamVec),1);


for v = 1:length(ParamVec)
    
    TextCell = regexp( fileread('noisy_OB_PC_params_2CG.txt'), '\n', 'split');
    wstring = [pstring,num2str(ParamVec(v))];
    
    TextCell{change_line} = sprintf('%s',wstring);
    fid = fopen('noisy_OB_PC_params_2CG.txt', 'w');
    fprintf(fid, '%s\n', TextCell{:});
    fclose(fid);

    % ALWAYS CHECK tsim AND tfinal !!!!!!!
    [OSN Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] ...
        = noisy_VLFP_2CG(numtp, numtrials, input_file, Delay, Wfrac, Hd);
    
    mitFFT = fft(detrend(MitLFPs(trim:end-100,1),'constant'),NFFT)/L;
    
    absmitFFT = 2*abs(mitFFT(1:NFFT/2+1));
    maxpwrmit(v) = max(absmitFFT(ROI));
    
    maxind = find(absmitFFT == maxpwrmit(v));
    Fmax(v) = f(maxind);

    disp([num2str(length(ParamVec) - v),' more iterations'])
    
end

