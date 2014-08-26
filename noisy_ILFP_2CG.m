function [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitILFPs GraProxILFPs GraDistILFPs FfoILFPs PyrILFPs] = noisy_ILFP_2CG(numtp, numtrials, input_file,Delay,Wfrac,Hd)


% Simulate LFP activity using membrane voltage

% INPUTS
%
% numtp       -      # time points in simulation
% numtrials   -      # times to run simulation

% OUTPUTS



MitILFPs.GradistMit = zeros(numtp,numtrials);
MitILFPs.GloMit = zeros(numtp,numtrials);
GraDistILFPs.MitGradist = zeros(numtp,numtrials);
GraProxILFPs.PyrGra = zeros(numtp,numtrials);
GraProxILFPs.GraGra = zeros(numtp,numtrials);
PyrILFPs.MitPyr = zeros(numtp,numtrials);
PyrILFPs.PyrPyr = zeros(numtp,numtrials);
PyrILFPs.FbaPyr = zeros(numtp,numtrials);
PyrILFPs.FfoPyr = zeros(numtp,numtrials);
FfoILFPs = zeros(numtp,numtrials);

%%% Must sum ALL currents into each neuron type
for ii = 1:numtrials
    [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback OSNsource param InputCurrent] = noisy_OB_PC_network_2CG_modified(input_file,Delay,Wfrac);
    % Mitral LFP
    nmit  = length(Mitral);
    % Filter input currrents into Mit (ignore respiration)
    Gradistmit_filtered = zeros(nmit,length(Mitral{1}.V));
%     Glomit_filtered = zeros(nmit,length(Mitral{1}.V));
    for n = 1:nmit
        Gradistmit_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Igradistmit(n,:));
        Gradistmit_filtered(n,:) = smoothts(Gradistmit_filtered(n,:),'b');
%         Glomit_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Iglomit(n,:));
%         Glomit_filtered(n,:) = smoothts(Glomit_filtered(n,:),'b');
    end
    MitILFPs.GradistMit(:,ii) = sum(Gradistmit_filtered,1);
%     MitILFPs.GloMit(:,ii) = sum(Glomit_filtered,1);
    
    
%     % Proximal Granule LFP
%     ngra  = length(GraProximal);
%     % Filter input currrents into proximal Gra
%     Mitgraprox_filtered = zeros(ngra,length(GraProximal{1}.V));
%     Pyrgraprox_filtered = zeros(ngra,length(GraProximal{1}.V));
%     Gragraprox_filtered = zeros(ngra,length(GraProximal{1}.V));
%     for n = 1:ngra
%         Mitgraprox_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Imitgraprox(n,:));
%         Mitgraprox_filtered(n,:) = smoothts(Mitgraprox_filtered(n,:),'b');
%         Pyrgraprox_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Ipyrgraprox(n,:));
%         Pyrgraprox_filtered(n,:) = smoothts(Pyrgraprox_filtered(n,:),'b');
%         Gragraprox_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Igragra(n,:));
%         Gragraprox_filtered(n,:) = smoothts(Gragraprox_filtered(n,:),'b');
%     end
%     GraILFPs.MitGra(:,ii) = sum(Mitgraprox_filtered,1);
%     GraILFPs.PyrGra(:,ii) = sum(Pyrgraprox_filtered,1);
%     GraILFPs.GraGra(:,ii) = sum(Gragraprox_filtered,1);
    
    % Distal Granule LFP
    ngra  = length(GraDistal);
    MitGradist_filtered = zeros(ngra,length(GraDistal{1}.V));
    for n = 1:ngra
        MitGradist_filtered(n,:) = filtfilt(Hd.Numerator,1,GraDistal{n}.V);
        MitGradist_filtered(n,:) = smoothts(MitGradist_filtered(n,:),'b');
    end
    GraDistILFPs.MitGradist(:,ii) = sum(MitGradist_filtered,1);
    

    
    % Pyramidal LFP
%     npyr  = length(Pyramidal);
%     MitPyr_filtered = zeros(npyr,length(Pyramidal{1}.V));
%     PyrPyr_filtered = zeros(npyr,length(Pyramidal{1}.V));
%     FbaPyr_filtered = zeros(npyr,length(Pyramidal{1}.V));
%     FfoPyr_filtered = zeros(npyr,length(Pyramidal{1}.V));
%     for n = 1:npyr
%         MitPyr_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Imitpyr(n,:));
%         MitPyr_filtered(n,:) = smoothts(MitPyr_filtered(n,:),'b');
%         PyrPyr_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Ipyrpyr(n,:));
%         PyrPyr_filtered(n,:) = smoothts(PyrPyr_filtered(n,:),'b');
%         FbaPyr_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Ifbapyr(n,:));
%         FbaPyr_filtered(n,:) = smoothts(FbaPyr_filtered(n,:),'b');
%         FfoPyr_filtered(n,:) = filtfilt(Hd.Numerator,1,InputCurrent.Iffopyr(n,:));
%         FfoPyr_filtered(n,:) = smoothts(FfoPyr_filtered(n,:),'b');
%     end
%     PyrILFPs.MitPyr(:,ii) = sum(MitPyr_filtered,1);
%     PyrILFPs.PyrPyr(:,ii) = sum(PyrPyr_filtered,1);
%     PyrILFPs.FbaPyr(:,ii) = sum(FbaPyr_filtered,1);
%     PyrILFPs.FfoPyr(:,ii) = sum(FfoPyr_filtered,1);
    
    
        % Feedforward LFP
%     nffo  = length(Feedforward);
%     Feedforward_filtered = zeros(nffo,length(Feedforward{1}.V));
%     for n = 1:nffo
%         Feedforward_filtered(n,:) = filtfilt(Hd.Numerator,1,Feedforward{n}.V);
%         Feedforward_filtered(n,:) = smoothts(Feedforward_filtered(n,:),'b');
%     end
%     FfoLFPs(:,ii) = sum(Feedforward_filtered,1)';
end

end
    