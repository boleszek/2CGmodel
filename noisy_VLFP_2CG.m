function [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent MitLFPs GraProxLFPs GraDistLFPs FfoLFPs PyrLFPs] = noisy_VLFP_2CG(numtp, numtrials, input_file,Delay,Wfrac,Hd)


% Simulate LFP activity using membrane voltage

% INPUTS
%
% numtp       -      # time points in simulation
% numtrials   -      # times to run simulation

% OUTPUTS



MitLFPs = zeros(numtp,numtrials);
GraProxLFPs = zeros(numtp,numtrials);
GraDistLFPs = zeros(numtp,numtrials);
FfoLFPs = zeros(numtp,numtrials);
PyrLFPs = zeros(numtp,numtrials);


for ii = 1:numtrials
    [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback OSNsource param InputCurrent] = noisy_OB_PC_network_2CG(input_file,Delay,Wfrac);
    % Mitral LFP
    nmit  = length(Mitral);
    Mitral_filtered = zeros(nmit,length(Mitral{1}.V));
    for n = 1:nmit
        Mitral_filtered(n,:) = filtfilt(Hd.Numerator,1,Mitral{n}.V);
    end
    MitLFPs(:,ii) = sum(Mitral_filtered,1);
    MitLFPs(:,ii) = smoothts(MitLFPs(:,ii),'b');
    
    % Proximal Granule LFP
    ngra  = length(GraProximal);
    GraProximal_filtered = zeros(ngra,length(GraProximal{1}.V));
    for n = 1:ngra
        GraProximal_filtered(n,:) = filtfilt(Hd.Numerator,1,GraProximal{n}.V);
        GraProximal_filtered(n,:) = smoothts(GraProximal_filtered(n,:),'b');
    end
    GraProxLFPs(:,ii) = sum(GraProximal_filtered,1);
    
    % Distal Granule LFP
    ngra  = length(GraProximal);
    GraDistal_filtered = zeros(ngra,length(GraDistal{1}.V));
    for n = 1:ngra
        GraDistal_filtered(n,:) = filtfilt(Hd.Numerator,1,GraDistal{n}.V);
        GraDistal_filtered(n,:) = smoothts(GraDistal_filtered(n,:),'b');
    end
    GraDistLFPs(:,ii) = sum(GraDistal_filtered,1);
    
    % Feedforward LFP
    nffo  = length(Feedforward);
    Feedforward_filtered = zeros(ngra,length(Feedforward{1}.V));
    for n = 1:nffo
        Feedforward_filtered(n,:) = filtfilt(Hd.Numerator,1,Feedforward{n}.V);
        Feedforward_filtered(n,:) = smoothts(Feedforward_filtered(n,:),'b');
    end
    FfoLFPs(:,ii) = sum(Feedforward_filtered,1)';
    
    % Pyramidal LFP
    npyr  = length(Pyramidal);
    Pyramidal_filtered = zeros(ngra,length(Pyramidal{1}.V));
    for n = 1:npyr
        Pyramidal_filtered(n,:) = filtfilt(Hd.Numerator,1,Pyramidal{n}.V);
        Pyramidal_filtered(n,:) = smoothts(Pyramidal_filtered(n,:),'b');
    end
    PyrLFPs(:,ii) = sum(Pyramidal_filtered,1);
end

end
    