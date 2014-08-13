
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