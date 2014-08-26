function [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback OSNsource param InputCurrent] = noisy_OB_PC_network_2CG(inputFile,Delay,Wfrac)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is an adaption of two programs, bulbmain.m and piriformmain.m,
% originally developed by Licurgo de Almeida.
%
% description of the original model can be found in....
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
% inputfile - txt file with parameters
% Delay     - [LOTdelay OPdelay MitGABAdelay GraAMPAdelay GraGABAdelay PyrGABAdelay PyrAMPAdelay IntAMPAdelay]
%                in ms
% Wfrac     - a vector of length ntimepoints specifying how weights
%               increase over time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUTS:
% Mitral, GraProximal, GraDistal, Feedforward, Pyramidal,
% Feedback  - data structures including parameters and simulated neuronal activity
% OSNsource  -  odor matrix (simulates odor by stimulating subset of OSN)
% param      -  model parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Starting program

tic;

if strcmpi(inputFile(end - 2:end),'txt') % if we use a txt as input, the
    % program reads the txt file and create the variables, if we use a mat file,
    % the program loads the variables.
    
    % Set new rand seed.
%     s = RandStream.create('mt19937ar','seed',sum(100*clock));
%     RandStream.setDefaultStream(s);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Reading parameters from input file and creating neurons
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param.inputFile = 'noisy_OB_PC_params_2CG.txt';
    
    
    % Open the input file for reading
    fid1 = fopen(param.inputFile,'r');
    if fid1 == -1
        msgbox('Could not open the input file! Make sure the filname and path are correct.','ERROR');
        return;
    end
    
    str = fgetl(fid1);
    
    % Get network parameters.
    while ~strcmpi(str,'neurons')
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % network parameters
            otherwise
                param = SetNetworkParameters(param,str,Delay,Wfrac);
        end
        
        str = fgetl(fid1);
    end
    
    
    % Create cells

        OSNfile = strcat(param.outputPath,param.OSNsource);
        [param,OSNsource,Mitral,GraProximal,GraDistal,Feedforward,Pyramidal,Feedback] = CreateCells(param,OSNfile);

    
    str = fgetl(fid1);
    
    % Get cell parameters.
    while ~strcmpi(str,'end')
        
        switch lower(str)
            case '' % It's possible to use blank lines to organize the
                % neuronal parameters
            case 'mitral'
                celltype = 'mitral';
            case 'graproximal'
                celltype = 'graproximal';
            case 'gradistal'
                celltype = 'gradistal';
            case 'feedforward'
                celltype = 'feedforward';
            case 'pyramidal'
                celltype = 'pyramidal';
            case 'feedback'
                celltype = 'feedback';
            otherwise
                switch celltype
                    case 'mitral'
                        Mitral = SetNeuronParameters(Mitral,param.nMitral,str);
                    case 'graproximal'
                        GraProximal = SetNeuronParameters(GraProximal,param.nGraprox,str);
                    case 'gradistal'
                        GraDistal = SetNeuronParameters(GraDistal,param.nGradist,str);    
                    case 'feedforward'
                        Feedforward = SetNeuronParameters(Feedforward,param.nPyramidal,str);
                    case 'pyramidal'
                        Pyramidal = SetNeuronParameters(Pyramidal,param.nPyramidal,str);
                    case 'feedback'
                        Feedback = SetNeuronParameters(Feedback,param.nFeedback,str);
                end
                
        end
        
        str = fgetl(fid1);
    end

    
    fclose(fid1); % Close input file
    fname = inputFile(1:end - 3);
    fname = strcat(fname,'mat');
    save(fname,'Mitral','GraProximal','GraDistal','Feedforward','Pyramidal','Feedback','OSNsource','param');
    
elseif strcmpi(inputFile(end - 2:end),'mat') %if the input file is .mat
    % we have to
    load(inputFile,'Mitral','GraProximal','GraDistal','OSNsource','param');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Neuronal activity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent] = NeuroActivity(Mitral,GraProximal,GraDistal,Feedforward,Pyramidal,Feedback,param,OSNsource);


toc;
end

function param = SetNetworkParameters(param,str,Delay,Wfrac)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%-----------OB--------------OB----------------OB--------------------------
%
% This function sets the parameters for the different neurons
%
% Modified by Boleszek Osinski on 06/14/2013
%
% Licurgo de Almeida
% 12/20/2010
% Information not related with the parameters of different neuros.
% Path: path where we save input and output file
% dt: timestep (in ms)
% tsim: simulation time (in ms)
% tinit: stimulus begin (in ms)
% tfinal: stimulus end (in ms)
% nMitral: number of mitral cells
% nGradist: number of granule distal synapses
% nGraprox: number of granule cell soma
% GraGracon = if true, granule cells connect to each other
% DistalON = if true graded inhibitory distal Granule dendrites are present
% ProximalON = if true spiking inhibitory proximal Granule soma are present
% BulbH = Bulb height (in distance units)
% BulbW = Bulb width (in distance units)
% NoiseMit = Mitral cell noise std
% NoiseGraprox = Granule proximal dendrite noise std
% NoiseGradist = Granule distal dendrite noise std
% mFactor = multiplicative factor. number of granule cells = number of
% mitral cells * mFactor
% OSNsource = source of the OSN inputs.
% Odorant = odortant number.
% Grasource = source of the granule cells (used only for comparison
% (add the name of new parameters here)
% mitAHP = if true, mitral cells show AHP
% graAHP = if true, granule cells show AHP
% graLLD = if true, granule cells show LLD
% Respiration = if true, the input is modulated by an oscillation
% representing the respiration
% RespFreq = Respiratory frequency
% Iext = external current to OSNs
% Inoise = noise fraction of OSN input
% randfrac = fraction of uniform input to MCs
% SpikeV = Spike voltage
% CChanceGraMit = Chance of connection between Gra and Mit
% CChanceGraGra = Chance of connection between Gra
%
%
%-------------PC-------------PC-------------PC-------------PC-------------
%
%
% Modified by Boleszek Osinski on 06/14/2013
%
% Licurgo de Almeida
% 02/25/2011
% Information not related with the parameters of different neuros.
% * nPyramidal = number of pyramindal cells (and feedforward inhibitory neurons)
% * nFeedback = number of feedback inhibitory cells
% * CChancePyrPyr = Chance of connection between two different Pyramidal
% cells.
% * CChancePyrFba = Chance of connection between Pyramidal and Feedback
% cells.
% * CChanceFbaPyr = Chance of connection between Feedback and Pyramidal
% cells.
% * CChanceFfoPyr = Chance of connection between Feedforward cells and Pyramidal
% cells.
% * NoiseFfo = Feedforward neuron noise std
% * NoisePyr = Pyramidal cell noise std
% * NoiseFba = Feedback cell noise std
% * NoiseParam = true if we want a slighly variation on the parameters
% * Noiselevel = percetage of variation in each parameter
% * pyrAHP = if true, pyrAHP is ON
%
%
%-----------OBPC-----------OBPC-----------OBPC-----------OBPC-------------
%
% Boleszek Osinski 06/14/2013
%
% flagLOT    = if true Mi-Pyr connection is in tact
% LOTvar   = variability of LOT conduction delay (in ms)
% flagOP     = if true then Pyr-Gr connection is in tact
% OPvar    = variability of OP conduction delay (in ms)
% flagWLOTvar    = if true mit-pyr weights are scaled at each timestep with Wfrac
% flagWOPvar     = if true pyr-gra weights are scaled at each timestep with Wfrac
% flagWGRAMITvar = if true gra-mit weights are scaled at each timestep with Wfrac
% flagWMitGravar = if true mit-gra weights are scaled at each timestep with Wfrac
% CChancePyrGra  = Chance of connection between Pyramidal and Granule cells
% CChanceMitPyr  = Chance of connection between Mitral cells and Pyramidal cells
% CChanceMitFfo  = Chance of connection between Mitral cells and Feedforward cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.Wfrac = Wfrac; % Weight fraction for scaling weights over time
param.Delay = Delay; % Synaptic latencies

str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    
    % OB parameters
    case 'path'
        param.outputPath = ParValue; %path
    case 'dt'
        param.dt = str2num(ParValue); %step
    case 'tsim'
        param.tsim = str2num(ParValue); %simulation time
    case 'tinit'
        param.tinit = str2num(ParValue); %stimulus begin
    case 'tfinal'
        param.tfinal = str2num(ParValue); %stimulus end
    case 'nmitral'
        param.nMitral = str2num(ParValue); %number of Mitral cells
    case 'ngradist'
        param.nGradist = str2num(ParValue); %number of granule distal synapses
    case 'ngraprox'
        param.nGraprox = str2num(ParValue); %number of Granule cell soma
    case 'gragracon'
        param.GraGracon = str2num(ParValue); %Gra-Gra connections
    case 'distalon'
        param.DistalON = str2num(ParValue); %flag distal gra dendrites
    case 'proximalon'
        param.ProximalON = str2num(ParValue); %flag proximal gra dendrites
    case 'bulbh'
        param.BulbH = str2num(ParValue); %Bulb height
    case 'bulbw'
        param.BulbW = str2num(ParValue); %Bulb width
    case 'noisemit'
        param.noisemit = str2num(ParValue); %noise Mitral
    case 'noisegraprox'
        param.noisegraprox = str2num(ParValue); %noise granule prox
    case 'noisegradist'
        param.noisegradist = str2num(ParValue); %noise granule dist
    case 'preset'
        param.flagpreset = str2num(ParValue); %flag preset positions
    case 'mfactor'
        param.mFactor = str2num(ParValue); % multiplicative factor
    case 'osnsource'
        param.OSNsource = ParValue; %name of the file source for OSN information
    case 'odorant'
        param.Odorant = str2num(ParValue); %odorant number
    case 'grasource'
        param.Grasource = ParValue; %name of the file source for granule cells information
    case 'mitahp'
        param.mitAHP = str2num(ParValue); %flag presence of mitral AHP
    case 'graahp'
        param.graAHP = str2num(ParValue); %flag presence of granule AHP
    case 'gralld'
        param.graLLD = str2num(ParValue); %flag presence of granule LLD
    case 'respiration'
        param.flagRespiration = str2num(ParValue); %flag respiratory modulation
    case 'respfreq'
        param.RespFreq = str2num(ParValue); %respiratory frequency
    case 'iext'
        param.Iext = str2num(ParValue); %External current to OSNs
    case 'inoise'
        param.Inoise = str2num(ParValue); %noise fraction of OSN input
    case 'randfrac'
        param.randfrac = str2num(ParValue); %fraction of uniform input (scaling input to MCs)
    case 'spikev'
        param.SpikeV = str2num(ParValue); %Spike voltage
    case 'cchancegramit'
        param.CChanceGraMit = str2num(ParValue); %chance of connection between Granule and Mitral cells
    case 'cchancegragra'
        param.CChanceGraGra = str2num(ParValue); %chance of connection between Granule cells
            
        
        % PC parameters
    case 'npyramidal'
        param.nPyramidal = str2num(ParValue); %number of Pyramidal cells (and Feedforward)
    case 'nfeedback'
        param.nFeedback = str2num(ParValue); %number of Feedback cells
    case 'cchancepyrpyr'
        param.CChancePyrPyr = str2num(ParValue); %chance of connection between two Pyramidal cells
    case 'cchancepyrfba'
        param.CChancePyrFba = str2num(ParValue); %chance of connection between Pyramidal and Feedback cells
    case 'cchancefbapyr'
        param.CChanceFbaPyr = str2num(ParValue); %chance of connection between Feedback and Pyramidal cells
    case 'cchanceffopyr'
        param.CChanceFfoPyr = str2num(ParValue); %chance of connection between Feedforward and Pyramidal cells
    case 'ffoneurons'
        param.FfoNeurons = str2num(ParValue); %flag distal ffo neurons
    case 'noiseffo'
        param.noiseffo = str2num(ParValue); %noise Feedforward neurons
    case 'noisepyr'
        param.noisepyr = str2num(ParValue); %noise Pyramidal cells
    case 'noisefba'
        param.noisefba = str2num(ParValue); %noise Feedback neurons
    case 'noiseparam'
        param.NoiseParam = str2num(ParValue); %flag noise on the parameters
    case 'noiselevel'
        param.NoiseLevel = str2num(ParValue); %noise value
    case 'pyrahp'
        param.pyrAHP = str2num(ParValue); %flag use of pyrAHP

        
        % OB-PC parameters...
    case 'flaglot'
        param.flagLOT = str2num(ParValue); %if true Mi-Pyr connections are in tact
    case 'lotvar'
        param.LOTvar = str2num(ParValue); %variability in LOT conduction delay
    case 'flagop'
        param.flagOP = str2num(ParValue); %if true Pyr-Gr connections are in tact
    case 'opvar'
        param.OPvar = str2num(ParValue); %variability in OP conduction delay
    case 'flagwlotvar'
        param.flagWLOTvar = str2num(ParValue); %if true mit-gra weights scale with time by Wfrac     
    case 'flagwopvar'
        param.flagWOPvar = str2num(ParValue); %if true pyr-gra weights scale with time by Wfrac
    case 'flagwgramitvar'
        param.flagWGRAMITvar = str2num(ParValue); %if true gramit weights scale with time by Wfrac    
    case 'flagwmitgravar'
        param.flagWMitGravar = str2num(ParValue); %if true mitgra weights scale with time by Wfrac        
    case 'cchancepyrgra'
        param.CChancePyrGra = str2num(ParValue); %chance of connection between Pyramidal and Granule cells
    case 'cchancemitpyr'
        param.CChanceMitPyr = str2num(ParValue); %chance of connection between Mitral and Pyramidal cells
    case 'cchancemitffo'
        param.CChanceMitFfo = str2num(ParValue); %chance of connection between Mitral and Feedforward cells
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end


function [param,OSNsource,Mitral,GraProximal,GraDistal,Feedforward,Pyramidal,Feedback] = CreateCells(param,OSNfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates neurons based on the OSN file source (excel file)
%
% Modified by Boleszek Osinski on 06/16/2013 and 08/13/2014
% Licurgo de Almeida
% 12/21/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% input pattern to MCs
    aux_inputs = load(OSNfile); %read file with OSN information
%    OSNsource.inputs = ones(param.nMitral,1);
      OSNsource.inputs = aux_inputs.inputmatrand(param.Odorant,:);
      OSNsource.inputs(OSNsource.inputs<param.randfrac) = 1;
%    OSNsource.inputs = aux_inputs.inputmat(param.Odorant,:);
    %OSNsource.inputs = OSNsource.inputs / max(OSNsource.inputs); %%normalize data
    
    
    Mitral = cell(param.nMitral,1);
    for ii = 1:param.nMitral
        Mitral{ii}.input = []; %no input for now
        Mitral{ii}.label = 'Mitral';
    end
    
    GraProximal = cell(param.nGraprox,1);
    GraDistal = cell(param.nGradist,1);
    for ii = 1:param.nGraprox
        GraProximal{ii}.input = [];
        GraProximal{ii}.label = 'GraProximal';
    end
    for ii = 1:param.nGradist
        GraDistal{ii}.input = [];
        GraDistal{ii}.label = 'GraDistal';
    end


Pyramidal = cell(param.nPyramidal,1);
Feedforward = cell(param.nPyramidal,1);

for ii = 1:param.nPyramidal
    Pyramidal{ii}.input = []; %no input for now
    Pyramidal{ii}.label = 'Pyramidal';
    Feedforward{ii}.input = [];
    Feedforward{ii}.label = 'Feedforward'; %the number of pyramidal and
    % feedforward cells is always the same
end

Feedback = cell(param.nFeedback,1);

for ii = 1:param.nFeedback
    Feedback{ii}.input = []; %no input for now
    Feedback{ii}.label = 'Feedback';
end
end



function N = SetNeuronParameters(N,ncells,str)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the parameters for the different neurons
% Modified by Boleszek Osinski on 06/18/2013
%
%----------OB&PC-------------OB&PC----------OB&PC----------OB&PC----------
%
% OB and PC Neurons can present the following set of parameters:
% * tau = charging time constant of the neuron (ms). ms (not s) is the basic time unit in this program.
% * R = Membrane resistance (in ohms)
% * Fthresh = Average firing threshold (in V)
% * Vrest = resting potential
% * Vhyper = hyperpolarization potential
% * EAMPA = AMPA's Nernst potential.
% * EGABA = GABA's Nernst potential.
% * EAHP = AHP's Nernst potential.
% * ELLD = LLD's Nernst potential.
% * tauAMPA1 = AMPA's rising tau.
% * tauAMPA2 = AMPA's falling tau.
% * tauGABA1 = GABA's rising tau.
% * tauGABA2 = GABA's falling tau.
% * tauAHP = AHP time constant.
% * tauLLD = LLD time constant.
%
%
%-------------OB-------------OB-------------OB-------------OB-------------
%
% Licurgo de Almeida
% 11/02/2010
% OB Neurons can present the following set of parameters:
%
% * gmaxAMPA = AMPA's max conductance
% * gmaxGABA = GABA's max conductance
% * gmaxGABAP = Periglomerula GABA's max conductance (for Mitral cells
% only)
% * gmaxAHP = AHP's max conductance
% * IACh = Addition current when ACh is ON in non-spiking cells.
% * nGracon = number of granule cells connected to mitral or granule cells
% * CellRadius = radius of mitral and granule cells' dendrites
% * tauPROX1 = proximal Pyr-Gra synapse rising tau (Gra cell only)
% * tauPROX2 = proximal Pyr-Gra synapse falling tau (Gra cell only)
% * tauAMPA1 = distal Mit-Gra AMPA receptor rising tau (Gra cell only)
% * tauAMPA2 = distal Mit-Gra AMPA receptor falling tau (Gra cell only)
% * tauNMDA1 = distal Mit-Gra NMDA receptor rising tau (Gra cell only)
% * tauNMDA2 = distal Mit-Gra NMDA receptor falling tau (Gra cell only)
% * wAMPAMI = excitatory synaptic weight from Mitral cell to Granule cell AMPA
% * wNMDAMI = excitatory synaptic weight from Mitral cell to Granule cell NMDA
% * wGABAGR = inhibitory synaptic weight from Granule cell to Gra or Mit
% * wAMPAGL = excitatory synaptic weight from Glomerulus cell to Mitral cell
%-------------PC-------------PC-------------PC-------------PC-------------
%
% Licurgo de Almeida
% 02/28/2011
% PC Neurons cah present the following set of parameters:
%
% * wAMPA = excitatory synaptic weight
% * wAMPAPY = excitatory synaptic weight for Pyramidal association fibers
% (for Pyramidal cells only)
% * wGABA = inhibitory synaptic weight
% * wGABAFF = inhibitory synaptic weight from Feedforward cells (for
% Pyramidal cells only)
% * wAHP = AHP synaptic weight used to calculate AHP activation curves.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
str_aux = 1;
% Find parameter name
while str(str_aux) ~= ' '
    str_aux = str_aux + 1;
end

ParName = str(1:str_aux - 1); % name of the parameter
ParValue = str(str_aux + 1:end); % value of the parameter

switch lower(ParName)
    % OB and PC
    case 'tau'
        for ii = 1:ncells
            N{ii}.tau = str2num(ParValue); %time in ms
        end
    case 'r'
        for ii = 1:ncells
            N{ii}.R = str2num(ParValue); %resistance in ohms
        end
    case 'fthresh'
        for ii = 1:ncells
            N{ii}.FThresh = str2num(ParValue); %potential in volts
        end
    case 'vrest'
        for ii = 1:ncells
            N{ii}.Vrest = str2num(ParValue); %potential in volts
        end
    case 'vhyper'
        for ii = 1:ncells
            N{ii}.Vhyper = str2num(ParValue); %potential in volts
        end
    case 'noise'
        for ii = 1:ncells
            N{ii}.Noise = str2num(ParValue);
        end
    case 'eampa'
        for ii = 1:ncells
            N{ii}.EAMPA = str2num(ParValue); %potential in volts
        end
    case 'egaba'
        for ii = 1:ncells
            N{ii}.EGABA = str2num(ParValue); %potential in volts
        end
    case 'eahp'
        for ii = 1:ncells
            N{ii}.EAHP = str2num(ParValue); %potential in volts
        end
    case 'elld'
        for ii = 1:ncells
            N{ii}.ELLD = str2num(ParValue); %potential in volts
        end
    case 'tauampa1'
        for ii = 1:ncells
            N{ii}.tauAMPA1 = str2num(ParValue); %time in ms
        end
    case 'tauampa2'
        for ii = 1:ncells
            N{ii}.tauAMPA2 = str2num(ParValue); %time in ms
        end
    case 'taugaba1'
        for ii = 1:ncells
            N{ii}.tauGABA1 = str2num(ParValue); %time in ms
        end
    case 'taugaba2'
        for ii = 1:ncells
            N{ii}.tauGABA2 = str2num(ParValue); %time in ms
        end
    case 'tauahp'
        for ii = 1:ncells
            N{ii}.tauAHP = str2num(ParValue); %time in ms
        end
    case 'taulld'
        for ii = 1:ncells
            N{ii}.tauLLD = str2num(ParValue); %time in ms
        end
        
        % OB
    case 'gmaxampa'
        for ii = 1:ncells
            N{ii}.gmaxAMPA = str2num(ParValue); %AMPA channel
            % conductance in siemens
        end
    case 'gmaxgaba'
        for ii = 1:ncells
            N{ii}.gmaxGABA = str2num(ParValue); %GABA channel
            % conductance in siemens
        end
    case 'gmaxahp'
        for ii = 1:ncells
            N{ii}.gmaxAHP = str2num(ParValue); %"AHP" channel (for the simulation)
            % conductance in siemens
        end
    case 'iach'
        for ii = 1:ncells
            N{ii}.IACh = str2num(ParValue); %Additional current when ACh is ON
            % in non-spiking neurons
        end
    case 'gmaxgabap'
        for ii = 1:ncells
            N{ii}.gmaxGABAP = str2num(ParValue); %GABA channel from PG cells
            % conductance in siemens
        end
    case 'ngracon'
        for ii = 1:ncells
            N{ii}.NGraCon = str2num(ParValue); %number of granule cells
            % connected to a mitral or granule cell
        end
    case 'cellradius'
        for ii = 1:ncells
            N{ii}.CellRadius = str2num(ParValue); % (in distance units) If 
            % param.conntype is 'spatial' the user must set the radius of
            % mitral and granule cells
        end
    case 'tauprox1'
        for ii = 1:ncells
            N{ii}.tauPROX1 = str2num(ParValue); %time in ms
        end
    case 'tauprox2'
        for ii = 1:ncells
            N{ii}.tauPROX2 = str2num(ParValue); %time in ms
        end
    case 'taudist1'
        for ii = 1:ncells
            N{ii}.tauAMPA1 = str2num(ParValue); %time in ms
        end
    case 'taudist2'
        for ii = 1:ncells
            N{ii}.tauAMPA2 = str2num(ParValue); %time in ms
        end
    case 'taunmda1'
        for ii = 1:ncells
            N{ii}.tauNMDA1 = str2num(ParValue); %time in ms
        end
    case 'taunmda2'
        for ii = 1:ncells
            N{ii}.tauNMDA2 = str2num(ParValue); %time in ms
        end
    case 'wampami'
        for ii = 1:ncells
            N{ii}.wAMPAMI = str2num(ParValue); % excitatory synaptic weight from Mitral to Granule AMPA
        end
    case 'wnmdami'
        for ii = 1:ncells
            N{ii}.wNMDAMI = str2num(ParValue); % excitatory synaptic weight from Mitral to Granule NMDA
        end
    case 'wgabagr'
        for ii = 1:ncells
            N{ii}.wGABAGR = str2num(ParValue); % inhibitory synaptic weight from Granule to Gra or Mit
        end
    case 'wampagl'
        for ii = 1:ncells
            N{ii}.wAMPAGL = str2num(ParValue); % excitatory synaptic weight from Glom to Mit
        end
        
        % PC
    case 'wampa'
        for ii = 1:ncells
            N{ii}.wAMPA = str2num(ParValue); % excitatory synaptic weight
        end
    case 'wampapy'
        for ii = 1:ncells
            N{ii}.wAMPAPY = str2num(ParValue); % excitatory synaptic weight
            % from associative connections (pyramidal cells only)
        end
    case 'wgaba'
        for ii = 1:ncells
            N{ii}.wGABA = str2num(ParValue); % inhibitory synaptic weight
        end
    case 'wgabaff'
        for ii = 1:ncells
            N{ii}.wGABAFF = str2num(ParValue); % inhibitory synaptic weight
            % from feedforward cells (pyramidal cells only)
        end
    case 'wahp'
        for ii = 1:ncells
            N{ii}.wAHP = str2num(ParValue); % AHP amplitude
        end
        
        % New parameters must be added here...
        
    otherwise
        disp(['parameter ' ParName ' does not exist']);
        
end
end


function Cell = SetCellPosition(Cell,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the position on different neurons when no preset is
% presented.
% The main function of this program is neurogenesismain.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncells = length(Cell);

auxposy = param.BulbH / floor(sqrt(ncells));

if sqrt(ncells) == floor(sqrt(ncells))
    auxposx = param.BulbW / sqrt(ncells);
else
    auxposx = param.BulbW / (sqrt(ncells) + 1);
end

countcell = 0;
x = 0;
y = 0;

while countcell < ncells
    countcell = countcell + 1;
    if y >= param.BulbH
        y = 0;
        x = x + auxposx;
    end
    Cell{countcell}.X = mean([x,(x + auxposx)]);
    Cell{countcell}.Y = mean([y,(y + auxposy)]);
    y = y + auxposy;
end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mitral GraProximal GraDistal Feedforward Pyramidal Feedback param InputCurrent] = NeuroActivity(Mitral,GraProximal,GraDistal,Feedforward,Pyramidal,Feedback,param,OSNsource)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates the neuronal activity
%
% Licurgo de Almeida
% 11/03/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create initial setup -- OB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if param.flagRespiration == true
    Respiration = CreateRespFreq(param);
else
    Respiration = ones(1,round(param.tsim / param.dt));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%   OB   %%%%%

nds = param.nGradist/param.nGraprox; % number of distal synapses per granule cell


    % random connectivity
        MatGradistMit = SetConnections(param.nMitral,param.nGradist,param.CChanceGraMit);
        MatMitGradist = MatGradistMit'; % reciprocal synapses
        MatProxProx = SetConnections(param.nGraprox,param.nGraprox,param.CChanceGraGra);
        % There are nds distal synapses to each granule cell soma. Note: nGradist must = nds*nGraprox 
        % (this is hard coded)
        MatProxDist = zeros(param.nGradist,param.nGraprox);
        for ii = 1:param.nGraprox
            MatProxDist((nds*(ii-1)+1):nds*ii,ii) = 1;
        end

%%%%%   PC   %%%%%

    MatPyrFba = SetConnections(param.nFeedback,param.nPyramidal,param.CChancePyrFba);
    MatFbaPyr = SetConnections(param.nPyramidal,param.nFeedback,param.CChanceFbaPyr);
    MatFfoPyr = SetConnections(param.nPyramidal,param.nPyramidal,param.CChanceFfoPyr);
    MatPyrPyr = SetConnections(param.nPyramidal,param.nPyramidal,param.CChancePyrPyr);


%%%%%   OB-PC   %%%%%

    MatPyrGraprox = SetConnections(param.nGraprox,param.nPyramidal,param.CChancePyrGra);
    MatMitPyr = SetConnections(param.nPyramidal,param.nMitral,param.CChanceMitPyr);
    MatMitFfo = SetConnections(param.nPyramidal,param.nMitral,param.CChanceMitFfo);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set weights matrix -- PC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%   OB   %%%%%

    wGloMit = eye(param.nMitral)*Mitral{1}.wAMPAGL;
    wPglMit = wGloMit;



    wGABAMit = zeros(param.nMitral,1);
    for ii = 1:param.nMitral
        wGABAMit(ii) = Mitral{ii}.wGABAGR;
    end
    wGradistMit = SetWeights(MatGradistMit,wGABAMit,param);
    
    wAMPAGradist = zeros(param.nGradist,1);
    wNMDAGradist = zeros(param.nGradist,1);
    wGABAGraProx = zeros(param.nGraprox,1);
    for ii = 1:param.nGradist
        wAMPAGradist(ii) = GraDistal{ii}.wAMPAMI;
        wNMDAGradist(ii) = GraDistal{ii}.wNMDAMI;
    end
    for ii = 1:param.nGraprox
        wGABAGraProx(ii) = GraProximal{ii}.wGABAGR;
    end
    wMitGradistAMPA = SetWeights(MatMitGradist,wAMPAGradist,param);
    wMitGradistNMDA = SetWeights(MatMitGradist,wNMDAGradist,param);
    wProxProx = SetWeights(MatProxProx,wGABAGraProx,param);



%%%%%   PC   %%%%%

    wGABAPyr = zeros(param.nPyramidal,1);
    wAMPAPYPyr = wGABAPyr;
    wGABAFFPyr = wGABAPyr;
    for ii = 1:param.nPyramidal
        wGABAPyr(ii) = Pyramidal{ii}.wGABA;
        wAMPAPYPyr(ii) = Pyramidal{ii}.wAMPAPY;
        wGABAFFPyr(ii) = Pyramidal{ii}.wGABAFF;
    end
    wFbaPyr = SetWeights(MatFbaPyr,wGABAPyr,param);
    wPyrPyr = SetWeights(MatPyrPyr,wAMPAPYPyr,param);
    wFfoPyr = SetWeights(MatFfoPyr,wGABAFFPyr,param);

    wAMPAFba = zeros(param.nFeedback,1);
    for ii = 1:param.nFeedback
        wAMPAFba(ii) = Feedback{ii}.wAMPA;
    end
    wPyrFba = SetWeights(MatPyrFba,wAMPAFba,param);



%%%%%   OB-PC   %%%%%

    wAMPAPyr = zeros(param.nPyramidal,1);
    wAMPAFfo = wAMPAPyr;
    for ii = 1:param.nPyramidal
        wAMPAPyr(ii) = Pyramidal{ii}.wAMPA;
        wAMPAFfo(ii) = Feedforward{ii}.wAMPA;
    end
    wMitPyr = SetWeights(MatMitPyr,wAMPAPyr,param);
    wMitFfo = SetWeights(MatMitFfo,wAMPAFfo,param);
            
    wAMPAGraprox = zeros(param.nGraprox,1);
    for ii = 1:param.nGraprox
        wAMPAGraprox(ii) = GraProximal{ii}.wAMPAPY;
    end
    
    wPyrGraprox = SetWeights(MatPyrGraprox,wAMPAGraprox,param);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set OSN parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OSN external input
% Iextosn = SetExtInput(param,OSN);

Iextosn = zeros(param.nMitral,round(param.tsim / param.dt));

for ii = 1:param.nMitral
    Iextosn(ii,round(param.tinit / param.dt) : round(param.tfinal / param.dt)) = param.Iext * OSNsource.inputs(ii);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Stores the voltage of each Mitral cell at a given time
Vmit = zeros(param.nMitral,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Smit = Vmit;
VAHPmit = Vmit;

refracmit = 3; % Refractory period after firing (ms)


Restmit = zeros(param.nMitral,1); % resting potential
Threshmit = Restmit; % current firing threshold
Hypermit = Restmit; % hyperpolarization potential
Rmit = Restmit; % membrane resistance
taumit = Restmit; % tau neuron
countrefracmit = Restmit; % Refractory period counter. This is assumed to
% be the same for all Mitral cells so we don't need to add it to the for
% below.
gmaxAMPAmit = Restmit; % Max AMPA conductance
tauAMPA1mit = Restmit; % AMPA's rising tau.
tauAMPA2mit = Restmit; % AMPA's falling tau.
EAMPAmit = Restmit; % AMPA's Nernst potential.
gmaxGABAmit = Restmit; % Max GABA conductance
tauGABA1mit = Restmit; % GABA's rising tau.
tauGABA2mit = Restmit; % GABA's falling tau.
EGABAmit = Restmit; % GABA's Nernst potential.

if param.mitAHP == true
    gmaxAHPmit = Restmit; % Max AHP conductance
    tauAHPmit = Restmit; % AHP tau.
    EAHPmit = Restmit; % AHP's Nernst potential.
end


% ...New parameters should be added here and inside the 'for'...


for ii = 1:param.nMitral
    Restmit(ii) = Mitral{ii}.Vrest;
    Hypermit(ii) = Mitral{ii}.Vhyper;
    Rmit(ii) = Mitral{ii}.R;
    taumit(ii) = Mitral{ii}.tau;
    gmaxAMPAmit(ii) = Mitral{ii}.gmaxAMPA;
    tauAMPA1mit(ii) = Mitral{ii}.tauAMPA1;
    tauAMPA2mit(ii) = Mitral{ii}.tauAMPA2;
    EAMPAmit(ii) = Mitral{ii}.EAMPA;
    gmaxGABAmit(ii) = Mitral{ii}.gmaxGABA;
    tauGABA1mit(ii) = Mitral{ii}.tauGABA1;
    tauGABA2mit(ii) = Mitral{ii}.tauGABA2;
    EGABAmit(ii) = Mitral{ii}.EGABA;
    if param.mitAHP == true
        gmaxAHPmit(ii) = Mitral{ii}.gmaxAHP;
        tauAHPmit(ii) = Mitral{ii}.tauAHP;
        EAHPmit(ii) = Mitral{ii}.EAHP;
    end
    

    Threshmit(ii) = Mitral{ii}.FThresh;
    
    
    Mitral{ii}.Connections = MatGradistMit(ii,:);
    Mitral{ii}.ConWeights = wGradistMit(ii,:);
end

% Initialize Mitral cells potentials
Vmit(:,1) = Restmit;
Vmit_nospike = Vmit;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0mit = zeros(param.nMitral,1) - 10000000;

Iglomit = zeros(param.nMitral,1); % Input coming from Glo
Iglomit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));
maxgAMPAmit = getmaxg(param.dt,tauAMPA1mit,tauAMPA2mit); % Get max conductance
% amplitude

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Mitral cell can be
% connected to a varied number of granule cells
tGABA0mit = zeros(param.nGradist,round(param.tsim / param.dt)) - 10000000; %NOTE!!! for 2CG model tGABA0mit exists for each dendrite!
Igradistmit = zeros(param.nMitral,1); % Input coming from Granule cells
Igradistmit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));
% only used when distal gra dendrites are removed
Igraproxmit = zeros(param.nMitral,1); % Input coming from Granule cells
Igraproxmit_matrix = zeros(param.nMitral,round(param.tsim / param.dt));

maxgGABAmit = getmaxg(param.dt,tauGABA1mit,tauGABA2mit); % Get max conductance
% amplitude


% AHP time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0.
if param.mitAHP == true
    tAHP0mit = zeros(param.nMitral,1) - 10000000;
%     maxgAHPmit = getmaxg(param.dt,tauAHP1mit,tauAHP2mit); % Get max conductance
% not using the rising and falling taus for AHP current
end
Iahpmit = zeros(param.nMitral,1); % AHP current
% amplitude


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Mitral cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Proximal Granule cell parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.ProximalON == true
    
% Stores the voltage of each Granule cell at a given time
Vgraprox = zeros(param.nGraprox,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sgraprox = Vgraprox;
% Probability of each cell firing
Pfiregraprox = Vgraprox;

VAHPgraprox = Vgraprox;
VLLDgraprox = Vgraprox;

refracgraprox = 500; % Refractory period after firing (ms)

Restgraprox = zeros(param.nGraprox,1); % resting potential
Threshgraprox = Restgraprox; % current firing threshold
Hypergraprox = Restgraprox; % hyperpolarization potential
Rgraprox = Restgraprox; % membrane resistance
taugraprox = Restgraprox; % tau neuron
countrefracgraprox = Restgraprox; % Refractory period counter. This is assumed to
% be the same for all Granule cells so we don't need to add it to the for below.

gmaxAMPAgraprox = Restgraprox; % Max AMPA conductance

tauPROX1 = Restgraprox; % AMPA's rising tau.
tauPROX2 = Restgraprox; % AMPA's falling tau.

EAMPAgraprox = Restgraprox; % AMPA's Nernst potential.
gmaxGABAgraprox = Restgraprox; % Max GABA conductance
tauGABA1graprox = Restgraprox; % GABA's rising tau.
tauGABA2graprox = Restgraprox; % GABA's falling tau.
EGABAgraprox = Restgraprox; % GABA's Nernst potential.

if param.graAHP == true
    gmaxAHPgraprox = Restgraprox; % Max AHP conductance
    tauAHPgraprox = Restgraprox; % AHP tau.
    EAHPgraprox = Restgraprox; % AHP's Nernst potential.
end
if param.graLLD == true
    tauLLDgraprox = Restgraprox; % LLD tau.
    ELLDgraprox = Restgraprox; % LLD's Nernst potential.
end


MahpCongra = eye(param.nGraprox); % matrix with AHP connection
Wahpgra = MahpCongra; % AHP weights


% ...New parameters should be added here and inside the 'for'...

for ii = 1:param.nGraprox
    Restgraprox(ii) = GraProximal{ii}.Vrest;
    Hypergraprox(ii) = GraProximal{ii}.Vhyper;
    Rgraprox(ii) = GraProximal{ii}.R;
    taugraprox(ii) = GraProximal{ii}.tau;
    gmaxAMPAgraprox(ii) = GraProximal{ii}.gmaxAMPA;
    tauPROX1(ii) = GraProximal{ii}.tauPROX1;
    tauPROX2(ii) = GraProximal{ii}.tauPROX2;
    EAMPAgraprox(ii) = GraProximal{ii}.EAMPA;
    gmaxGABAgraprox(ii) = GraProximal{ii}.gmaxGABA;
    tauGABA1graprox(ii) = GraProximal{ii}.tauGABA1;
    tauGABA2graprox(ii) = GraProximal{ii}.tauGABA2;
    EGABAgraprox(ii) = GraProximal{ii}.EGABA;
    if param.graAHP == true
        gmaxAHPgraprox(ii) = GraProximal{ii}.gmaxAHP;
        tauAHPgraprox(ii) = GraProximal{ii}.tauAHP;
        EAHPgraprox(ii) = GraProximal{ii}.EAHP;
    end
    if param.graLLD == true
        tauLLDgraprox(ii) = GraProximal{ii}.tauLLD;
        ELLDgraprox(ii) = GraProximal{ii}.ELLD;
    end
    
    Threshgraprox(ii) = GraProximal{ii}.FThresh;
    
    GraProximal{ii}.PyrgraConnections = MatPyrGraprox(ii,:);
end

% Initialize Granule cells potentials
Vgraprox(:,1) = Restgraprox;
Vgraprox_nospike = Vgraprox;

% bulbar input current to graprox
Idistprox = zeros(param.nGraprox,1); % Input coming from dist to prox
Idistprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));
Imitgraprox = zeros(param.nGraprox,1);
Imitgraprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));

% input voltage (for direct distprox summing)
Vdistprox = zeros(param.nGraprox,1); % V input coming from dist to prox
Vdistprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));

% pyramidal input current to graprox
Ipyrgraprox = zeros(param.nGraprox,1); % Input coming from pyr to prox
Ipyrgraprox_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0

% tAMPA0pyr_gra is a matrix that stores tau values delayed by OP delay time
tAMPA0pyr_gra = zeros(param.nPyramidal,round(param.tsim / param.dt)) - 10000000;
% generate vector of OP delay variabilities from each pyr neuron within time
% resolution of dt (in ms)
OP_delay_vars = param.dt*ceil((2*param.OPvar*rand(param.nPyramidal,1)-param.OPvar)/param.dt);
% tAMPA0mit_gra is a matrix that stores tau values delayed by mit-gra delay time
tAMPA0mit_gra = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;


maxgAMPAgraprox = getmaxg(param.dt,tauPROX1,tauPROX2); % Get max conductance amplitude

% maxgAMPAgraprox = 1; % Set gmax = 1 because it will be normalized anyways

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Granule cell can be
% connected to a varied number of other granule cells
tGABA0graprox = zeros(param.nGraprox,round(param.tsim / param.dt)) - 10000000;
Igragra = zeros(param.nGraprox,1); % Input coming from other granule cells
Igragra_matrix = zeros(param.nGraprox,round(param.tsim / param.dt));

maxgGABAgraprox = getmaxg(param.dt,tauGABA1graprox,tauGABA2graprox); % Get max conductance amplitude

% maxgGABAgraprox = 1; % Set gmax = 1 because it will be normalized anyways

% AHP time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0.
if param.graAHP == true
    tAHP0graprox = zeros(param.nGraprox,1) - 10000000;
%     maxgAHPgra = getmaxg(param.dt,tauAHP1gra,tauAHP2gra); % Get max
%     conductance
end
if param.graLLD == true
    tLLD0graprox = zeros(param.nGraprox,1) - 10000000;
%     maxgLLDgra = getmaxg(param.dt,tauLLD1gra,tauLLD2gra); % Get max
%     conductance
end

% !!!!!!!!!!! Get rid of AHP current since it will be calculated in main
% for loop
% Iahpgraprox = zeros(param.nGranule,1); % AHP current
% % amplitude

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Proximal Granule cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Distal Granule cell parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param.DistalON == true
% Stores the voltage of each Granule cell at a given time
Vgradist = zeros(param.nGradist,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sgradist= Vgradist;


Pfiregradist = Vgradist; %Pfiregradist is the probability of firing when 
                             %we use graded inhibition from distal gra dendrites

refracgradist = 2; % Refractory period after firing (ms)

Restgradist = zeros(param.nGradist,1); % resting potential
Threshgradist = Restgradist; % current firing threshold
Hypergradist = Restgradist; % hyperpolarization potential
Rgradist = Restgradist; % membrane resistance
taugradist = Restgradist; % tau neuron
countrefracgradist = Restgradist; % Refractory period counter. This is assumed to
% be the same for all Granule cells so we don't need to add it to the for below.
gmaxAMPAgradist = Restgradist; % Max AMPA conductance

tauAMPA1 = Restgradist; % AMPA's rising tau.
tauAMPA2 = Restgradist; % AMPA's falling tau.
tauNMDA1 = Restgradist; % NMDA's rising tau.
tauNMDA2 = Restgradist; % NMDA's falling tau.

EAMPAgradist = Restgradist; % AMPA's Nernst potential.


% ...New parameters should be added here and inside the 'for'...

% Graded inhibition??????????????

for ii = 1:param.nGradist
    Restgradist(ii) = GraDistal{ii}.Vrest;
    Hypergradist(ii) = GraDistal{ii}.Vhyper;
    Rgradist(ii) = GraDistal{ii}.R;
    taugradist(ii) = GraDistal{ii}.tau;
    gmaxAMPAgradist(ii) = GraDistal{ii}.gmaxAMPA;
    tauAMPA1(ii) = GraDistal{ii}.tauAMPA1;
    tauAMPA2(ii) = GraDistal{ii}.tauAMPA2;
    tauNMDA1(ii) = GraDistal{ii}.tauNMDA1;
    tauNMDA2(ii) = GraDistal{ii}.tauNMDA2;
    EAMPAgradist(ii) = GraDistal{ii}.EAMPA;
    ENMDAgradist(ii) = GraDistal{ii}.EAMPA; % nmda reversal potential is the same as AMPA
    % including prox parameters here because they are used for pyr input 
    % onto distal dendrites when proximal compartment is removed.
    %NOTE: This code expects both prox and dist to be ON, so this is
    %unnecessary
%     tauPROX1(ii) = GraProximal{ii}.tauPROX1;
%     tauPROX2(ii) = GraProximal{ii}.tauPROX2;
%     EAMPAgraprox(ii) = GraProximal{ii}.EAMPA;
    

    Threshgradist(ii) = GraDistal{ii}.FThresh;
    
    GraDistal{ii}.Connections = MatMitGradist(ii,:);
end

% Initialize Granule cells potentials
Vgradist(:,1) = Restgradist;

% proximal granule or pyramidal input currents to gradist
Iproxdist = zeros(param.nGradist,1); % Input coming from prox to dist
Iproxdist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));
Ipyrgradist = zeros(param.nGradist,1); % Input coming from pyr to dist
Ipyrgradist_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

% mitral input to gradist AMPA and NMDA receptors
ImitgradistAMPA = zeros(param.nGradist,1); % Input coming from mitral cells
ImitgradistAMPA_matrix = zeros(param.nGradist,round(param.tsim / param.dt));

ImitgradistNMDA = zeros(param.nGradist,1); % Input coming from mitral cells
ImitgradistNMDA_matrix = zeros(param.nGradist,round(param.tsim / param.dt));


% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0gradist = zeros(param.nMitral,1) - 10000000;

% tAMPA0pyr_gra is a matrix that stores tau values delayed by OP delay time
tAMPA0pyr_gra = zeros(param.nPyramidal,round(param.tsim / param.dt)) - 10000000;
% generate vector of OP delay variabilities from each pyr neuron within time
% resolution of dt (in ms)
OP_delay_vars = param.dt*ceil((2*param.OPvar*rand(param.nPyramidal,1)-param.OPvar)/param.dt);
% tAMPA0mit_gra is a matrix that stores tau values delayed by mit-gra delay time
tAMPA0mit_gra = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;
% Wmitgradist = ones(param.nGranule,param.nMitral); % Synaptic weights

maxgAMPAgradist = getmaxg(param.dt,tauAMPA1,tauAMPA2); % Get max conductance amplitude for AMPA
maxgNMDAgradist = getmaxg(param.dt,tauNMDA1,tauNMDA2); % Get max conductance amplitude for NMDA

% Using proximal parameter for pyr input onto distal dendrites when proximal compartment is removed
maxgAMPAgraprox = getmaxg(param.dt,tauPROX1,tauPROX2); % Get max conductance amplitude

% maxgAMPAgradist = 1; % Set gmax = 1 because it will be normalized anyways

% Initialize vesicle fraction variable
% Vtau = 6;
% Vreplenish = 1./(exp(-(-5:0.1:5)/0.7)+1); % sigmoid function curve models the replenishing population of vesicles
% VFRACgradist = zeros(param.nGradist,round(param.tsim / param.dt) + length(Vreplenish));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Distal Granule cell synapse parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Feedforward neurons parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Feedforward cell at a given time
Vffo = zeros(param.nPyramidal,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sffo = Vffo;
Cmitffo = Vffo;

refracffo = 30; % Refractory period after firing (ms)
countrefracffo = zeros(param.nPyramidal,1);

Restffo = ExtractParameter(Feedforward,'Vrest',param);
Hyperffo = ExtractParameter(Feedforward,'Vhyper',param);
tauffo = ExtractParameter(Feedforward,'tau',param);
Threshffo = ExtractParameter(Feedforward,'FThresh',param);
tauAMPA1ffo = ExtractParameter(Feedforward,'tauAMPA1',param);
tauAMPA2ffo = ExtractParameter(Feedforward,'tauAMPA2',param);
EAMPAffo = ExtractParameter(Feedforward,'EAMPA',param);

% ...New parameters should be added here...

% Initialize Feedforward cells potentials
Vffo(:,1) = Restffo;
% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. tAMPA0ffo stores t0 delayed
% by LOT delay
tAMPA0ffo = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;
Imitffo = zeros(param.nPyramidal,1); % Input coming from Mitral cells
Imitffo_matrix = zeros(param.nPyramidal,round(param.tsim / param.dt));
maxgAMPAffo = getmaxg(param.dt,tauAMPA1ffo,tauAMPA2ffo); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Feedforward neurons parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Pyramidal cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Pyramidal cell at a given time
Vpyr = zeros(param.nPyramidal,round(param.tsim / param.dt));

% Binary matrix recording only the spikes
Spyr = Vpyr;
Cmitpyr = Vpyr;
Cpyrpyr = Vpyr;
Cfbapyr = Vpyr;
Cffopyr = Vpyr;

% Stores the AHP values over time
VAHPpyr = Vpyr;

refracpyr = 2; % Refractory period after firing (ms)
countrefracpyr = zeros(param.nPyramidal,1);

Restpyr = ExtractParameter(Pyramidal,'Vrest',param);
Hyperpyr = ExtractParameter(Pyramidal,'Vhyper',param);
taupyr = ExtractParameter(Pyramidal,'tau',param);
Threshpyr = ExtractParameter(Pyramidal,'FThresh',param);
tauAMPA1pyr = ExtractParameter(Pyramidal,'tauAMPA1',param);
tauAMPA2pyr = ExtractParameter(Pyramidal,'tauAMPA2',param);
EAMPApyr = ExtractParameter(Pyramidal,'EAMPA',param);
tauGABA1pyr = ExtractParameter(Pyramidal,'tauGABA1',param);
tauGABA2pyr = ExtractParameter(Pyramidal,'tauGABA2',param);
EGABApyr = ExtractParameter(Pyramidal,'EGABA',param);
tauAHPpyr = ExtractParameter(Pyramidal,'tauAHP',param);
EAHPpyr = ExtractParameter(Pyramidal,'EAHP',param);

% ...New parameters should be added here...

% Initialize Pyramidal cells potentials
Vpyr(:,1) = Restpyr;
Vpyr_nospike = Vpyr; % voltage without added spikes is used for calculating input currents

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. tAMPA0pyr is a matrix that
% will hold values of t0 delayed by LOT delay time
tAMPA0pyr = zeros(param.nMitral,round(param.tsim / param.dt)) - 10000000;
% generate vector of LOT delay variabilities from each mit neuron within time
% resolution of dt (in ms)
LOT_delay_vars = param.dt*ceil((2*param.LOTvar*rand(param.nMitral,1)-param.LOTvar)/param.dt);

Imitpyr = zeros(param.nPyramidal,1); % Input coming from Mitral cells
Imitpyr_matrix = zeros(param.nFeedback,round(param.tsim / param.dt));
maxgAMPApyr = getmaxg(param.dt,tauAMPA1pyr,tauAMPA2pyr); % Get max conductance
% amplitude

% AMPAPY time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPAPY0pyr = zeros(param.nPyramidal,round(param.tsim / param.dt)) - 10000000;
Ipyrpyr = zeros(param.nPyramidal,1); % Input coming from Pyramidal cells
Ipyrpyr_matrix = zeros(param.nFeedback,round(param.tsim / param.dt));

% GABA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Pyramidal cell can be
% connected to a varied number of Feedback cells
tGABA0pyr = zeros(param.nFeedback,round(param.tsim / param.dt)) - 10000000;
Ifbapyr = zeros(param.nPyramidal,1); % Input coming from Feedback cells
Ifbapyr_matrix = zeros(param.nFeedback,round(param.tsim / param.dt));
maxgGABApyr = getmaxg(param.dt,tauGABA1pyr,tauGABA2pyr); % Get max conductance
% amplitude

% GABAFF time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. Each Pyramidal cell can be
% connected to one Feedforward cell
tGABAFF0pyr = zeros(param.nPyramidal,round(param.tsim / param.dt)) - 10000000;
Iffopyr = zeros(param.nPyramidal,1); % Input coming from Feedforward cells
Iffopyr_matrix = zeros(param.nPyramidal,round(param.tsim / param.dt));

% AHP time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0. The AHP in Pyramidal cells
% depends on their own fire.
tAHP0pyr = zeros(param.nPyramidal,1) - 10000000;
% maxgAHPpyr = getmaxg(param.dt,tauCA1pyr,tauCA2pyr); % Get max conductance
% % amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Pyramidal cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Feedback cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stores the voltage of each Feedback cell at a given time
Vfba = zeros(param.nFeedback,round(param.tsim / param.dt));
% Binary matrix recording only the spikes
Sfba = Vfba;
Cpyrfba = Vfba;

refracfba = 2; % Refractory period after firing (ms)
countrefracfba = zeros(param.nFeedback,1);

Restfba = ExtractParameter(Feedback,'Vrest',param);
Hyperfba = ExtractParameter(Feedback,'Vhyper',param);
taufba = ExtractParameter(Feedback,'tau',param);
Threshfba = ExtractParameter(Feedback,'FThresh',param);
tauAMPA1fba = ExtractParameter(Feedback,'tauAMPA1',param);
tauAMPA2fba = ExtractParameter(Feedback,'tauAMPA2',param);
EAMPAfba = ExtractParameter(Feedback,'EAMPA',param);

% ...New parameters should be added here...

% Initialize Feedback cells potentials
Vfba(:,1) = Restfba;
Vfba_nospike = Vfba;

% AMPA time counter. This variable starts with a very negative value just
% to make sure that the currents will be = 0
tAMPA0fba = zeros(param.nPyramidal,round(param.tsim / param.dt)) - 10000000;
Ipyrfba = zeros(param.nFeedback,1); % Input coming from Pyramidal cells
Ipyrfba_matrix = zeros(param.nFeedback,round(param.tsim / param.dt));
maxgAMPAfba = getmaxg(param.dt,tauAMPA1fba,tauAMPA2fba); % Get max conductance
% amplitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Feedback cells parameters and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin neuron simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Start loop
for tt = 2:round(param.tsim / param.dt)
    t = tt * param.dt; % current time
    

    % Mitral Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Glo input to Mitral cells. This current is applied directly in
    % Vmit, no channel conductance.
%     Iglomit(:) = diag(wGloMit) .* Pfireglo(:,tt) .* gmaxAMPAmit; % glo - mit input

%     Iglomit(:) = diag(wGloMit) .* Pfireosn(:,tt) .* gmaxAMPAmit .* ... 
%                 (ones(param.nMitral,1) + param.Inoise*randn(param.nMitral,1)); % direct OSN - Mit input
%     
%     Iglomit_matrix(:,tt) = Iglomit(:);

% No OSN, just plain old direct inupt to Mit cells
Iglomit(:) = Respiration(tt) .* Iextosn(:,tt) .* (ones(param.nMitral,1) + param.Inoise*randn(param.nMitral,1));
Iglomit_matrix(:,tt) = Iglomit(:);
    
    
    % Do weights change with time? If yes scale by Wfrac
    if param.flagWGRAMITvar == true
        wGradistMit_scaled = param.Wfrac(tt-1) * wGradistMit;
        else
        wGradistMit_scaled = wGradistMit;
    end
    
    % Scale Gradist-Mat weights to simulate vesicle depletion due to Gra
    % spikes
%     Igs = Vgraprox(:,tt-1) == param.SpikeV; % index of gra cells that just spiked
%     VFRACgradist(Igs,tt-1) = 0;
%     VFRACgradist(:,tt) = VFRACgradist(:,tt-1) + (param.dt/Vtau)*VFRACgradist(:,tt-1);
%     wGradistMit_scaled = repmat(VFRACgradist(:,tt)',param.nMitral,1) .* wGradistMit;
    % NOTE: we want wGradistMit_scaled to have the same dimensions as
    % wGradistMit, namely [nMitral x nProximal]
            
    % Get Granule graded inputs to Mitral cells
    Igradistmit(:) = SetInoSpike_GraMit(gmaxGABAmit,Pfiregradist(:,tt-1),EGABAmit,Vmit_nospike(:,tt - 1),...
        wGradistMit_scaled);
    Igradistmit_matrix(:,tt) = Igradistmit(:);
    
    
    % Get AHP to Mitral cells
    if param.mitAHP == true
        VAHPmit(:,tt) = VAHPmit(:,tt - 1) + (param.dt / tauAHPmit(1)) *...
        (EAHPmit(1) * Smit(:,tt -1) - VAHPmit(:,tt - 1));    
    end
    
    Vnoisemit = param.noisemit .* randn(param.nMitral,1);
    
    % Mitral cell potential
    
    % Forwards Euler
    Vmit(:,tt) = Vmit(:,tt - 1) + (param.dt ./ taumit(:)) .* (Rmit(:) .*...
        (Iglomit(:) + Igradistmit(:)) - Vmit(:,tt - 1) + Restmit(:) + VAHPmit(:,tt) + Vnoisemit);
    
    % Backwards Euler
    
%     g_gramit = (((tauGABA1mit(1) * tauGABA2mit(1)) / (tauGABA1mit(1) - tauGABA2mit(1))) *...
%     (exp(-(t - tGABA0mit(:,tt)) / tauGABA1mit(1)) - exp(-(t - tGABA0mit(:,tt)) / tauGABA2mit(1)))) / maxgGABAmit(1);
% 
%     Vmit(:,tt) = (Vmit(:,tt - 1) + (param.dt ./ taumit(1)) .* (Rmit(:) .*(Iglomit(:) + (wGraMit_scaled * g_gramit) .* EGABAmit(1) + Restmit(:))))...
%         ./ (1 + (param.dt ./ taumit(1)) .* (Rmit(:) .* (wGraMit_scaled * g_gramit) + 1));
%     
    
    
%     AHPmit(:,tt) = Iahpmit(:);
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vmit(:,tt - 1) == param.SpikeV;
    Vmit(I,tt) = Hypermit(I);
    Vmit_nospike(:,tt) = Vmit(:,tt);
    
    % I is a vector of 1s or 0s of length ncells
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Granule variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This variable is used for input onto distal gra dendrites
        Mit_delayed_spike_time = t + param.Delay(4);
        tAMPA0mit_gra(I,ceil(Mit_delayed_spike_time/param.dt):end) = Mit_delayed_spike_time;
        tNMDA0mit_gra = tAMPA0mit_gra;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pyramidal variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This variable is set later in the code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if param.mitAHP == true
        tAHP0mit(I) = t; % new t0 for the Mitral cell AHP
    end
    
    countrefracmit(I) = refracmit / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracmit > 0;
    countrefracmit(I) = countrefracmit(I) - 1;
    
    
    I = find(countrefracmit == 0); % if countrefracmit = 0 the neuron can fire
    % again
    
%     if ~isempty(I)
%         if param.flagnoisemit == true
%             aux_J = SpikeNoise(Restmit(I),Threshmit(I),param,Vmit(I,tt),'Mitral');
%             J = find(aux_J);
%         else
%             J = find(Vmit(I,tt) >= Threshmit(I));
%         end
%         if ~isempty(J)
%             Vmit(I(J),tt) = param.SpikeV; % Action potential
%             Smit(I(J),tt) = 1; % Record spike time
%         end
%     end
    
    if ~isempty(I)
    J = find(Vmit(I,tt) >= Threshmit(I));
    if ~isempty(J)
            Vmit(I(J),tt) = param.SpikeV; % Action potential
            Smit(I(J),tt) = 1; % Record spike time
    end
    end
    
    
    % Granule Distal Dendrites
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Only evaluate if distal granule dendritic synapses are included
if param.DistalON == true
        
    % Get Mitral input to Distal Granule AMPA receptors
    ImitgradistAMPA(:) = SetI(tauAMPA1(1),tauAMPA2(1),t,tAMPA0mit_gra(:,tt),maxgAMPAgradist(1),...
        wMitGradistAMPA,EAMPAgradist(1),Vgradist(:,tt - 1));
    ImitgradistAMPA_matrix(:,tt) = ImitgradistAMPA(:);
    
    % Mitral input to Distal Granule NMDA receptors
    ImitgradistNMDA = SetI_NMDA(tauNMDA1(1),tauNMDA2(1),t,tNMDA0mit_gra(:,tt),...
        wMitGradistNMDA,ENMDAgradist(1),Vgradist(:,tt - 1));
    ImitgradistNMDA_matrix(:,tt) = ImitgradistNMDA(:);
    
    % Do Pyr-gra weights change with time? If yes scale by Wfrac
    if param.flagWOPvar == true
        wPyrGradist_scaled = param.Wfrac(tt-1) * wPyrGraprox;
        else
        wPyrGradist_scaled = wPyrGraprox;
    end
    
    % Get Proximal Granule or Pyramidal input to Distal Granule Dendrites depending on 
    % wether or not proximal exist.
    if param.ProximalON == true
        % Prox-dist current is applied directly in Vgradist, no ion channel conductance.
        % Sgraprox(:,tt-1) is [nGraprox x 1]
        % MatProxDist is [nGradist x nGraprox]
        % gmaxAMPAgradist is [nGradist x 1]
        % we want Iproxdist to be [nGradist x 1] so we need to transpose all 3 terms!
        Iproxdist(:) = Sgraprox(:,tt-1)'*MatProxDist' .* gmaxAMPAgradist';
        Iproxdist_matrix(:,tt) = Iproxdist(:);
    else
        % spiking input from Pyr
        if param.flagOP == true
            % use matrix of delayed tAMPA0pyr_gra values
            Ipyrgradist(:) = SetI(tauPROX1(1),tauPROX2(1),t,tAMPA0pyr_gra(:,tt),maxgAMPAgraprox(1),...
            wPyrGradist_scaled,EAMPAgraprox(1),Vgradist(:,tt - 1));
            Ipyrgradist_matrix(:,tt) = Ipyrgradist(:);
            % NOTE: Still using rise and fall time constants of proximal
            % dendrites for pyr input
        end
    end
    
    
    Vnoisegradist = param.noisegradist .* randn(param.nGradist,1);
    
    % Distal Granule cell potential
    % Forward Euler
    Vgradist(:,tt) = Vgradist(:,tt - 1) + (param.dt ./ taugradist(:)) .* (Rgradist(:) .*...
        (ImitgradistAMPA(:) + ImitgradistNMDA(:) + Iproxdist(:) + Ipyrgradist(:))...
        - Vgradist(:,tt - 1) + Restgradist(:) + Vnoisegradist);
    
    
    % "firing probability" for Graded inhibition from distal gra dendrites
    Pfiregradist(:,tt) = NSpikeP(Vgradist(:,tt),Restgradist,Threshgradist);
    
end
    
    % Granule Proximal Soma
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if param.ProximalON == true
    % Get bulbar input to granule soma
    
    
    if param.DistalON == true
        % Distal Granule dendrite input to Proximal Granule. This current is applied directly in
        % Vgraprox, no channel conductance (Nernst potetial not used).
        % Pfiregradist is [nGradist x 1]
        % MatProxDist is [nGradist x nGraprox]
        % gmaxAMPAgraprox is [nGraprox x 1]
        % we want Idistprox to be [nGraprox x 1] so we need to transpose only 1st and 3rd terms!
        Idistprox(:) = Pfiregradist(:,tt)'*MatProxDist .* gmaxAMPAgraprox';
%         Vdistprox(:) = (Vgradist(:,tt)'*MatProxDist)./nds; % Avg dendritic membrane potentials per soma
%         Idistprox(:) = Idistprox(:)./Rgraprox(:);
        Idistprox_matrix(:,tt) = Idistprox(:);
    else
        % if distal dendrites are removed, mit synapse directly onto gra
        Imitgraprox(:) = SetI(tauPROX1(1),tauPROX2(1),t,tAMPA0mit_gra(:,tt),maxgAMPAgraprox(1),...
        wMitGradistAMPA,EAMPAgraprox(1),Vgraprox_nospike(:,tt - 1));
        Imitgraprox_matrix(:,tt) = Imitgraprox(:);
    end
    
    
    % Do Pyr-gra weights change with time? If yes scale by Wfrac
    if param.flagWOPvar == true
        wPyrGraprox_scaled = param.Wfrac(tt-1) * wPyrGraprox;
        else
        wPyrGraprox_scaled = wPyrGraprox;
    end
    
    % Get Pyramidal input to Granule Proximal Dendrites
    if param.flagOP == true
        % use matrix of delayed tAMPA0pyr_gra values
        Ipyrgraprox(:) = SetI(tauPROX1(1),tauPROX2(1),t,tAMPA0pyr_gra(:,tt),maxgAMPAgraprox(1),...
        wPyrGraprox_scaled,EAMPAgraprox(1),Vgraprox_nospike(:,tt - 1));
        Ipyrgraprox_matrix(:,tt) = Ipyrgraprox(:);
    end
    
    % Get Granule input to Granule cells
    if param.GraGracon == true
        Igragra(:) = SetI(tauGABA1graprox(1),tauGABA2graprox(1),t,tGABA0graprox(:,tt),maxgGABAgraprox(1),...
            wProxProx,EGABAgraprox(1),Vgraprox_nospike(:,tt - 1));
        Igragra_matrix(:,tt) = Igragra(:);
    end
    
    % Get AHP to Proximal Granule cells
    if param.graAHP == true
        VAHPgraprox(:,tt) = VAHPgraprox(:,tt - 1) + (param.dt / tauAHPgraprox(1)) *...
        (EAHPgraprox(1) * Sgraprox(:,tt -1) - VAHPgraprox(:,tt - 1));    
    end
    
    % Get LLD to Proximal Granule cells
    if param.graLLD == true
        VLLDgraprox(:,tt) = VLLDgraprox(:,tt - 1) + (param.dt / tauLLDgraprox(1)) *...
        (ELLDgraprox(1) * Sgraprox(:,tt -1) - VLLDgraprox(:,tt - 1));    
    end
    
    Vnoisegraprox = param.noisegraprox .* randn(param.nGraprox,1);
    
    % Proximal Granule cell potential
    % Forward Euler
    Vgraprox(:,tt) = Vgraprox(:,tt - 1) + (param.dt ./ taugraprox(:)) .* (Rgraprox(:) .*...
        (Igragra(:) + Idistprox(:) + Imitgraprox(:) + Ipyrgraprox(:))...
        - Vgraprox(:,tt - 1) + Restgraprox(:) + VAHPgraprox(:,tt) + VLLDgraprox(:,tt) + Vnoisegraprox);

    % Backwards Euler
%     g_mitgra = (((tauPROX1(1) * tauPROX2(1)) / (tauPROX1(1) - tauPROX2(1))) *...
%     (exp(-(t - tAMPA0mit_gra(:,tt)) / tauPROX1(1)) - exp(-(t - tAMPA0mit_gra(:,tt)) / tauPROX2(1)))) / maxgAMPAgraprox(1);
% 
%     g_pyrgra = (((tauPROX1(1) * tauPROX2(1)) / (tauPROX1(1) - tauPROX2(1))) *...
%     (exp(-(t - tAMPA0pyr_gra(:,tt)) / tauPROX1(1)) - exp(-(t - tAMPA0pyr_gra(:,tt)) / tauPROX2(1)))) / maxgAMPAgraprox(1);
% 
%     Vgraprox(:,tt) = (Vgraprox(:,tt - 1) + (param.dt ./ taugraprox(1)) .* (Rgraprox(:) .*((wMitGradist * g_mitgra) .* EAMPAgraprox(1) + (wPyrGra_scaled * g_pyrgra) .* EAMPAgraprox(1) + Restgraprox(:))))...
%         ./ (1 + (param.dt ./ taugraprox(1)) .* (Rgraprox(:) .* (wMitGradist * g_mitgra + wPyrGra_scaled * g_pyrgra) + 1));

    
    % NOTE Gra-Gra and AHP are ignored in backwards euler for now
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vgraprox(:,tt - 1) == param.SpikeV;
    Vgraprox(I,tt) = Hypergraprox(I);
    Vgraprox_nospike(:,tt) = Vgraprox(:,tt);
    
    
    % Proximal Granule Dendrite variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GraGra_delayed_spike_time = t + param.Delay(5); % new t0 for the Distal Granule cell GABA current
    tGABA0graprox(I,ceil(GraGra_delayed_spike_time/param.dt):end) = GraGra_delayed_spike_time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mitral Cell Variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE!!! For the 2CG model tGABA0mit is replicated for each of the nds
    % dendrites belonging to each gra cell soma, so it has dimensions
    % [nGradist x dt*tsim]
        Idist = zeros(param.nGradist,1);
        for ii = 1:param.nGraprox
        Idist((nds*(ii-1)+1):nds*ii) = I(ii);
        end
        % new t0 for the Mitral cell GABA current is delayed by MitGABAdelay
        Gra_delayed_spike_time = t + param.Delay(3);
        tGABA0mit(logical(Idist),ceil(Gra_delayed_spike_time/param.dt):end) = Gra_delayed_spike_time;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if param.graAHP == true
        tAHP0graprox(I) = t; % new t0 for the Proximal Granule Dendrite AHP   
    end
    if param.graLLD == true
        tLLD0graprox(I) = t; % new t0 for the Proximal Granule Dendrite LLD   
    end
    
    countrefracgraprox(I) = refracgraprox / param.dt; % neurons that just fired get
    % into the refractory period
    I = countrefracgraprox > 0;
    countrefracgraprox(I) = countrefracgraprox(I) - 1;
    
    I = find(countrefracgraprox == 0); % if countrefracgra = 0 the neuron can fire
    % again
    
%     if ~isempty(I)
%         if param.flagnoisegraprox == true
%             aux_J = SpikeNoise(Restgraprox(I),Threshgraprox(I),param,Vgraprox(I,tt),'GraProximal');
%             J = find(aux_J);
%         else
%             J = find(Vgraprox(I,tt) >= Threshgraprox(I));
%         end
%         if ~isempty(J)
%             Vgraprox(I(J),tt) = param.SpikeV; % Action potential
%             Sgraprox(I(J),tt) = 1; % Record spike time
%         end
%     end

    if ~isempty(I)
        J = find(Vgraprox(I,tt) >= Threshgraprox(I));
        
        if ~isempty(J)
            Vgraprox(I(J),tt) = param.SpikeV; % Action potential
            Sgraprox(I(J),tt) = 1; % Record spike time
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PC neurons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% get t0 variables for Pyr and FF neurons
I = Vmit(:,tt - 1) == param.SpikeV;
    
    
    if param.flagLOT == true
        % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Ivals = 1:param.nMitral;
        Ivals = Ivals(I); % Need actual values of indices for spiking neurons. I is just [1 0 1 0 0 1...]
        LOT_delayed_spike_times = t + param.Delay(1) + LOT_delay_vars; % LOT delay is in ms
        % NOTE: each LOT fiber has it's own latency so I just loop
                % through the spiking neurons to set latencies individually
        for ii = Ivals
            tAMPA0pyr(ii,ceil(LOT_delayed_spike_times(ii)/param.dt):end) = LOT_delayed_spike_times(ii);
        end

        % Feedforward cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if param.FfoNeurons == true
            Ivals = 1:param.nMitral;
            Ivals = Ivals(I); % Need actual values of indices for spiking neurons. I is just [1 0 1 0 0 1...]
            LOT_delayed_spike_times = t + param.Delay(1) + LOT_delay_vars; % LOT delay is in ms
            % NOTE: each LOT fiber has it's own latency so I just loop
                % through the spiking neurons to set latencies individually
            for ii = Ivals
                tAMPA0ffo(ii,ceil(LOT_delayed_spike_times(ii)/param.dt):end) = LOT_delayed_spike_times(ii);
            end
        end
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    % Feedforward Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Mitral input to Feedforward cells
    
    % Feedforward neurons are only simulated if FfoNeurons = true
    if param.FfoNeurons == true
        
    % Do weights change with time? If yes scale by Wfrac
    if param.flagWLOTvar == true
        wMitFfo_scaled = param.Wfrac(tt-1) * wMitFfo;
        else
        wMitFfo_scaled = wMitFfo;
    end
    
    if param.flagLOT == true
    Imitffo(:) = SetI(tauAMPA1ffo(1),tauAMPA2ffo(1),t,tAMPA0ffo(:,tt),maxgAMPAffo(1),...
        wMitFfo_scaled,EAMPAffo(1),Vffo(:,tt - 1));
    Imitffo_matrix(:,tt) = Imitffo(:);
    end

    % Feedforward cell potential
    Vffo(:,tt) = Vffo(:,tt - 1) + (param.dt ./ tauffo(:)) .* (Imitffo(:) - ...
        Vffo(:,tt - 1) + Restffo(:));
    
    Cmitffo(:,tt) = Imitffo(:);
    
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vffo(:,tt - 1) == param.SpikeV;
    Vffo(I,tt) = Hyperffo(I);
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FF_delayed_spike_time = t + param.Delay(6);
    tGABAFF0pyr(I,ceil(FF_delayed_spike_time/param.dt):end) = FF_delayed_spike_time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracffo(I) = refracffo / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracffo > 0;
    countrefracffo(I) = countrefracffo(I) - 1;
    
    I = find(countrefracffo == 0); % if countrefracffo = 0 the neuron can fire
    % again
    if ~isempty(I)
        if param.flagnoiseffo == true
            aux_J = SpikeNoise(Restffo(I),Threshffo(I),param,Vffo(I,tt),'Feedforward');
            J = find(aux_J);
        else
            J = find(Vffo(I,tt) >= Threshffo(I));
        end
        if ~isempty(J)
            Vffo(I(J),tt) = param.SpikeV; % Action potential
            Sffo(I(J),tt) = 1; % Record spike time
        end
    end
    
    end
    
    % Pyramidal Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Mitral input to Pyramidal cells
    
    % Do weights change with time? If yes scale by Wfrac
    if param.flagWLOTvar == true
        wMitPyr_scaled = param.Wfrac(tt-1) * wMitPyr;
        else
        wMitPyr_scaled = wMitPyr;
    end
    
    if param.flagLOT == true
        % delayed conductance uses tAMPA0pyr matrix that stores all delayed
        % t0 values
        Imitpyr(:) = SetI(tauAMPA1pyr(1),tauAMPA2pyr(1),t,tAMPA0pyr(:,tt),maxgAMPApyr(1),...
        wMitPyr_scaled,EAMPApyr(1),Vpyr_nospike(:,tt - 1));
        Imitpyr_matrix(:,tt) = Imitpyr(:);
    end
    
    % Get Pyramidal input to Pyramidal cells
    Ipyrpyr(:) = SetI(tauAMPA1pyr(1),tauAMPA2pyr(1),t,tAMPAPY0pyr(:,tt),maxgAMPApyr(1),...
        wPyrPyr,EAMPApyr(1),Vpyr_nospike(:,tt - 1));
    Ipyrpyr_matrix(:,tt) = Ipyrpyr(:);
    
    % Get Feedback input to Pyramidal cells
    Ifbapyr(:) = SetI(tauGABA1pyr(1),tauGABA2pyr(1),t,tGABA0pyr(:,tt),maxgGABApyr(1),...
        wFbaPyr,EGABApyr(1),Vpyr_nospike(:,tt - 1));
    Ifbapyr_matrix(:,tt) = Ifbapyr(:);
    
    % Get Feedforward input to Pyramidal cells
    Iffopyr(:) = SetI(tauGABA1pyr(1),tauGABA2pyr(1),t,tGABAFF0pyr(:,tt),maxgGABApyr(1),...
        wFfoPyr,EGABApyr(1),Vpyr_nospike(:,tt - 1));
    Iffopyr_matrix(:,tt) = Iffopyr(:);
    
    if param.pyrAHP == true
        VAHPpyr(:,tt) = VAHPpyr(:,tt - 1) + (param.dt / tauAHPpyr(1)) *...
            (EAHPpyr(1) * Spyr(:,tt -1) - VAHPpyr(:,tt - 1));
    end

    
    Vnoisepyr = param.noisepyr .* randn(param.nPyramidal,1);
    
    % Pyramidal cell potential
    % Forward Euler
    Vpyr(:,tt) = Vpyr(:,tt - 1) + (param.dt ./ taupyr(:)) .* ((Imitpyr(:)...
         + Ipyrpyr(:) + Ifbapyr(:) + Iffopyr(:)) - Vpyr(:,tt - 1)...
         + Restpyr(:) + VAHPpyr(:,tt) + Vnoisepyr);
     
%   Backwards Euler
%     g_mitpyr = (((tauAMPA1pyr(1) * tauAMPA2pyr(1)) / (tauAMPA1pyr(1) - tauAMPA2pyr(1))) *...
%     (exp(-(t - tAMPA0pyr(:,tt)) / tauAMPA1pyr(1)) - exp(-(t - tAMPA0pyr(:,tt)) / tauAMPA2pyr(1)))) / maxgAMPApyr(1);
% 
%     g_pyrpyr = (((tauAMPA1pyr(1) * tauAMPA2pyr(1)) / (tauAMPA1pyr(1) - tauAMPA2pyr(1))) *...
%     (exp(-(t - tAMPA0pyr(:,tt)) / tauAMPA1pyr(1)) - exp(-(t - tAMPA0pyr(:,tt)) / tauAMPA2pyr(1)))) / maxgAMPApyr(1);
% 
%     g_fbapyr = (((tauGABA1pyr(1) * tauGABA2pyr(1)) / (tauGABA1pyr(1) - tauGABA2pyr(1))) *...
%     (exp(-(t - tGABA0pyr(:,tt)) / tauGABA1pyr(1)) - exp(-(t - tGABA0pyr(:,tt)) / tauGABA2pyr(1)))) / maxgGABApyr(1);
% 
%     Vpyr(:,tt) = (Vpyr(:,tt - 1) + (param.dt ./ taupyr(1)) .* ...
%         ((wMitPyr_scaled * g_mitpyr) .* EAMPApyr(1) + (wPyrPyr * g_pyrpyr) .* EAMPApyr(1) + (wFbaPyr * g_fbapyr) .* EGABApyr(1) + Restpyr(:)))...
%         ./ (1 + (param.dt ./ taupyr(1)) .* ((wMitPyr_scaled * g_mitpyr + wPyrPyr * g_pyrpyr + wFbaPyr * g_fbapyr) + 1));
%     


     Cmitpyr(:,tt) = Imitpyr(:);
     Cpyrpyr(:,tt) = Ipyrpyr(:);
     Cfbapyr(:,tt) = Ifbapyr(:);
     Cffopyr(:,tt) = Iffopyr(:);
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vpyr(:,tt - 1) == param.SpikeV;
    Vpyr(I,tt) = Hyperpyr(I);
    Vpyr_nospike(:,tt) = Vpyr(:,tt);
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new t0 for the Pyramidal cell AMPA current
    PyrPyr_delayed_spike_time = t + param.Delay(7);
    tAMPAPY0pyr(I,ceil(PyrPyr_delayed_spike_time/param.dt):end) = PyrPyr_delayed_spike_time;
    
    tAHP0pyr(I) = t; % new t0 for the Pyramidal cell AHP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Feedback cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new t0 for the Feedback cell AMPA current
    PyrFba_delayed_spike_time = t + param.Delay(8);
    tAMPA0fba(I,ceil(PyrFba_delayed_spike_time/param.dt):end) = PyrFba_delayed_spike_time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Granule variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This variable is used for input onto either proximal or distal gra
    % dendrites, depending on which of them is present
    if param.flagOP == true
        Ivals = 1:param.nPyramidal;
        Ivals = Ivals(I); % Need actual values of indices for spiking neurons. I is just [1 0 1 0 0 1...]
            OP_delayed_spike_times = t + param.Delay(2) + OP_delay_vars; % OP delay is in ms
            % NOTE: each OP fiber has it's own latency so I just loop
                % through the spiking neurons to set latencies individually
            for ii = Ivals
                tAMPA0pyr_gra(ii,ceil(OP_delayed_spike_times(ii)/param.dt):end) = OP_delayed_spike_times(ii);
            end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracpyr(I) = refracpyr / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracpyr > 0;
    countrefracpyr(I) = countrefracpyr(I) - 1;
    
    I = find(countrefracpyr == 0); % if countrefracpyr = 0 the neuron can fire
    % again
%     if ~isempty(I)
%         if param.noisepyr == true
%             aux_J = SpikeNoise(Restpyr(I),Threshpyr(I),param,Vpyr(I,tt),'Pyramidal');
%             J = find(aux_J);
%         else
%             J = find(Vpyr(I,tt) >= Threshpyr(I));
%         end
%         if ~isempty(J)
%             Vpyr(I(J),tt) = param.SpikeV; % Action potential
%             Spyr(I(J),tt) = 1; % Record spike time
%         end
%     end
    if ~isempty(I)
        J = find(Vpyr(I,tt) >= Threshpyr(I));

        if ~isempty(J)
            Vpyr(I(J),tt) = param.SpikeV; % Action potential
            Spyr(I(J),tt) = 1; % Record spike time
        end
    end
    
    
    % Feedback Cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get Pyramidal input to Feedback cells
    Ipyrfba(:) = SetI(tauAMPA1fba(1),tauAMPA2fba(1),t,tAMPA0fba(:,tt),maxgAMPAfba(1),...
        wPyrFba,EAMPAfba(1),Vfba_nospike(:,tt - 1));
    Ipyrfba_matrix(:,tt) = Ipyrfba(:);

    Vnoisefba = param.noisefba .* randn(param.nFeedback,1);
    
    % Feedback cell potential
    % Forward Eueler
    Vfba(:,tt) = Vfba(:,tt - 1) + (param.dt ./ taufba(:)) .*...
        (Ipyrfba(:) - Vfba(:,tt - 1) + Restfba(:) + Vnoisefba);
    
    % Backwards Euler
%     g_pyrfba = (((tauAMPA1fba(1) * tauAMPA2fba(1)) / (tauAMPA1fba(1) - tauAMPA2fba(1))) *...
%     (exp(-(t - tAMPA0fba(:,tt)) / tauAMPA1fba(1)) - exp(-(t - tAMPA0fba(:,tt)) / tauAMPA2fba(1)))) / maxgAMPAfba(1);
% 
%     Vfba(:,tt) = (Vfba(:,tt - 1) + (param.dt ./ taufba(1)) .* ((wPyrFba * g_pyrfba) .* EAMPAfba(1) + Restfba(:)))...
%         ./ (1 + (param.dt ./ taupyr(1)) .* ((wPyrFba * g_pyrfba) + 1));
%     
%     
%     Cpyrfba(:,tt) = Ipyrfba(:);
%     
    
    % If the neuron fired last cycle, neuron potential hyperpotentializes
    I = Vfba(:,tt - 1) == param.SpikeV;
    Vfba(I,tt) = Hyperfba(I);
    Vfba_nospike(:,tt) = Vfba(:,tt);
    
    % Pyramidal cell variable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % new t0 for the Pyramidal cell GABA current
    Fba_delayed_spike_time = t + param.Delay(6);
    tGABA0pyr(I,ceil(Fba_delayed_spike_time/param.dt):end) = Fba_delayed_spike_time;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    countrefracfba(I) = refracfba / param.dt; % neurons that just fired get into the
    % refractory period
    I = countrefracfba > 0;
    countrefracfba(I) = countrefracfba(I) - 1;
    
    I = find(countrefracfba == 0); % if countrefracfba = 0 the neuron can fire
    % again
%     if ~isempty(I)
%         if param.noisefba == true
%             aux_J = SpikeNoise(Restfba(I),Threshfba(I),param,Vfba(I,tt),'Feedback');
%             J = find(aux_J);
%         else
%             J = find(Vfba(I,tt) >= Threshfba(I));
%         end
%         if ~isempty(J)
%             Vfba(I(J),tt) = param.SpikeV; % Action potential
%             Sfba(I(J),tt) = 1; % Record spike time
%         end
%     end
    
if ~isempty(I)
        J = find(Vfba(I,tt) >= Threshfba(I));
        if ~isempty(J)
            Vfba(I(J),tt) = param.SpikeV; % Action potential
            Sfba(I(J),tt) = 1; % Record spike time
        end
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OSN and Mitral cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nMitral
    Mitral{ii}.V = Vmit(ii,:); % Save neuronal activity
    Mitral{ii}.VAHP = VAHPmit(ii,:);
    Mitral{ii}.S = Smit(ii,:); % Save spike time
end

% Granule cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nGraprox
    if param.ProximalON == true
        GraProximal{ii}.V = Vgraprox(ii,:); % Save neuronal activity
        GraProximal{ii}.S = Sgraprox(ii,:); % Save spike time
    else
        GraProximal{ii}.V = zeros(1,round(param.tsim / param.dt)); % 0 vector
        GraProximal{ii}.S = zeros(1,round(param.tsim / param.dt)); % 0 vector
    end
end
for ii = 1:param.nGradist
    if param.DistalON == true
        GraDistal{ii}.V = Vgradist(ii,:); % Save neuronal activity
        GraDistal{ii}.S = Sgradist(ii,:); % Save spike time
    else
        GraDistal{ii}.V = zeros(1,round(param.tsim / param.dt)); % 0 vector
        GraDistal{ii}.S = zeros(1,round(param.tsim / param.dt)); % 0 vector
    end
end

% Save OSN file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- OSN removed

% Feedforward and Pyramidal cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nPyramidal
    Feedforward{ii}.V = Vffo(ii,:); % Save neuronal activity
    Feedforward{ii}.S = Sffo(ii,:); % Save spike time
    Feedforward{ii}.Cmit = Cmitffo(ii,:);
    Pyramidal{ii}.V = Vpyr(ii,:); % Save neuronal activity
    Pyramidal{ii}.AHP = VAHPpyr(ii,:); % Save AHP activity
    Pyramidal{ii}.S = Spyr(ii,:); % Save spike time
    Pyramidal{ii}.Cmit = Cmitpyr(ii,:);
    Pyramidal{ii}.Cffo = Cffopyr(ii,:);
    Pyramidal{ii}.Cpyr = Cpyrpyr(ii,:);
    Pyramidal{ii}.Cfba = Cfbapyr(ii,:);
    Feedforward{ii}.MitCon = MatMitFfo(ii,:); % Save Mitral connections
    Pyramidal{ii}.MitCon = MatMitPyr(ii,:); % Save Mitral connections
    Pyramidal{ii}.FfoCon = MatFfoPyr(ii,:); % Save Feedforward connections
    Pyramidal{ii}.FbaCon = MatFbaPyr(ii,:); % Save Feedback connections
    Pyramidal{ii}.PyrCon = MatPyrPyr(ii,:); % Save Pyramidal connections
end

% Feedback cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:param.nFeedback
    Feedback{ii}.V = Vfba(ii,:); % Save neuronal activity
    Feedback{ii}.S = Sfba(ii,:); % Save spike time
    Feedback{ii}.Cpyr = Cpyrfba(ii,:);
    Feedback{ii}.PyrCon = MatPyrFba(ii,:); % Save Pyramidal connections
end

% Save OB and PC connections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Savefile = strcat(param.outputPath,'ConnData');
    save(Savefile,'MatGradistMit','MatProxProx','MatProxDist','-append');
    Savefile = strcat(param.outputPath,'LearnData'); % learned weights
    save(Savefile,'wGradistMit','wProxProx','wGloMit','-append');

    Savefile = strcat(param.outputPath,'ConnData');
    save(Savefile,'MatPyrFba','MatFbaPyr','MatFfoPyr','MatPyrPyr','-append');
    Savefile = strcat(param.outputPath,'LearnData'); % learned weights
    save(Savefile,'wPyrPyr','wFbaPyr','wFfoPyr','wPyrFba','-append');


    Savefile = strcat(param.outputPath,'ConnData');
    save(Savefile,'MatPyrGraprox','MatMitPyr','MatMitFfo','-append');
    Savefile = strcat(param.outputPath,'LearnData'); % learned weights
    save(Savefile,'wPyrGraprox','wMitPyr','wMitFfo','-append');


% Save input currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to proximal granule
if param.ProximalON == true
    if param.DistalON == true
        InputCurrent.Idistprox = Idistprox_matrix;
    else
        InputCurrent.Imitgraprox = Imitgraprox_matrix;
    end
    InputCurrent.Ipyrgraprox = Ipyrgraprox_matrix;
    InputCurrent.Igragra = Igragra_matrix;
end

% to distal granule
if param.DistalON == true
    if param.ProximalON == true
        InputCurrent.Iproxdist = Iproxdist_matrix;
    else
        InputCurrent.Ipyrgradist = Ipyrgradist_matrix;
    end
    InputCurrent.ImitgradistAMPA = ImitgradistAMPA_matrix;
    InputCurrent.ImitgradistNMDA = ImitgradistNMDA_matrix;
end

% to mitral
if param.DistalON == true
    InputCurrent.Igradistmit = Igradistmit_matrix;
else
    InputCurrent.Igraproxmit = Igraproxmit_matrix;
end
InputCurrent.Iglomit = Iglomit_matrix;

% to pyramidal
InputCurrent.Imitpyr = Imitpyr_matrix;
InputCurrent.Ipyrpyr = Ipyrpyr_matrix;
InputCurrent.Ifbapyr = Ifbapyr_matrix;
InputCurrent.Iffopyr = Iffopyr_matrix;

% to feedworward
InputCurrent.Imitffo = Imitffo_matrix;

% to feedback
InputCurrent.Ipyrfba = Ipyrfba_matrix;

InputCurrent.Pfiregradist = Pfiregradist;
InputCurrent.Pfiregraprox = Pfiregraprox;

end



function Iext = SetExtInput(param,Neuron)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the inputs to OSN neurons
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/02/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Iext = zeros(param.nMitral,round(param.tsim / param.dt));

for ii = 1:param.nMitral
    Iext(ii,round(param.tinit / param.dt) : round(param.tfinal / param.dt)) = Neuron{ii}.input;
end

Iext = Iext * param.Iext;

end

function R = CreateRespFreq(param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function creates an artificial respiration to modulate the bulb
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 04/06/2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = param.dt:param.dt:param.tsim;
time = time / 1000; %converting the time to seconds
R = -cos(2 * pi * param.RespFreq * time);
R = R + 1;
R = R / 2;
end

function maxg = getmaxg(dt,tau1,tau2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function finds the max amplitude for g, so we can use this value to
% normalize the curve
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 11/10/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 100;

x = 0:dt:tmax;
maxg = zeros(length(tau1),1);

for ii = 1: length(tau1)
    if ii == 1
        y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii))) * (exp(-(x)...
            / tau1(ii)) - exp(-(x) / tau2(ii)));
        maxg(ii) = max(y);
    else
        if tau1(ii) == tau1(ii -1) && tau2(ii) == tau2(ii -1)
            maxg(ii) = maxg(ii - 1);
        else
            y = ((tau1(ii) * tau2(ii)) / (tau1(ii) - tau2(ii)))...
                * (exp(-(x) / tau1(ii)) - exp(-(x) / tau2(ii)));
            maxg(ii) = max(y);
        end
    end
end
end

function Ic = SetInoSpike(gmax,P,E,V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns the nonspike pglo input for a given excitation
%
% Licurgo de Almeida
% 08/04/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    g = P .* gmax;
    Ic = g .* (E - V);
end

function Ic = SetInoSpike_GraMit(gmax,P,E,V,W)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function returns the distal granule summed graded input to mitral cells
% (weight matrix W is not necessarily square)
%
%
% Modifed by Boleszek Osinski on 07/16/2013 to allow for different numbers
% of neurons (particularly for nonspiking input from Gradist to Mitral)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: if W  = wGraMit then it has dimensions nMitral x nGranule

%     g = P .* gmax;
Ic = zeros(size(W,1),1);
    for ii = 1:size(W,1)
    Ic(ii) = sum(W(ii,:)' .* P * gmax(ii) * (E(ii) - V(ii)));
    end
    
end

function INMDA = SetI_NMDA(tau1,tau2,t,t0,W,E,V)
% NMDA Synapse
g_norm = 1;
Mg_conc = 1;
gamma = 0.004;
eta = 0.33;
E_Mg = -70; % assuming that Vrest = -70

% t1 = 60;
% t2 = 10;
% 
% t = 0:1:200;
% 
% v = (E_Mg-0.1):0.001:(E_Mg+0.1);
% 
% g = zeros(length(t), length(v));
% for i = 1:length(t)
%     for j = 1:length(v)
%         g(i,j) = g_norm * (exp(-t(i)/t1) - exp(-t(i)/t2))/(1+eta*Mg_conc*exp(-(v(j)-E_Mg)/gamma));
%     end
% end
% figure(4321)
% surf(v,t,g, 'FaceColor', 'interp', 'edgecolor', 'none', 'FaceLighting', 'phong')
% colorbar;
% camlight left;
% xlabel('V_M (V)');ylabel('t (ms)');zlabel('Conductance')
% 
% % playing with gamma
% gamma = 0.005;
% Mg_block = 1./(1+eta*Mg_conc*exp(-(v-E_Mg)/gamma));
% plot(v,Mg_block)



% % variables for debugging
% tau1 = tauNMDA1(1);
% tau2 = tauNMDA2(1);
% t0 = tNMDA0mit_gra(:,tt);
% W = wMitGradist;
% E = ENMDAgradist(1);
% V = Vgradist(:,tt - 1);

% NOTE: If W is Wmitgradist then it has dimensions ngradist x nmit

g = g_norm * (exp(-(t - t0)/tau1) - exp(-(t - t0)/tau2)); % dim 100 x 1
Mg_block = 1./(1+eta*Mg_conc*exp(-(V-E_Mg)/gamma)); % dim 800 x 1
INMDA = (W * g) .* (V-E) .* Mg_block;

% NOTE!!! All voltages are in V, not mV, so parameters such as gamma have
% to be scaled accordingly
end


function I = SpikeNoise(Rest,Thres,param,V,tnet)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the spiking chance
% The main function of this program is NeuroActivity.m
%
% Licurgo de Almeida
% 12/20/2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(tnet,'GraProximal')
    Vlimit = 0e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 3;
    chance = 1;
elseif strcmp(tnet,'GraDistal')
    Vlimit = 0.1e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 2;
    chance = 1;
elseif strcmp(tnet,'Mitral')
    Vlimit = 0.1e-3; % Limits where firing chances varied between 0 and 1
    pfactor = 3;
    chance = 1;
elseif strcmp(tnet,'Feedforward')
    Vlimit = 0e-3;
    pfactor = 10;
    chance = 1;
elseif strcmp(tnet,'Pyramidal')
    Vlimit = 0e-3;
    pfactor = 7;
    chance = 1;
elseif strcmp(tnet,'Feedback')
    Vlimit = 0e-3;
    pfactor = 3;
    chance = 1;
else
    pfactor = 5;
end


Rest = Rest - Vlimit;
%Thres = Thres + Vlimit;
y = (V - Rest) ./ (Thres -  Rest);
J = y <= 0;
y(J) = 0;
J = y > 1;
y(J) = 1;
y = y.^pfactor;
I = rand(length(Rest),1) <= y * chance;

end

function P = NSpikeP(V,R,T)
% matrix of probabilities forced to be within range 0 - 1
P = (V - R) ./ (T - R);
J = P <= 0;
P(J) = 0;
J = P > 1;
P(J) = 1;

end


function Ic = SetI(tau1,tau2,t,t0,normg,W,E,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function calculates the currents at each time step.
% This function is called from NeuroActivity.m
%
% Licurgo de Almeida
% 11/18/2010
%
% Ic: channel current
% gmax: maximum conductance
% tau1: channel's rising time
% tau2: channel's falling time
% t: step time
% t0: last time the pre synaptic neuron fired
% normg: normalizes the conductance curve between 0 and 1
% W: synaptic weights
% E: Nernst potential
% V: neuron's potential from last timestep
% Mcon: connection matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g = (((tau1 * tau2) / (tau1 - tau2)) *...
    (exp(-(t - t0) / tau1) - exp(-(t - t0) / tau2))) / normg;

Ic = (W * g) .* (E - V);

end



function Mat = SetConnections(cell1,cell2,cchance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the connection between different neurons
% The main function of this program is piriformmain.m
%
% Modified by Boleszek Osinski on 07/15/2013 from the original
%
% Licurgo de Almeida
% 03/01/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: this function only sets fixed connections
    
        Mat = zeros(cell1,cell2);
            for ii = 1:cell1
                convector = randperm(cell2);
                convector = convector(1:round(cell2 * cchance));
                Mat(ii,convector) = 1;
            end

end

function w = SetWeights(Mat,weights,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the synaptic weights between different neurons
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/02/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = zeros(size(Mat));
for ii = 1:length(weights)
    w(ii,:) = Mat(ii,:) * weights(ii);
end

if param.NoiseParam == true
    for ii = 1:length(weights)
        signvar = rand(1,size(Mat,2));
        I = signvar < 0.5;
        signvar(I) = -1;
        I = signvar >= 0.5;
        signvar(I) = 1;
        
        changeval = (rand(1,size(Mat,2)) * param.NoiseLevel * weights(ii)) .* Mat(ii,:);
        w(ii,:) = w(ii,:) + (signvar .* changeval);
    end
end
end

function parameter = ExtractParameter(neuron,ParName,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function extract the parameters from the neuronfile
% The main function of this program is piriformmain.m
%
% Licurgo de Almeida
% 03/01/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameter = zeros(length(neuron),1);

switch lower(ParName)
    case 'vrest'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.Vrest;
        end
    case 'vhyper'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.Vhyper;
        end
    case 'r'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.R;
        end
    case 'tau'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tau;
        end
    case 'tauahp'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauAHP;
        end
    case 'fthresh'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.FThresh;
        end
    case 'gmaxampa'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxAMPA;
        end
    case 'gmaxampapy'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxAMPAPY;
        end
    case 'gmaxgaba'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxGABA;
        end
    case 'gmaxgabaff'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.gmaxGABAFF;
        end
    case 'tauampa1'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauAMPA1;
        end
    case 'tauampa2'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauAMPA2;
        end
    case 'taugaba1'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauGABA1;
        end
    case 'taugaba2'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauGABA2;
        end
    case 'tauca1'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauCA1;
        end
    case 'tauca2'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.tauCA2;
        end
    case 'eampa'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.EAMPA;
        end
    case 'egaba'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.EGABA;
        end
    case 'eahp'
        for ii = 1:length(neuron)
            parameter(ii) = neuron{ii}.EAHP;
        end
end

if param.NoiseParam == true
    for ii = 1:length(parameter)
        signvar = rand;
        if signvar < 0.5
            signivar = -1;
        else
            signivar = 1;
        end
        parameter(ii) = parameter(ii) + (signivar * parameter(ii) * rand* param.NoiseLevel);
    end
end
end

