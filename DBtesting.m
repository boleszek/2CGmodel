
% load DB
DBload


%% Save parameters in DB by category


% Define network parameter names
InputParams_cols = {'tsim','Iext','Inoise','RespON','RespFQ','randfrac'};

OBParams_cols = {'NMC','NdGC','NpGC','NoiseMC','NoisedGC','NoisepGC'...
    ,'pGCLLD','CMG'};

PCParams_cols = {'NP','NFB','NoiseP','NoiseFB','CPP','CPF','CFP'};

OBPCParams_cols = {'flagLOT','LOTvar','flagOP','OPvar','CMP','CPG'};

% Define cell parameter names
MCParams_cols = {'Vrest','Fthresh','WGM'};
dGCParams_cols = {'tNMDA1','tNMDA2','gmax','Vrest','Fthresh','WMG'};
pGCParams_cols = {'tLLD','ELLD','gmax','Vrest','Fthresh','WPG'};


% Extract parameters values
input_file = 'noisy_OB_PC_params_2CG.txt';
[names, vals] = textread(input_file,'%s %s',162);

inds = [3 6 7 8 9 10];
InputParams = cell(1,length(inds));
for ii = 1:length(inds)
    InputParams(ii) = {vals{inds(ii)}};
    % for some wierd reason I had to put {} on the rhs of = so that they
    % would show up as comumns in DB
end

inds = [14 15 16 19 20 21 24 30];
OBParams = cell(1,length(inds));
for ii = 1:length(inds)
    OBParams(ii) = {vals{inds(ii)}};
end

inds = [32 33 36 37 39 40 41];
PCParams = cell(1,length(inds));
for ii = 1:length(inds)
    PCNetworkParams(ii) = {vals{inds(ii)}};
end

inds = [43 44 45 46 49 50];
OBPCParams = cell(1,length(inds));
for ii = 1:length(inds)
    OBPCNetworkParams(ii) = {vals{inds(ii)}};
end

inds = [67 69 71];
MCParams = cell(1,length(inds));
for ii = 1:length(inds)
    MCParams(ii) = {vals{inds(ii)}};
end

inds = [77 78 81 82 84 86];
dGCParams = cell(1,length(inds));
for ii = 1:length(inds)
    dGCParams(ii) = {vals{inds(ii)}};
end

inds = [94 98 100 103 105 107];
pGCParams = cell(1,length(inds));
for ii = 1:length(inds)
    pGCParams(ii) = {vals{inds(ii)}};
end


% create new tables
% sqlquery = ['CREATE TABLE bo.inputparams(tsim varchar,Iext varchar,Inoise varchar,'...
%     'RespON varchar,RespFQ varchar,randfrac varchar)'];
% exec(conn,sqlquery);

test = fetch(conn,'select * from bo.inputparams;');


% insert data into table
tablename = 'bo.inputparams';
insert(conn,tablename,InputParams_cols,InputParams);



%% Save all parameters in DB in one table


% % create the table if it doesn't yet exist
% NOTE: This should only be run ONCE to initialize the table!
%
% sqlquery = ['CREATE TABLE bo.allparams(tsim varchar,Iext varchar,Inoise varchar,'...
%     'RespON varchar,RespFQ varchar,randfrac varchar,NMC varchar,NdGC varchar,NpGC varchar,'...
%     'NoiseMC varchar,NoisedGC varchar,NoisepGC varchar,pGCLLD varchar,CMG varchar,'...
%     'NP varchar,NFB varchar,NoiseP varchar,NoiseFB varchar,CPP varchar,CPF varchar,'...
%     'CFP varchar,flagLOT varchar,LOTvar varchar,flagOP varchar,OPvar varchar,CMP varchar,'...
%     'CPG varchar,MC_Vrest varchar,MC_Fthresh varchar,MC_WGM varchar,dGC_tNMDA1 varchar,dGC_tNMDA2 varchar,'...
%     'dGC_gmax varchar,dGC_Vrest varchar,dGC_Fthresh varchar,dGC_WMG varchar,pGC_tLLD varchar,pGC_ELLD varchar,'...
%     'pGC_gmax varchar,pGC_Vrest varchar,pGC_Fthresh varchar,pGC_WPG varchar,Date varchar,Comments varchar)'];
% 
% exec(conn,sqlquery);

% new params: tau1Ca tau2Ca VCAthresh

AllParams_cols = {'tsim','Iext','Inoise','RespON','RespFQ','randfrac',...
    'NMC','NdGC','NpGC','NoiseMC','NoisedGC','NoisepGC','pGCLLD','CMG',...
    'NP','NFB','NoiseP','NoiseFB','CPP','CPF','CFP','flagLOT','LOTvar',...
    'flagOP','OPvar','CMP','CPG','MC_Vrest','MC_Fthresh','MC_WGM','dGC_tNMDA1',...
    'dGC_tNMDA2','dGC_tCA1','dGC_tCA2','dGC_gmax','dGC_Vrest','dGC_Fthresh','dGC_VCAthresh','dGC_WMG','pGC_tLLD',...
    'pGC_ELLD','pGC_gmax','pGC_Vrest','pGC_Fthresh','pGC_WPG','Date','Comments'};


input_file = 'noisy_OB_PC_params_2CG.txt';
[names, vals] = textread(input_file,'%s %s',162);

inds = [3 6 7 8 9 10 14 15 16 19 20 21 24 30 32 33 36 37 39 40 41 43 44 45 ...
    46 49 50 67 69 71 77 78 79 80 83 84 86 87 89 97 101 103 106 108 110];
% inds of all relevant parameters

AllParams = cell(1,length(inds));
for ii = 1:length(inds)
    AllParams(ii) = {vals{inds(ii)}};
    % for some wierd reason I had to put {} on the rhs of = so that they
    % would show up as comumns in DB
end
dstr = datestr(clock);
AllParams(ii+1) = {dstr};
comm = ['testing stuff'];
AllParams(ii+2) = {comm};





% insert data into table
tablename = 'bo.allparams';
insert(conn,tablename,AllParams_cols,AllParams);


test = fetch(conn,'select * from bo.allparams;');






% delete tables
sqlquery = ['DROP TABLE bo.allparams'];
exec(conn,sqlquery);

