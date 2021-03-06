function writeparams(input_file,mode,comment,conn,fname)

% This function allows you to write a subset of parameters from the 
% input_file to the DB or to a text file, depending on the mode
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS
%
% input_file   -    name of parameter txt file (noisy_OB_PC_params_2CG.txt)
% mode         -    0: insert parameters into db table
%                   1: write parameters as .dat file
% comment      -    add a comment to append to the parameter table if desired
% conn         -    DB connection address (only used if mode = 0)
% fname        -    if mode = 0 fname is table name (bo.allparams)
%              -    if mode = 1 fname is file name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % create the table if it doesn't yet exist
% NOTE: This should only be run ONCE to initialize the table!
%
% sqlquery = ['CREATE TABLE bo.allparams(tsim varchar,Iext varchar,Inoise varchar,'...
%     'RespON varchar,RespFQ varchar,randfrac varchar,NMC varchar,NdGC varchar,NpGC varchar,'...
%     'NoiseMC varchar,NoisedGC varchar,NoisepGC varchar,pGCLLD varchar,CMG varchar,'...
%     'NP varchar,NFB varchar,NoiseP varchar,NoiseFB varchar,CPP varchar,CPF varchar,'...
%     'CFP varchar,flagLOT varchar,LOTvar varchar,flagOP varchar,OPvar varchar,CMP varchar,'...
%     'CPG varchar,MC_Vrest varchar,MC_Fthresh varchar,MC_WGM varchar,dGC_tNMDA1 varchar,dGC_tNMDA2 varchar,dGC_tCA1 varchar,dGC_tCA2 varchar,'...
%     'dGC_gmax varchar,dGC_Vrest varchar,dGC_Fthresh varchar,dGC_VCAthresh varchar,dGC_WMG varchar,pGC_tLLD varchar,pGC_ELLD varchar,'...
%     'pGC_gmax varchar,pGC_Vrest varchar,pGC_Fthresh varchar,pGC_WPG varchar,Date varchar,Comments varchar)'];
% 
% exec(conn,sqlquery);


AllParams_cols = {'tsim','Iext','Inoise','RespON','RespFQ','randfrac',...
    'NMC','NdGC','NpGC','NoiseMC','NoisedGC','NoisepGC','pGCLLD','CMG',...
    'NP','NFB','NoiseP','NoiseFB','CPP','CPF','CFP','flagLOT','LOTvar',...
    'flagOP','OPvar','CMP','CPG','MC_Vrest','MC_Fthresh','MC_WGM','dGC_tNMDA1',...
    'dGC_tNMDA2','dGC_tCA1','dGC_tCA2','dGC_gmax','dGC_Vrest','dGC_Fthresh','dGC_VCAthresh','dGC_WMG','pGC_tLLD',...
    'pGC_ELLD','pGC_gmax','pGC_Vrest','pGC_Fthresh','pGC_WPG','Date','Comments'};



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




if mode == 0
    % insert data into table
    
    insert(conn,name,AllParams_cols,AllParams);
    
elseif mode == 1
    % write data to txt file
    
    names = [names;'Comment'];
    vals = [vals;num2str(comment)];
    nrows = length(names);
    ca = [names, vals]; % create cell array
    filename = ['params/' fname '.dat'];
    fileID = fopen(filename,'w'); % create new file for writing 'w'
    formatSpec = '%s %s\n';
    
    for row = 1:nrows
        fprintf(fileID,formatSpec,ca{row,:});
    end
    
    fclose(fileID);
end
    
    


