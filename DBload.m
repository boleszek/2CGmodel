
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% before running DBload.m you must ssh into the imb-db server by typing the
% following 
%
% ssh boleszek@cronusx.uchicago.edu -L 5433:cronusx:5433
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

clc
driver='org.postgresql.Driver';
% url='jdbc:postgresql://localhost/';
url='jdbc:postgresql://localhost:5433/imb-db';

username= 'boleszek';
password= '8TFO0z4h';

conn=database('imb-db',username,password,driver,url)
%               'driver',driver,'URL',url)
isconnection(conn)
assert(strcmp(get(conn,'AutoCommit'),'on'))
setdbprefs('DataReturnFormat','cellarray');
dbp=setdbprefs;
assert(strcmp(dbp.DataReturnFormat,'cellarray'))
