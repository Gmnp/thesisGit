%% Setting up defaults.

close all hidden
clear all

randn('state',50);
rand('state',50);

height = 2/3;
set(0,...
'defaultfigureposition',[180 100 800 800*height],...
'defaultaxeslinewidth',1,...
'defaultaxesfontsize',11,...%16
'defaultlinelinewidth',2,... %2
'defaultpatchlinewidth',1,... %1
'defaultlinemarkersize',6,... %10
'defaulttextinterpreter','tex');
