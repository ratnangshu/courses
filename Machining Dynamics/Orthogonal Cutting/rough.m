clc;
clear;

%% Experimental Data
T = readtable('Book1.xlsx'); % table containing experimental data 
Feed_arr = table2array(T(:,1));
Ft_arr = table2array(T(:,2));
Ff_arr = table2array(T(:,3));
Chipt_arr = table2array(T(:,4));
figure(1);
plot(Feed_arr,Ft_arr,'b*'); hold on;
plot(Feed_arr,Ff_arr,'r*');

ft = fit(Feed_arr,Ft_arr,'exp1');
plot(ft,'b')
ft = fit(Feed_arr,Ff_arr,'exp1');
plot(ft,'r')