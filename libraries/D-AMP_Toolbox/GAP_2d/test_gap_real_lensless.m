
clear all
close all
clc

infile = 'pic_sim64k_books';
% for seperate reconstruction file
% infile = 'pic_meas64dualwh_1107'

vec_len = round(65536*0.125);
load(infile);
tmp1 = Yp(1:vec_len,1:3)-Yn(1:vec_len,1:3);
tmp2 = Yp(1:vec_len,4:6)-Yn(1:vec_len,4:6);
N = 2^16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CS RECONSTRUCTION %%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Making measurement matrix Phi
load perm64k;

