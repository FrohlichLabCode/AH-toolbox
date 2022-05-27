function [median_sem] = medianSem(mat)
% This function takes in a matrix of value (eg.power or PLV), normalize
% based on row mean and calculate median and sem of each column 

mat_normed = diag(1./mean(mat,2))*mat.*mean(mat(:));
median_sem = [nanmedian(mat_normed,1);
             nanstd(mat_normed,1)/sqrt(size(mat_normed,1))];

