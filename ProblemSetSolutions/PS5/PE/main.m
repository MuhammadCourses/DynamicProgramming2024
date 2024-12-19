
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file helps make any changes to parameter choices for partial equilibrium.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

p = define_parameters();
% loop over two values of r and call main_PE.m
r = [0.035, 0.04];
for i = 1:2
    p.r = r(i);
    main_PE
end