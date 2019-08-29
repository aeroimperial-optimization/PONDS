function [] = UODESys(n,g,input_file,output_file,solver,verbose)
% Implementation of UODESys package tool developed by Mayur Venkatram Lakshmi
% Reduces a finite-dimensional system into a smaller uncertain system
%
% n             Size of the reduced uncertain system
% g             [g1 g2 g3] to be used by UODESys package tool
% input_file    name of the .mat file for data input to UODESys
% output_file   name of the .mat output file
%
% Written by Mario Lino Valencia (June 2019)
% Imperial College London - Department of Aeronautics

Transformation(input_file);
Bounds_Calculation(n,g,output_file,solver,verbose);