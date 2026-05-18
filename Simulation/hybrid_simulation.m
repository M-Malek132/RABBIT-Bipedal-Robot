clear; clc; close all;
% Get the directory of the current script
current_dir = fileparts(mfilename('fullpath'));

% Get the parent directory
parent_dir = fileparts(current_dir);

% Add parent directory and all its subfolders to MATLAB path
addpath(genpath(parent_dir));
