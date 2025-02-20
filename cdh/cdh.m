clear;clc;close all

% Mars Hand Lens Imager
num_pixels = 1280*720;
fps = 2;
bits_per_pixel = 24;
mhli_kbps = fps*num_pixels*bits_per_pixel / 1000;

% REMS
% Each hour, every sol, REMS will record 5 minutes of data at 1 Hz for all sensors
sol = 24*60*60 + 39*60; % sol to seconds
total_data_size_per_sol = (14.73 *8e6)*1e-3; % bits/sol
rems_kbps = total_data_size_per_sol/sol; % kbps
% rate = 1; % Hz

% HiRise
uncompressed_image_size = 28* % kbit