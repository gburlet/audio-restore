% Gregory Burlet
% MUMT 605
% December 8, 2011
%
% This script is a test bed for the sinusoid + AR residual click removal
% algorithm.

clear all
clc

% INIT vars
noiseyAudio = 'mussorsky';
N = 2048;        
hopSize = 1024;   
p = 31;                   
q = 31;
detThresh = 4;                      
detStretch = 4;
numIter = 5;

[x, fs, nbits] = wavread([noiseyAudio, '.wav']);  

% DEBUG
% shorten y for testing: mono for x secs
%secs = 3;
%x = x(1:fs*secs,1);

% only mono
%x = x(:,1);

tic;
[x_clean, ~] = sinARdeclick(x, p, q, N, hopSize, detThresh, detStretch, numIter);
tElapsed = toc

display('done');
display('writing to file');

wavwrite(x_clean, fs, nbits, [noiseyAudio, '_clean.wav']);