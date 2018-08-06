% this script calls the different test cases to reproduce all numerical
% results from the following paper:
%
% T. Casper, U. Roemer, H. De Gersem, S. Schoeps. Coupled Simulation of
% Transient Heat Flow and Electric Currents in Thin Wires: Application to
% Bond Wires in Microelectronic Chip Packaging. Computers and Mathematics
% with Applications, submitted.
%
% authors:
% Thorben Casper, Ulrich Roemer, Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

clear;
close all;
addpath(genpath('src'));
verbose = [0 0];                                                           % triggers console outputs and plots
refinement = 'Fine';                                                       % {'Coarse','Fine'}

% if Octave is used, add compatOctave to path for compatibility
if isOctave, addpath('compatOctave'); end

% run individual test cases
resistor2D(refinement,verbose);
resistor3Dstraight(refinement,verbose);
resistor3Dbent(refinement,verbose);
chip(refinement,verbose);