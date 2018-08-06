function bool = isOctave()
% ISOCTAVE determines whether this function was called in Octave.
%
% authors:
% Thorben Casper, David Duque, Victoria Heinz, Abdul Moiz,
% Herbert De Gersem, Sebastian Schoeps
% Institut fuer Theorie Elektromagnetischer Felder
% Graduate School of Computational Engineering
% Technische Universitaet Darmstadt

bool = exist('OCTAVE_VERSION', 'builtin') ~= 0;

end