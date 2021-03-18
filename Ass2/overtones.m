function [overtones] = overtones(notes,levels)
%overtones Returns a matrix of the overtones of the given vector of notes
%up to the levels-th overtone. Matrix returned is of size n by levels, where notes has
%length n.
A = transpose(16.35*(2^(1/12)).^notes);
B = 1:levels;
C = A*B;
overtones = log(C/16.35) / log(2^(1/12));
end

