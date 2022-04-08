function out = array_to_rowLW(in,pixel_on)
%array_to_rowLW converts vector input and active pixel coordinates into a
% row matrix containing three-dimensional values.

% Initialisation
[ROWS,COLS] = getConstantsLW;
npon = length(pixel_on);
out = zeros(npon,3);

% Matrix generalisation
for p = 1:npon
    rc = v2a(pixel_on(p),ROWS,COLS);
    out(p,:) = [rc(1) rc(2) in(p)];
end