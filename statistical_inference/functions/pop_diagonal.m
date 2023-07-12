function diagonal = pop_diagonal(TGM, varargin)
%POP_DIAGONAL pulls out the diagonal from a Temporal Generalisation Matrix
%  
%   Input:
%   TGM         time-point by time-point square Temporal Generalisation
%               Matrix
%
%   Output:
%   diagonal    1 by time-point time-series corresponding to the input's
%               diagonal

    [x, y] = size(TGM);
    assert(x==y, 'Inputted TGM is not square. Ensure that rows and columns are the same size.')
    diagonal = squeeze(max(eye(x).*TGM));
end