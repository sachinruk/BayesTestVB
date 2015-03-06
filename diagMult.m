function res=diagMult(A,d,side)
%returns the result of right/left multiplying by a diagonal
%inputs:
%A: non diagonal matrix
%d: vector of diagonal elements
%side: 'l' or 'r'
if strcmp(side,'l') %left multiply the DIAGONAL
    res=bsxfun(@times,d,A);
else %right multiply the DIAGONAL
    res=bsxfun(@times,A,d');
end
    
    
    