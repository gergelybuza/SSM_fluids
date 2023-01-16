function vec = readd(W,vec)
% readds zeros in place of the removed lines
% (which is the correct solution)
bar = W.bar;
vec = [vec(1:W.N,:); zeros(2*W.N,size(vec,2)); vec(W.N+1:W.Q*bar-2*W.N,:);...
    zeros(bar,size(vec,2)); vec(W.Q*bar-2*W.N+1:end,:)];
end