function Uout = adjustphi(U,Win,Wout)
% adjusts U according to changes in discretization
% U here is in complex coefficient form; Win is the previous W

Uout = zeros((Wout.dim-2)/2,1); % complex dim
bar = Wout.bar;
barin = Win.bar;
Q = min([Wout.Q Win.Q]);
omFx = U(end-1:end);
U = tocomp(U(1:end-2));
if Win.N < Wout.N
    for j = 1:Q
        for i = 1:Wout.nvar
        Uout((j-1)*bar+(i-1)*Wout.N+1:(j-1)*bar+(i-1)*Wout.N+Win.N) = ...
            U((j-1)*barin+(i-1)*Win.N+1:(j-1)*barin+i*Win.N);
        end
    end
else
    for j = 1:Q
        for i = 1:Wout.nvar
        Uout((j-1)*bar+(i-1)*Wout.N+1:(j-1)*bar+i*Wout.N) = ...
            U((j-1)*barin+(i-1)*Win.N+1:(j-1)*barin+(i-1)*Win.N+Wout.N);
        end
    end
end
Uout = toreal(Uout);
Uout = [Uout; omFx];
end