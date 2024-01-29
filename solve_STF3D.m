function u = solve_STF3D(m_Gamma, m_Omega, k0, k1, dir, a1, tol)

Omega = dom(m_Omega, 4);
Gamma = dom(m_Gamma, 3);

Uh = fem(m_Gamma, 'P1');
I = integral(Gamma, Uh, Uh);
one = ones(1, length(m_Gamma.vtx));

% kO = one * I * kx(m_Gamma.vtx) / sum(ndv(m_Gamma));
kO = k1/sqrt(a1);
% kO = one * I / sum(ndv(m_Gamma));
% kO = 6;
% beta = @(X) kx(X).^2 - kO^2;

%%
Uh = fem(m_Omega, 'P1');


A            = build_STF3D(m_Gamma, k0, kO,a0, ax, tol);
pots         = build_potentials3D(m_Gamma, m_Omega, kO, tol);
[F, Fbnd]    = build_rhs3D(m_Gamma, m_Omega, k0, dir); 


% A = A0+A1;
% traceAv = [Ad; An];
% traceAv = [An; Ad];


% pots = [D -S];
Id = integral(Omega, Uh, Uh);

traces = A \ Fbnd;

sol = - pots * traces;

sol = Id \ sol;
% 
% M = bmm({A0.blk{1,1} + A1.blk{1,1}, A0.blk{1,2} + A1.blk{1,2}, An; ...
%          A0.blk{2,1} + A1.blk{2,1}, A0.blk{2,2} + A1.blk{2,2}, Ad; ...
%          D, -S, Id-Av});
%           
% M = [(A0.blk{1,1}+A1.blk{1,1}) (A0.blk{1,2} + A1.blk{1,2}) An; ...
%          (A0.blk{2,1} + A1.blk{2,1}) (A0.blk{2,2}+A1.blk{2,2}) Ad; ...
%          D -S Id-Av];



% 
% M = [A traceA; ...
%      pots Id-Av];

% M = [A traceAv;pots Id-Av];
% M = bmm({A, traceAv; pots, Id-Av});

% M0 = bmm({A, sparse(0*traceAv);sparse(0*pots), Id});
% 
% A0 = [A0.blk{1,1} A0.blk{1,2};A0.blk{2,1} A0.blk{2,2}];
% A1 = [A1.blk{1,1} A1.blk{1,2};A1.blk{2,1} A1.blk{2,2}];
    

 
% M0 = bmm({A, sparse(size(Ad, 1)+size(An,1), size(Ad, 2)); ...
%           sparse(size(S, 1), size(S,2)+ size(D, 2)), Id});


% sol = M \ F;
% 
% phi     = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
% lambda  = sol(1:size(m_Gamma, 1)); sol(1:size(m_Gamma, 1)) = [];
u       = sol;


end