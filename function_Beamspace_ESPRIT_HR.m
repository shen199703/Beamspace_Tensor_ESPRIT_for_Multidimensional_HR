% function: matlab code for beamspace ESPRIT algorithm by Fuxi, 2020/04/30

% Input
% Ubs:     signal subspace estimated from beamspace measurements
% B:       precoder or combiner matrix
% numB:    number of beams
% numS:	   number of sources

% Output
% estimated_frequency: estimated HR frequency

function estimated_frequency = function_Beamspace_ESPRIT_HR(Ubs,B,numB,numS)

B1 = B(1:end-1,:);
B2 = B(2:end,:);

F = pinv(B2)*B1;

Bh = B';
tt = [Bh(:,end) F'*Bh(:,1)];
[Qt,~] = qr(tt,0);

Q = eye(numB,numB) - Qt*Qt'; 

E1 = Q*Ubs;
E2 = Q*F'*Ubs;

EE = [E1';E2']*[E1 E2];

[E,~] = eigs(EE,2*numS);

E12 = E(1:numS,numS+1:end);
E22 = E(numS+1:end,numS+1:end);

[~,EigVal] = eigs(-E12*pinv(E22));  % eigen value in matrix form

Eval = diag(EigVal).';  % eigen values (row vector)

estimated_frequency = wrapToPi(angle(Eval))/pi; % angle of the eigen values | Wrap angle in radians to [?pi pi]