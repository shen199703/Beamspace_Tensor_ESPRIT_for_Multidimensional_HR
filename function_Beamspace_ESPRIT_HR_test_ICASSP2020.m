% test function_Beamspace_ESPRIT_HR
% 2020/04/30

clear
clc

tic

I1 = 8; % 1st dim (sensors)

I2 = 8; % 2nd dim (sensors)

I3 = 8; % 3rd dim (sensors)

% % ----------- Fig 2 ------------
% f1 = [0.55 0.719 0.906];  % 1st dim [-pi,pi]
% f2 = [0.41 0.777 0.276];  % 2nd dim [-pi,pi]  
% f3 = [0.34 0.906 0.358];  % 3rd dim [-pi,pi]

% % ----------- Fig 3 ------------
f1 = [0.2 0.2 0.6 0.8]; % 1st dim [-pi,pi]
f2 = [0.9 0.4 0.4 0.6];
f3 = [0.1 0.2 0.8 0.8];

% % ----------- Fig 4 ------------
% f1 = [0.1 0.3 0.8]; % 1st dim [-pi,pi]
% f2 = [0.9 0.4 0.2];
% f3 = [0.4 0.1 0.7];

R = length(f1); % source number is known

A1 = zeros(I1,R);
A2 = zeros(I2,R);
A3 = zeros(I3,R);

for r = 1:R
    z1 = exp(f1(r)*1i*pi);
    for i = 1:I1 % 1st dim
        A1(i,r) = z1^(i-1);
    end
    z2 = exp(f2(r)*1i*pi);
    for i = 1:I2 % 2nd dim
        A2(i,r) = z2^(i-1);
    end
    z3 = exp(f3(r)*1i*pi);
    for i = 1:I3 % 3rd dim
        A3(i,r) = z3^(i-1);
    end
end

%------------------------------------------------------
type = 'structured'; % r: random, s: structured | beams
%------------------------------------------------------
mum_sam = 10; % number of measurements
S = ones(mum_sam,R) + 1*randn(mum_sam,R);
 
% attenuation_vector = [1 0.85 0.75 0.5]; % attenuation for different sources
attenuation_vector = rand(1,R); % random attenuation
for r = 1:R % source index
    S(:,r) = attenuation_vector(r)*S(:,r);
end

beams = 6; % number of beams
J1 = beams; % 1st dim (beamspace)
J2 = beams; % 2nd dim (beamspace)
J3 = beams; % 3rd dim (beamspace)

runs = 20; % runs

f1_est = zeros(runs,R);
f2_est = zeros(runs,R);
f3_est = zeros(runs,R);

for it = 1:runs
    disp(it)
    switch lower(type)
        case 'random' % random beams 
            B1 = randn(I1,J1)+1j*randn(I1,J1); % random is simple
            B2 = randn(I2,J2)+1j*randn(I2,J2); %
            B3 = randn(I3,J3)+1j*randn(I3,J3); %
            
        case 'structured' % structured beams
            B1 = exp(1j*(0:I1-1).'*(1:J1)/(J1+1)*pi); % uniform grids
            B2 = exp(1j*(0:I2-1).'*(1:J2)/(J2+1)*pi);
            B3 = exp(1j*(0:I3-1).'*(1:J3)/(J3+1)*pi); 

        otherwise % error
            disp('accepted types: random or structured')
    end
    
    A1b = B1'*A1; % manifold in beamspace: 1st
    A2b = B2'*A2; % manifold in beamspace: 2nd
    A3b = B3'*A3; % manifold in beamspace: 3rd
    
    T = cpdgen({A1b,A2b,A3b,S}); % tensor without noise
    SNR = 20; % SNR(dB)
    Tnoisy = noisy(T,SNR); % noisy tensor measurments
 
    Uest = cpd(Tnoisy,R); % CP decomposition (NLS)
    
    A1fix = Uest{1}; % subspace: 1st dim
    A2fix = Uest{2}; % subspace: 2nd dim
    A3fix = Uest{3}; % subspace: 3rd dim
 
    for r = 1:R % one by one
        f1_est(it,r) = function_Beamspace_ESPRIT_HR(A1fix(:,r),B1,J1,1); % 1st dim
        f2_est(it,r) = function_Beamspace_ESPRIT_HR(A2fix(:,r),B2,J2,1); % 2nd dim
        f3_est(it,r) = function_Beamspace_ESPRIT_HR(A3fix(:,r),B3,J3,1); % 3rd dim
    end
    
end
       
figure % 3d plot
plot3(f1_est(:),f2_est(:),f3_est(:),'.'),grid
hold on,
plot3(f1,f2,f3,'rs','MarkerSize',10,'LineWidth',1.5) 
hold off
axis([0 1 0 1 0 1])
legend('Estimated Frequency','Actual Frequency')
xlabel('1st dim (\pi)','fontsize',12)
ylabel('2nd dim (\pi)','fontsize',12)
zlabel('3rd dim (\pi)','fontsize',12) 
 
for itr = 1:R
    text(f1(itr)+0.1,f2(itr)+0.1,f3(itr),['Source:',num2str(itr)],'HorizontalAlignment','left','FontSize',10);
end
 
toc