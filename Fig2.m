clear all;close all;clc

%% Exp-1
N = 100;
M = 200;

Rang = 20;

PT = zeros(Rang,Rang);

for K=1:1:Rang
    for J = 1:1:Rang
        fprintf('K = %d, J = %d\n',K,J);
        for Iter = 1:40
            % Dictionary A
            
%             % Gaussian
%             A = randn(N,M);
            
            % Fourier
            A = zeros(N,M);
            for k1 = 1:N
                Atemp = dftmtx(M);
                A(k1,:) = Atemp(randperm(M,1),:);
            end
            
            % Subspace matrix B
            B = dftmtx(N)/sqrt(N);
            B = B(:,1:K);
            % Linear sensing process Phi
            Phi = zeros(N,K*M);
            for j=1:M
                for i=1:K
                    Phi(:,K*(j-1)+i) = diag(B(:,i))*A(:,j);
                end
            end
            % Ground truth
            X0 = zeros(K,M);
            T = randperm(M,J);
            for k = 1:length(T)
                X0(:,T(k)) = randn(1)*randn(K,1);
            end
            % Observed signal
            y = Phi*reshape(X0,[K*M,1]);
            
            cvx_begin quiet
            variable X(K,M);
            minimize( sum( norms( X, 2, 1 ) ) )
            subject to
            y == Phi*reshape(X,[K*M,1]);
            cvx_end
            
            % Count the number of success recovery
            if norm(X-X0,'fro')/norm(X0,'fro') <= 1e-05
                PT(K,J) = PT(K,J) + 1;
            end
            
        end
        
    end
end

figure(1)
imagesc(PT/Iter);
h=colorbar;
ylabel(h, 'Recovery rate')
set(gca,'Ydir','normal')
xlabel('J from 1 to 20') 
ylabel('K from 1 to 20') 
xticks([2 4 6 8 10 12 14 16 18 20])
save('KJFourier.mat','PT')