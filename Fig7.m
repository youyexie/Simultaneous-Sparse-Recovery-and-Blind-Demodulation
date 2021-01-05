clear all;close all;clc

%% Exp-2
J = 5;
K = 5;
M = 200;
N = 100;
numIter = 40;

PT = zeros(numIter,9);

for eta = -3:5
    fprintf('eta = %d\n',eta);
    for Iter = 1:numIter
        % Dictionary A
        
        % Gaussian
        A = randn(N,M);
        
          % Subspace matrix B
        B = dftmtx(N)/sqrt(N);
        B = B(:,1:K);
        % Linear sensing process Phi
        Phi = zeros(N,K*M);
        for j=1:M
            for j2=1:K
                Phi(:,K*(j-1)+j2) = diag(B(:,j2))*A(:,j);
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
        
        % Noise
        x0norm = norm(X0,'fro');
        n_unjusted = randn([N,1])+randn([N,1])*1i;
        scale = (10^(-eta/2))*x0norm/norm(n_unjusted);
        n = scale.*n_unjusted;
        
        % Add noise to the observed signal
        y = y + n;
        
        cvx_begin quiet
        variable X(K,M);
        minimize( sum( norms( X, 2, 1 ) ) )
        subject to
        norms( Phi*reshape(X,[K*M,1]) - y,2 ) <= norm(n);
        cvx_end
        
        % Count the number of success recovery
        PT(Iter,eta+4) =  norm(X-X0,'fro')/norm(X0,'fro') ;
        
    end
    
end
 
%%
meanPT = mean(PT);
%logPT =  fliplr( 20*log10(PT));
logmeanPT = fliplr( 20*log10(meanPT) );
figure(4);
plot( logmeanPT ,'-' );
hold on

logstdPT = fliplr(20*log10(std(PT)+meanPT));

stdrange = abs(logstdPT - logmeanPT);

for k2 = 1: length( meanPT )
   line([k2,k2],[logmeanPT(k2)+stdrange(k2),logmeanPT(k2)-stdrange(k2)]);%[max(logPT(:,k2)),min(logPT(:,k2))]);
   line([k2-0.3,k2+0.3],[logmeanPT(k2)+stdrange(k2),logmeanPT(k2)+stdrange(k2)]);%[max(logPT(:,k2)),max(logPT(:,k2))])
   line([k2-0.3,k2+0.3],[logmeanPT(k2)-stdrange(k2),logmeanPT(k2)-stdrange(k2)]);%[min(logPT(:,k2)),min(logPT(:,k2))])
end

plot([1,2,3,4,5,6,7,8,9 ],[[-50,-40,-30,-20,-10,0,10,20,30]+36.3794],'--');
scatter(1:9,logmeanPT,'+','r','LineWidth',2);
hold off
grid on
%axis tight
xlabel('Noise to signal ratio (dB)')
xticks([1 2 3 4 5 6 7 8 9])
xticklabels({'-50','-40','-30','-20','-10','0','10','20','30'})
ylabel('Relative error (dB)')
save('20logNoisyGaussian.mat','PT')

%% 
% gc =  5*sqrt(6)+24*sqrt(5);             65.9131
% gamma = sqrt(2*200*log(2*5*200)+2*200+1)£»     58.6631
% P = log(4*sqrt(2*5)*gamma)/log(2)£»            9.5353
% fc = 5*sqrt(6)+24*sqrt(5*P)£»         177.9635
%
% 20*log10(gc)= 18.1897*2 = 36.3794
% 20*log10(fc)= 22.5033*2 = 45.0066