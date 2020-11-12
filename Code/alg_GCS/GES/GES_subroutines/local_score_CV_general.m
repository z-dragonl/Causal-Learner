
function [score] = local_score_CV_general(Data,Xi,PAi,parameters)
% calculate the local score
% using negative k-fold cross-validated log likelihood as the score
% based on a regression model in RKHS

T=size(Data,1);
X = Data(:,Xi);
lambda = parameters.lambda; % regularization parameter
k = parameters.kfold; % k-fold cross validation
n0 = floor(T/k);
gamma = 0.01;
Thresh = 1E-5;

if(~isempty(PAi))
    PA = Data(:,PAi);
    
    % set the kernel for X
    GX = sum((X.*X),2);
    Q = repmat(GX,1,T);
    R = repmat(GX',T,1);
    dists = Q + R - 2*X*X';
    dists = dists-tril(dists);
    dists=reshape(dists,T^2,1);
    width = sqrt(0.5*median(dists(dists>0))); % median value
    width = width*2; %%%
    theta = 1/(width^2);
    
    Kx = kernel([X], [X], [theta,1]); % Gaussian kernel
    
    H0 =  eye(T) - ones(T,T)/(T); % for centering of the data in feature space
    Kx = H0 * Kx * H0; % kernel matrix for X
    
    [eig_Kx, eix] = eigdec((Kx+Kx')/2, min(400, floor(T/2))); % /2
    IIx = find(eig_Kx > max(eig_Kx) * Thresh); eig_Kx = eig_Kx(IIx); eix = eix(:,IIx);
    mx = length(IIx);

    % set the kernel for PA
    Kpa = ones(T,T);
    for m = 1:size(PA,2)
        G = sum((PA(:,m).*PA(:,m)),2);
        Q = repmat(G,1,T);
        R = repmat(G',T,1);
        dists = Q + R - 2*PA(:,m)*PA(:,m)';
        dists = dists-tril(dists);
        dists=reshape(dists,T^2,1);
        width = sqrt(0.5*median(dists(dists>0))); % median value
        width = width*2; %%%
        theta = 1/(width^2);
        Kpa = Kpa.*kernel([PA(:,m)], [PA(:,m)], [theta,1]);
    end
    H0 =  eye(T) - ones(T,T)/T; % for centering of the data in feature space
    Kpa = H0 * Kpa * H0; % kernel matrix for PA
    
    
    CV = 0;
    for kk=1:k
        if(kk==1)
            Kx_te = Kx((kk-1)*n0+1:kk*n0,(kk-1)*n0+1:kk*n0);
            Kx_tr = Kx(kk*n0+1:T,kk*n0+1:T);
            Kx_tr_te = Kx(kk*n0+1:T,(kk-1)*n0+1:kk*n0);
            Kpa_te = Kpa((kk-1)*n0+1:kk*n0,(kk-1)*n0+1:kk*n0);
            Kpa_tr = Kpa(kk*n0+1:T,kk*n0+1:T);
            Kpa_tr_te = Kpa(kk*n0+1:T,(kk-1)*n0+1:kk*n0);
            nv = n0; % sample size of validated data
        end
        if(kk==k)
            Kx_te = Kx((kk-1)*n0+1:T,(kk-1)*n0+1:T);
            Kx_tr = Kx(1:(kk-1)*n0,1:(kk-1)*n0);
            Kx_tr_te = Kx(1:(kk-1)*n0,(kk-1)*n0+1:T);
            Kpa_te = Kpa((kk-1)*n0+1:T,(kk-1)*n0+1:T);
            Kpa_tr = Kpa(1:(kk-1)*n0,1:(kk-1)*n0);
            Kpa_tr_te = Kpa(1:(kk-1)*n0,(kk-1)*n0+1:T);
            nv = T-(kk-1)*n0;
        end
        if(kk<k & kk>1)
            Kx_te = Kx((kk-1)*n0+1:kk*n0,(kk-1)*n0+1:kk*n0);
            Kx_tr = Kx([1:(kk-1)*n0,kk*n0+1:T],[1:(kk-1)*n0,kk*n0+1:T]);
            Kx_tr_te = Kx([1:(kk-1)*n0,kk*n0+1:T],(kk-1)*n0+1:kk*n0);
            Kpa_te = Kpa((kk-1)*n0+1:kk*n0,(kk-1)*n0+1:kk*n0);
            Kpa_tr = Kpa([1:(kk-1)*n0,kk*n0+1:T],[1:(kk-1)*n0,kk*n0+1:T]);
            Kpa_tr_te = Kpa([1:(kk-1)*n0,kk*n0+1:T],(kk-1)*n0+1:kk*n0);
            nv = n0;
        end
        n1 = T-nv;
        tmp1 = pdinv(Kpa_tr + n1*lambda*eye(n1));
        tmp2 = tmp1*Kx_tr*tmp1;
        tmp3 = tmp1*pdinv(eye(n1) + n1*lambda^2/gamma*tmp2)*tmp1;
        
        A = (Kx_te + Kpa_tr_te'*tmp2*Kpa_tr_te - 2*Kx_tr_te'*tmp1*Kpa_tr_te ...
            - n1*lambda^2/gamma*Kx_tr_te'*tmp3*Kx_tr_te ...
            - n1*lambda^2/gamma*Kpa_tr_te'*tmp1*Kx_tr*tmp3*Kx_tr*tmp1*Kpa_tr_te ...
            + 2*n1*lambda^2/gamma*Kx_tr_te'*tmp3*Kx_tr*tmp1*Kpa_tr_te)/gamma;
        
        B = n1*lambda^2/gamma * tmp2 + eye(n1);
        L = chol(B)';
        C = sum(log(diag(L)));
%         CV = CV + (nv*nv*log(2*pi) + nv*C + nv*mx*log(gamma) + trace(A))/2;
                CV = CV + (nv*nv*log(2*pi) + nv*C + trace(A))/2;

    end
    CV =CV/k;
else
    % set the kernel for X
    GX = sum((X.*X),2);
    Q = repmat(GX,1,T);
    R = repmat(GX',T,1);
    dists = Q + R - 2*X*X';
    dists = dists-tril(dists);
    dists=reshape(dists,T^2,1);
    width = sqrt(0.5*median(dists(dists>0))); % median value
    width = width*2; %%%
    theta = 1/(width^2);
    
    Kx = kernel([X], [X], [theta,1]); % Gaussian kernel
    
    H0 =  eye(T) - ones(T,T)/(T); % for centering of the data in feature space
    Kx = H0 * Kx * H0; % kernel matrix for X
    
    %[eig_Kx, eix] = eigdec((Kx+Kx')/2, min(400, floor(T/2))); % /2
    [eig_Kx, eix] = eigdec((Kx+Kx'+1e-5*eye(size(Kx,1)))/2, min(400, floor(T/2)));
    
    
    IIx = find(eig_Kx > max(eig_Kx) * Thresh);
    mx = length(IIx);

    CV = 0;
    for kk=1:k
        if(kk==1)
            Kx_te = Kx((kk-1)*n0+1:kk*n0,(kk-1)*n0+1:kk*n0);
            Kx_tr = Kx(kk*n0+1:T,kk*n0+1:T);
            Kx_tr_te = Kx(kk*n0+1:T,(kk-1)*n0+1:kk*n0);
            nv = n0;
        end
        if(kk==k)
            Kx_te = Kx((kk-1)*n0+1:T,(kk-1)*n0+1:T);
            Kx_tr = Kx(1:(kk-1)*n0,1:(kk-1)*n0);
            Kx_tr_te = Kx(1:(kk-1)*n0,(kk-1)*n0+1:T);
            nv = T-(kk-1)*n0;
        end
        if(kk<k & kk>1)
            Kx_te = Kx((kk-1)*n0+1:kk*n0,(kk-1)*n0+1:kk*n0);
            Kx_tr = Kx([1:(kk-1)*n0,kk*n0+1:T],[1:(kk-1)*n0,kk*n0+1:T]);
            Kx_tr_te = Kx([1:(kk-1)*n0,kk*n0+1:T],(kk-1)*n0+1:kk*n0);
            nv = n0;
        end
        n1 = T-nv;
        A = (Kx_te - 1/(gamma*n1)*Kx_tr_te'*pdinv(eye(n1)+1/(gamma*n1)*Kx_tr)*Kx_tr_te)/gamma;
        
        B = 1/(gamma*n1)*Kx_tr + eye(n1);
        L = chol(B)';
        C = sum(log(diag(L)));
        
%         CV = CV + (nv*nv*log(2*pi) + nv*C + nv*mx*log(gamma) + trace(A))/2;
        CV = CV + (nv*nv*log(2*pi) + nv*C + trace(A))/2;

    end
    CV =CV/k;
end

score = CV; % negative cross-validated likelihood












