
function [score] = local_score_marginal_general(Data, Xi,PAi,parameters)
% calculate the local score by negative marginal likelihood
% based on a regression model in RKHS

T=size(Data,1);
X = Data(:,Xi);
dX = size(X,2);

% set the kernel for X
GX = sum((X.*X),2);
Q = repmat(GX,1,T);
R = repmat(GX',T,1);
dists = Q + R - 2*(X*X');
dists = dists-tril(dists);
dists=reshape(dists,T^2,1);
widthX = sqrt(0.5*median(dists(dists>0)));
widthX = widthX*2.5; % kernel width
theta = 1/(widthX^2);
H =  eye(T) - ones(T,T)/T;
Kx = kernel(X, X, [theta,1]);
Kx = H * Kx * H;

Thresh = 1E-5;
[eig_Kx, eix] = eigdec((Kx+Kx')/2, min(400, floor(T/4))); % /2
IIx = find(eig_Kx > max(eig_Kx) * Thresh);
eig_Kx = eig_Kx(IIx);
eix = eix(:,IIx);

if(~isempty(PAi))
    PA = Data(:,PAi);
    
    % set the kernel for PA
    for m = 1:size(PA,2)
        G = sum((PA(:,m).*PA(:,m)),2);
        Q = repmat(G,1,T);
        R = repmat(G',T,1);
        dists = Q + R - 2*PA(:,m)*PA(:,m)';
        dists = dists-tril(dists);
        dists=reshape(dists,T^2,1);
        widthPA(m) = sqrt(0.5*median(dists(dists>0)));
    end
    widthPA = widthPA*2.5; % kernel width
    
    covfunc = {'covSum', {'covSEard','covNoise'}};
    logtheta0 = [log(widthPA'); 0; log(sqrt(0.1))];
    [logtheta, fvals, iter] = minimize(logtheta0, 'gpr_multi_new', -300, covfunc, PA, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    
    [nlml dnlml] = gpr_multi_new(logtheta, covfunc, PA, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    
else
    
    covfunc = {'covSum', {'covSEard','covNoise'}};
    PA = zeros(T,1);
    logtheta0 = [100; 0; log(sqrt(0.1))];
    [logtheta, fvals, iter] = minimize(logtheta0, 'gpr_multi_new', -300, covfunc, PA, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    [nlml dnlml] = gpr_multi_new(logtheta, covfunc, PA, 2*sqrt(T) *eix * diag(sqrt(eig_Kx))/sqrt(eig_Kx(1)));
    
end

score = nlml; % negative log-likelihood











