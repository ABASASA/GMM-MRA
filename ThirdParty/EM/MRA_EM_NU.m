function[x, rho] = MRA_EM_NU(X, sigma, x,rho, tol, batch_niter)
% A fast mplementation of the afapted NonUniform EM.
% Using partial set of samples when N>5000
%
% Input:
% x -- the initial guess, if not appears then it is taken as random
% tol -- for alting criterion
% batch_niter -- number of iterations with a subset of samples
%

% X contains N observations, each of length L
[L, N] = size(X);

% Initial guess of the signal
if ~exist('x', 'var') || isempty(x)
    if isreal(X)
        x = randn(L, 1);
    else
        x = randn(L, 1) + 1i*randn(L, 1);
    end
end
x = x(:);
assert(length(x) == L, 'Initial guess x must have length L.');

% Tolerance to declare convergence
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
end


% In practice, we iterate with the DFT of the signal x
fftx = fft(x);
if ~exist('rho', 'var') || isempty(rho)

    rho  = rand(L,1); rho = rho/sum(rho);
end

% Precomputations on the observations
fftX = fft(X);
sqnormX = repmat(sum(abs(X).^2, 1), L, 1);

% If the number of observations is large, get started with iterations
% over only a sample of the observations
if N >= 5000
    
    if ~exist('batch_niter', 'var') || isempty(batch_niter)
        batch_niter = 3000;
    end
    batch_size = 1000;
    
    for iter = 1 : batch_niter
       % sample = randi(N, batch_size, 1);
        [fftx_new, rho_new] = EM_iteration(fftx, rho, fftX, sqnormX, sigma);
        fftx = fftx_new;
        rho = rho_new;
    end
end

% In any case, finish with full passes on the data
full_niter = 10000;
for iter = 1 : full_niter
    [fftx_new, rho_new] = EM_iteration(fftx, rho, fftX, sqnormX, sigma);
    if relative_error(ifft(fftx), ifft(fftx_new)) < tol
        break;
    end
    fftx = fftx_new;
    rho = rho_new;
end
fprintf('\t\tEM UN: %d full iterations\n', iter);
x = ifft(fftx);
end


% Execute one iteration of EM with current estimate of the DFT of the
% signal given by fftx, and DFT's of the observations stored in fftX, and
% squared 2-norms of the observations stored in sqnormX, and noise level
% sigma.
function[fftx_new, rho_new] = EM_iteration(fftx, rho, fftX, sqnormX, sigma)
% The magic is done here...
C = ifft(bsxfun(@times, conj(fftx), fftX));
% the distance
T = (2*C - sqnormX)/(2*sigma^2);
T = bsxfun(@minus, T, max(T, [], 1));
W = exp(T);
% update with rho as a weight Expectation
W = bsxfun(@times, rho, W);
W = bsxfun(@times, W, 1./sum(W, 1));
% thew new estimations Maximization
fftx_new = mean(conj(fft(W)).*fftX, 2);
rho_new = mean(W,2);
end





