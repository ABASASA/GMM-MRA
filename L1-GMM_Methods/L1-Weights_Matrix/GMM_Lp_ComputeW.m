function [W] = GMM_Lp_ComputeW(invarseCovariance, D)
% using the notation of the article "L_p properties..."
W2 = null(D');
W1 = (invarseCovariance) * D;

W = [W1, W2]';

end