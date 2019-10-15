% THIS FILE HAS BEEN WRITTEN BY pwfitobject#tomatlab.m %
%
% This file is part of GTMpw -- Piecewise polynomial model of the GTM
% published under the GNU General Public License v3.

alpha0 = 2.8119e-01;

%% Cx.alpha(alpha,varargin)
Cx.alpha1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 3.8732e-02 + 2.4362e-01.*alpha + 4.4516e+00.*alpha.^2 - 1.7394e+01.*alpha.^3 ;
Cx.alpha2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 1.8840e-02 - 1.3042e-01.*alpha + 1.6872e-01.*alpha.^2 - 2.2361e-02.*alpha.^3 ;

%% Cz.alpha(alpha,varargin)
Cz.alpha1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 1.6739e-02 - 5.2414e+00.*alpha - 1.8661e+00.*alpha.^2 + 2.8466e+01.*alpha.^3 ;
Cz.alpha2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 3.6475e-01 - 2.7116e+00.*alpha + 1.6470e+00.*alpha.^2 - 3.6917e-01.*alpha.^3 ;

%% Cm.alpha(alpha,varargin)
Cm.alpha1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 1.1918e-01 - 1.4654e+00.*alpha + 8.1294e+00.*alpha.^2 - 3.1983e+01.*alpha.^3 ;
Cm.alpha2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 2.4671e-01 - 2.8471e+00.*alpha + 2.7475e+00.*alpha.^2 - 1.1045e+00.*alpha.^3 ;

%% Cx.beta(alpha,beta,varargin)
Cx.beta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 1.1672e-02 + 1.2703e-02.*alpha + 1.5769e-03.*beta - 2.0488e+00.*alpha.^2 + 2.7274e-02.*alpha.*beta + 6.6188e-02.*beta.^2 + 9.8081e+00.*alpha.^3 + 2.4895e-01.*alpha.^2.*beta - 5.7230e-01.*alpha.*beta.^2 - 5.8807e-03.*beta.^3 - 1.0651e+01.*alpha.^4 - 1.1777e+00.*alpha.^3.*beta + 1.9362e+00.*alpha.^2.*beta.^2 - 7.5085e-03.*alpha.*beta.^3 - 4.2589e-02.*beta.^4 ;
Cx.beta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 8.7327e-03 + 1.0051e-01.*alpha + 2.5540e-03.*beta - 2.3882e-01.*alpha.^2 + 2.5486e-03.*alpha.*beta + 2.5159e-02.*beta.^2 + 1.9870e-01.*alpha.^3 - 7.3374e-03.*alpha.^2.*beta + 1.5967e-01.*alpha.*beta.^2 - 1.0826e-02.*beta.^3 - 5.4549e-02.*alpha.^4 + 2.5315e-03.*alpha.^3.*beta - 1.4795e-01.*alpha.^2.*beta.^2 + 1.0077e-02.*alpha.*beta.^3 - 4.2589e-02.*beta.^4 ;

%% Cz.beta(alpha,beta,varargin)
Cz.beta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 3.5300e-02 - 1.8145e-02.*alpha + 6.2692e-04.*beta + 2.9108e+00.*alpha.^2 + 1.3280e-03.*alpha.*beta + 2.3916e-02.*beta.^2 - 6.3277e+00.*alpha.^3 - 9.8216e-02.*alpha.^2.*beta + 4.0180e+00.*alpha.*beta.^2 - 1.0264e-04.*beta.^3 - 8.8046e+00.*alpha.^4 + 3.2797e-01.*alpha.^3.*beta - 8.6089e+00.*alpha.^2.*beta.^2 + 6.0262e-04.*alpha.*beta.^3 + 1.6490e-01.*beta.^4 ;
Cz.beta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 3.9843e-02 + 1.8333e-01.*alpha - 4.6875e-03.*beta - 2.6808e-01.*alpha.^2 + 3.8917e-02.*alpha.*beta - 2.9307e-01.*beta.^2 + 1.6848e-01.*alpha.^3 - 8.8197e-02.*alpha.^2.*beta + 3.3285e+00.*alpha.*beta.^2 + 8.5875e-05.*beta.^3 - 3.9043e-02.*alpha.^4 + 5.5971e-02.*alpha.^3.*beta - 2.1481e+00.*alpha.^2.*beta.^2 - 6.7790e-05.*alpha.*beta.^3 + 1.6490e-01.*beta.^4 ;

%% Cm.beta(alpha,beta,varargin)
Cm.beta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 6.7812e-02 - 2.0468e-01.*alpha - 1.5819e-03.*beta - 9.1129e+00.*alpha.^2 + 2.0618e-02.*alpha.*beta - 1.3963e+00.*beta.^2 + 5.4917e+01.*alpha.^3 + 1.6980e-01.*alpha.^2.*beta - 1.6326e+00.*alpha.*beta.^2 + 1.7124e-03.*beta.^3 - 8.3622e+01.*alpha.^4 - 9.8963e-01.*alpha.^3.*beta + 1.1760e+01.*alpha.^2.*beta.^2 - 2.7050e-02.*alpha.*beta.^3 + 1.1643e+00.*beta.^4 ;
Cm.beta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 3.5664e-01 - 2.9539e+00.*alpha - 9.6103e-02.*beta + 7.8023e+00.*alpha.^2 + 5.0970e-01.*alpha.*beta - 5.8900e-01.*beta.^2 - 7.6674e+00.*alpha.^3 - 7.3566e-01.*alpha.^2.*beta - 1.7400e+00.*alpha.*beta.^2 + 5.4340e-03.*beta.^3 + 2.4693e+00.*alpha.^4 + 2.9621e-01.*alpha.^3.*beta + 1.9314e+00.*alpha.^2.*beta.^2 - 4.0285e-02.*alpha.*beta.^3 + 1.1643e+00.*beta.^4 ;

%% Cx.xi(alpha,beta,xi,varargin)
Cx.xi1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 3.7723e-03 + 3.2585e-02.*alpha - 6.0340e-04.*beta + 4.5963e-07.*xi - 2.6170e-01.*alpha.^2 - 1.4504e-03.*alpha.*beta + 4.7915e-07.*alpha.*xi + 4.1911e-02.*beta.^2 + 1.4060e-02.*beta.*xi + 1.3546e-01.*xi.^2 + 1.1946e+00.*alpha.^3 + 9.4009e-02.*alpha.^2.*beta - 5.9831e-05.*alpha.^2.*xi - 1.4749e-01.*alpha.*beta.^2 + 1.6095e-02.*alpha.*beta.*xi - 1.2893e-01.*alpha.*xi.^2 + 1.0245e-04.*beta.^3 - 8.1436e-08.*beta.^2.*xi - 2.1709e-07.*beta.*xi.^2 - 4.1259e-07.*xi.^3 - 2.7321e+00.*alpha.^4 - 3.0757e-01.*alpha.^3.*beta + 1.9575e-04.*alpha.^3.*xi + 9.9276e-01.*alpha.^2.*beta.^2 - 5.7472e-01.*alpha.^2.*beta.*xi + 1.6144e-01.*alpha.^2.*xi.^2 - 6.0153e-04.*alpha.*beta.^3 + 4.7813e-07.*alpha.*beta.^2.*xi + 1.2746e-06.*alpha.*beta.*xi.^2 + 2.4224e-06.*alpha.*xi.^3 - 1.1052e-01.*beta.^4 - 2.3485e-03.*beta.^3.*xi - 4.8208e-02.*beta.^2.*xi.^2 - 1.2406e-02.*beta.*xi.^3 - 4.9920e-01.*xi.^4 ;
Cx.xi2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 7.8972e-02 - 6.1064e-01.*alpha - 1.7730e-03.*beta + 1.0651e-06.*xi + 1.4761e+00.*alpha.^2 + 6.7567e-03.*alpha.*beta - 4.2503e-06.*alpha.*xi + 1.3482e-01.*beta.^2 - 2.5397e-02.*beta.*xi - 8.3482e-02.*xi.^2 - 1.4893e+00.*alpha.^3 - 7.6132e-03.*alpha.^2.*beta + 4.8453e-06.*alpha.^2.*xi - 2.5315e-01.*alpha.*beta.^2 - 1.2616e-02.*alpha.*beta.*xi + 8.9678e-01.*alpha.*xi.^2 - 8.5719e-05.*beta.^3 + 6.8135e-08.*beta.^2.*xi + 1.8163e-07.*beta.*xi.^2 + 3.4520e-07.*xi.^3 + 5.3050e-01.*alpha.^4 + 2.6324e-03.*alpha.^3.*beta - 1.6754e-06.*alpha.^3.*xi + 1.9350e-01.*alpha.^2.*beta.^2 + 2.6409e-02.*alpha.^2.*beta.*xi - 7.1728e-01.*alpha.^2.*xi.^2 + 6.7667e-05.*alpha.*beta.^3 - 5.3786e-08.*alpha.*beta.^2.*xi - 1.4338e-07.*alpha.*beta.*xi.^2 - 2.7250e-07.*alpha.*xi.^3 - 1.1052e-01.*beta.^4 - 2.3485e-03.*beta.^3.*xi - 4.8208e-02.*beta.^2.*xi.^2 - 1.2406e-02.*beta.*xi.^3 - 4.9920e-01.*xi.^4 ;

%% Cz.xi(alpha,beta,xi,varargin)
Cz.xi1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 6.9065e-03 - 3.0066e-02.*alpha - 6.0422e-04.*beta + 4.5964e-07.*xi - 3.2187e-01.*alpha.^2 - 1.4524e-03.*alpha.*beta + 4.7916e-07.*alpha.*xi - 3.3158e-02.*beta.^2 - 1.3952e-01.*beta.*xi + 2.8801e-02.*xi.^2 + 3.0967e+00.*alpha.^3 + 9.4136e-02.*alpha.^2.*beta - 5.9831e-05.*alpha.^2.*xi - 3.2760e-01.*alpha.*beta.^2 + 3.8512e-01.*alpha.*beta.*xi - 2.4725e-01.*alpha.*xi.^2 + 1.0259e-04.*beta.^3 - 8.1437e-08.*beta.^2.*xi - 2.1709e-07.*beta.*xi.^2 - 4.1259e-07.*xi.^3 - 5.2893e+00.*alpha.^4 - 3.0799e-01.*alpha.^3.*beta + 1.9575e-04.*alpha.^3.*xi + 1.4223e-01.*alpha.^2.*beta.^2 + 1.6898e-01.*alpha.^2.*beta.*xi - 7.9508e-01.*alpha.^2.*xi.^2 - 6.0234e-04.*alpha.*beta.^3 + 4.7814e-07.*alpha.*beta.^2.*xi + 1.2746e-06.*alpha.*beta.*xi.^2 + 2.4224e-06.*alpha.*xi.^3 + 1.0922e-01.*beta.^4 + 2.7443e-02.*beta.^3.*xi - 5.8419e-02.*beta.^2.*xi.^2 + 1.9698e-01.*beta.*xi.^3 + 1.1799e-01.*xi.^4 ;
Cz.xi2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 3.8499e-02 - 3.3220e-01.*alpha - 1.7754e-03.*beta + 1.0652e-06.*xi + 7.6159e-01.*alpha.^2 + 6.7658e-03.*alpha.*beta - 4.2503e-06.*alpha.*xi - 4.3183e-02.*beta.^2 - 5.6235e-03.*beta.*xi - 2.7747e-01.*xi.^2 - 4.7996e-01.*alpha.^3 - 7.6235e-03.*alpha.^2.*beta + 4.8454e-06.*alpha.^2.*xi - 3.2569e-01.*alpha.*beta.^2 - 3.4344e-02.*alpha.*beta.*xi + 6.7629e-01.*alpha.*xi.^2 - 8.5835e-05.*beta.^3 + 6.8135e-08.*beta.^2.*xi + 1.8163e-07.*beta.*xi.^2 + 3.4520e-07.*xi.^3 + 5.4264e-02.*alpha.^4 + 2.6360e-03.*alpha.^3.*beta - 1.6754e-06.*alpha.^3.*xi + 2.6225e-01.*alpha.^2.*beta.^2 - 3.2709e-02.*alpha.^2.*beta.*xi - 2.0604e-01.*alpha.^2.*xi.^2 + 6.7758e-05.*alpha.*beta.^3 - 5.3786e-08.*alpha.*beta.^2.*xi - 1.4338e-07.*alpha.*beta.*xi.^2 - 2.7250e-07.*alpha.*xi.^3 + 1.0922e-01.*beta.^4 + 2.7443e-02.*beta.^3.*xi - 5.8419e-02.*beta.^2.*xi.^2 + 1.9698e-01.*beta.*xi.^3 + 1.1799e-01.*xi.^4 ;

%% Cm.xi(alpha,beta,xi,varargin)
Cm.xi1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 3.9411e-02 + 5.9004e-02.*alpha - 6.0496e-04.*beta + 4.5957e-07.*xi + 1.6299e-01.*alpha.^2 - 1.4542e-03.*alpha.*beta + 4.7909e-07.*alpha.*xi + 4.0073e-01.*beta.^2 - 1.3185e-02.*beta.*xi - 5.6380e-01.*xi.^2 + 2.5679e+00.*alpha.^3 + 9.4252e-02.*alpha.^2.*beta - 5.9822e-05.*alpha.^2.*xi - 1.5478e+00.*alpha.*beta.^2 + 6.8362e-02.*alpha.*beta.*xi - 4.3885e-01.*alpha.*xi.^2 + 1.0272e-04.*beta.^3 - 8.1425e-08.*beta.^2.*xi - 2.1706e-07.*beta.*xi.^2 - 4.1253e-07.*xi.^3 - 4.5943e+00.*alpha.^4 - 3.0837e-01.*alpha.^3.*beta + 1.9572e-04.*alpha.^3.*xi - 1.2951e+00.*alpha.^2.*beta.^2 - 4.7948e-02.*alpha.^2.*beta.*xi + 3.5598e-01.*alpha.^2.*xi.^2 - 6.0308e-04.*alpha.*beta.^3 + 4.7807e-07.*alpha.*beta.^2.*xi + 1.2744e-06.*alpha.*beta.*xi.^2 + 2.4221e-06.*alpha.*xi.^3 + 8.4562e-02.*beta.^4 + 8.2458e-02.*beta.^3.*xi + 8.4057e-02.*beta.^2.*xi.^2 + 3.8256e-02.*beta.*xi.^3 + 2.2011e+00.*xi.^4 ;
Cm.xi2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 6.1519e-02 + 6.3582e-01.*alpha - 1.7776e-03.*beta + 1.0650e-06.*xi - 1.6194e+00.*alpha.^2 + 6.7742e-03.*alpha.*beta - 4.2497e-06.*alpha.*xi - 5.1999e-01.*beta.^2 + 5.0694e-02.*beta.*xi - 6.1107e-02.*xi.^2 + 1.4360e+00.*alpha.^3 - 7.6329e-03.*alpha.^2.*beta + 4.8447e-06.*alpha.^2.*xi + 1.7270e+00.*alpha.*beta.^2 - 2.0177e-01.*alpha.*beta.*xi - 2.6120e+00.*alpha.*xi.^2 - 8.5941e-05.*beta.^3 + 6.8126e-08.*beta.^2.*xi + 1.8161e-07.*beta.*xi.^2 + 3.4515e-07.*xi.^3 - 4.3350e-01.*alpha.^4 + 2.6392e-03.*alpha.^3.*beta - 1.6751e-06.*alpha.^3.*xi - 1.2968e+00.*alpha.^2.*beta.^2 + 1.0484e-01.*alpha.^2.*beta.*xi + 1.7265e+00.*alpha.^2.*xi.^2 + 6.7842e-05.*alpha.*beta.^3 - 5.3779e-08.*alpha.*beta.^2.*xi - 1.4336e-07.*alpha.*beta.*xi.^2 - 2.7246e-07.*alpha.*xi.^3 + 8.4562e-02.*beta.^4 + 8.2458e-02.*beta.^3.*xi + 8.4057e-02.*beta.^2.*xi.^2 + 3.8256e-02.*beta.*xi.^3 + 2.2011e+00.*xi.^4 ;

%% Cx.eta(alpha,beta,eta,varargin)
Cx.eta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 8.1468e-03 - 3.5615e-02.*alpha - 5.2634e-03.*beta - 7.4808e-03.*eta - 2.9167e-01.*alpha.^2 + 6.0624e-03.*alpha.*beta + 1.4804e-01.*alpha.*eta + 4.8134e-03.*beta.^2 - 4.1260e-05.*beta.*eta - 1.0173e-01.*eta.^2 + 1.1731e+00.*alpha.^3 + 2.6599e-03.*alpha.^2.*beta - 4.1108e-01.*alpha.^2.*eta + 3.3333e-02.*alpha.*beta.^2 - 1.1125e-02.*alpha.*beta.*eta + 1.0500e-01.*alpha.*eta.^2 + 8.4957e-03.*beta.^3 - 1.0452e-03.*beta.^2.*eta - 4.7324e-03.*beta.*eta.^2 - 3.5844e-02.*eta.^3 ;
Cx.eta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 4.2133e-02 - 2.7600e-01.*alpha - 1.6612e-03.*beta - 1.0578e-02.*eta + 5.4562e-01.*alpha.^2 - 7.8063e-03.*alpha.*beta + 6.1256e-02.*alpha.*eta + 3.8804e-02.*beta.^2 - 4.0077e-03.*beta.*eta - 6.3339e-02.*eta.^2 - 2.9295e-01.*alpha.^3 + 6.4237e-03.*alpha.^2.*beta - 6.3275e-02.*alpha.^2.*eta - 8.7548e-02.*alpha.*beta.^2 + 2.9808e-03.*alpha.*beta.*eta - 3.1541e-02.*alpha.*eta.^2 + 8.4957e-03.*beta.^3 - 1.0452e-03.*beta.^2.*eta - 4.7324e-03.*beta.*eta.^2 - 3.5844e-02.*eta.^3 ;

%% Cz.eta(alpha,beta,eta,varargin)
Cz.eta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 5.6595e-03 + 4.0723e-02.*alpha - 1.6347e-02.*beta - 5.2605e-01.*eta + 2.0510e-01.*alpha.^2 + 7.3804e-02.*alpha.*beta + 4.9999e-02.*alpha.*eta - 1.1846e-02.*beta.^2 + 9.2472e-03.*beta.*eta + 3.6625e-03.*eta.^2 - 1.7695e+00.*alpha.^3 - 6.0034e-02.*alpha.^2.*beta + 8.1498e-01.*alpha.^2.*eta + 3.0598e-01.*alpha.*beta.^2 - 2.7588e-02.*alpha.*beta.*eta + 1.8077e-01.*alpha.*eta.^2 - 8.7928e-03.*beta.^3 + 2.9651e-01.*beta.^2.*eta + 1.2072e-02.*beta.*eta.^2 + 6.4452e-01.*eta.^3 ;
Cz.eta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 1.0014e-01 + 5.5943e-01.*alpha - 6.6225e-03.*beta - 5.9393e-01.*eta - 9.2611e-01.*alpha.^2 + 2.9030e-02.*alpha.*beta + 5.8514e-01.*alpha.*eta + 1.0999e-01.*beta.^2 + 1.0158e-03.*beta.*eta - 3.5189e-02.*eta.^2 + 4.5188e-01.*alpha.^3 - 2.3787e-02.*alpha.^2.*beta - 2.2963e-01.*alpha.^2.*eta - 1.2730e-01.*alpha.*beta.^2 + 1.6850e-03.*alpha.*beta.*eta + 3.1894e-01.*alpha.*eta.^2 - 8.7928e-03.*beta.^3 + 2.9651e-01.*beta.^2.*eta + 1.2072e-02.*beta.*eta.^2 + 6.4452e-01.*eta.^3 ;

%% Cm.eta(alpha,beta,eta,varargin)
Cm.eta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 2.8162e-02 + 3.1738e-02.*alpha - 5.3573e-02.*beta - 1.8513e+00.*eta - 1.7187e-01.*alpha.^2 + 2.0615e-01.*alpha.*beta - 1.7156e-01.*alpha.*eta - 1.9943e-01.*beta.^2 + 2.6264e-02.*beta.*eta - 1.7332e-01.*eta.^2 - 2.9994e-01.*alpha.^3 - 6.7944e-02.*alpha.^2.*beta + 5.1357e+00.*alpha.^2.*eta + 8.1622e-01.*alpha.*beta.^2 - 2.2685e-02.*alpha.*beta.*eta + 8.9728e-01.*alpha.*eta.^2 + 6.4297e-03.*beta.^3 + 6.9331e-01.*beta.^2.*eta + 1.1074e-01.*beta.*eta.^2 + 1.3240e+00.*eta.^3 ;
Cm.eta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 9.0086e-02 - 3.7200e-01.*alpha - 1.9152e-03.*beta - 1.8798e+00.*eta + 4.3560e-01.*alpha.^2 + 2.6226e-03.*alpha.*beta + 1.4697e+00.*alpha.*eta + 8.8271e-03.*beta.^2 + 1.9153e-02.*beta.*eta - 8.6518e-02.*eta.^2 - 1.3931e-01.*alpha.^3 + 2.5174e-03.*alpha.^2.*beta - 3.4114e-01.*alpha.^2.*eta + 7.5616e-02.*alpha.*beta.^2 + 2.6041e-03.*alpha.*beta.*eta + 5.8859e-01.*alpha.*eta.^2 + 6.4297e-03.*beta.^3 + 6.9331e-01.*beta.^2.*eta + 1.1074e-01.*beta.*eta.^2 + 1.3240e+00.*eta.^3 ;

%% Cx.zeta(alpha,beta,zeta,varargin)
Cx.zeta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 1.0868e-02 + 2.6173e-02.*alpha + 5.2633e-03.*beta + 1.9287e-07.*zeta + 3.5140e-01.*alpha.^2 - 6.0623e-03.*alpha.*beta + 3.9465e-07.*alpha.*zeta + 2.2125e-02.*beta.^2 + 9.0830e-02.*beta.*zeta - 4.9030e-03.*zeta.^2 - 2.7311e+00.*alpha.^3 - 2.6616e-03.*alpha.^2.*beta - 2.8751e-05.*alpha.^2.*zeta - 4.2296e-02.*alpha.*beta.^2 - 3.9616e-02.*alpha.*beta.*zeta + 1.3211e-02.*alpha.*zeta.^2 - 8.4957e-03.*beta.^3 - 3.9132e-08.*beta.^2.*zeta - 3.4438e-08.*beta.*zeta.^2 - 2.8813e-08.*zeta.^3 + 6.0591e+00.*alpha.^4 + 5.4926e-06.*alpha.^3.*beta + 9.4066e-05.*alpha.^3.*zeta + 5.2904e-02.*alpha.^2.*beta.^2 + 2.4491e-01.*alpha.^2.*beta.*zeta - 2.1104e-01.*alpha.^2.*zeta.^2 + 1.0743e-08.*alpha.*beta.^3 + 2.2975e-07.*alpha.*beta.^2.*zeta + 2.0219e-07.*alpha.*beta.*zeta.^2 + 1.6917e-07.*alpha.*zeta.^3 - 1.1023e-02.*beta.^4 + 5.1670e-03.*beta.^3.*zeta - 4.2279e-02.*beta.^2.*zeta.^2 - 1.2995e-01.*beta.*zeta.^3 - 5.5814e-02.*zeta.^4 ;
Cx.zeta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 1.6078e-02 - 1.4031e-01.*alpha + 1.6612e-03.*beta + 5.3527e-07.*zeta + 4.7847e-01.*alpha.^2 + 7.8061e-03.*alpha.*beta - 2.0609e-06.*alpha.*zeta + 8.1680e-02.*beta.^2 + 2.0346e-01.*beta.*zeta - 1.2186e-01.*zeta.^2 - 6.6942e-01.*alpha.^3 - 6.4236e-03.*alpha.^2.*beta + 2.3284e-06.*alpha.^2.*zeta - 3.1625e-01.*alpha.*beta.^2 - 4.6272e-01.*alpha.*beta.*zeta + 4.5616e-01.*alpha.*zeta.^2 - 8.4957e-03.*beta.^3 + 3.2740e-08.*beta.^2.*zeta + 2.8813e-08.*beta.*zeta.^2 + 2.4107e-08.*zeta.^3 + 2.9788e-01.*alpha.^4 - 4.7010e-08.*alpha.^3.*beta - 8.0507e-07.*alpha.^3.*zeta + 2.7397e-01.*alpha.^2.*beta.^2 + 3.2509e-01.*alpha.^2.*beta.*zeta - 3.0707e-01.*alpha.^2.*zeta.^2 - 1.2085e-09.*alpha.*beta.^3 - 2.5845e-08.*alpha.*beta.^2.*zeta - 2.2745e-08.*alpha.*beta.*zeta.^2 - 1.9030e-08.*alpha.*zeta.^3 - 1.1023e-02.*beta.^4 + 5.1670e-03.*beta.^3.*zeta - 4.2279e-02.*beta.^2.*zeta.^2 - 1.2995e-01.*beta.*zeta.^3 - 5.5814e-02.*zeta.^4 ;

%% Cz.zeta(alpha,beta,zeta,varargin)
Cz.zeta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 7.7290e-03 - 1.0611e-01.*alpha + 1.6347e-02.*beta + 1.9301e-07.*zeta + 2.5432e-01.*alpha.^2 - 7.3804e-02.*alpha.*beta + 3.9493e-07.*alpha.*zeta + 5.4555e-02.*beta.^2 + 5.4424e-02.*beta.*zeta + 1.1724e-01.*zeta.^2 + 3.4593e+00.*alpha.^3 + 6.0033e-02.*alpha.^2.*beta - 2.8771e-05.*alpha.^2.*zeta + 8.8630e-02.*alpha.*beta.^2 - 3.6667e-01.*alpha.*beta.*zeta + 6.2840e-02.*alpha.*zeta.^2 + 8.7928e-03.*beta.^3 - 3.9159e-08.*beta.^2.*zeta - 3.4462e-08.*beta.*zeta.^2 - 2.8833e-08.*zeta.^3 - 6.3974e+00.*alpha.^4 + 5.5951e-06.*alpha.^3.*beta + 9.4131e-05.*alpha.^3.*zeta - 9.8427e-01.*alpha.^2.*beta.^2 + 1.9337e+00.*alpha.^2.*beta.*zeta - 1.0490e+00.*alpha.^2.*zeta.^2 + 1.0943e-08.*alpha.*beta.^3 + 2.2991e-07.*alpha.*beta.^2.*zeta + 2.0233e-07.*alpha.*beta.*zeta.^2 + 1.6928e-07.*alpha.*zeta.^3 - 9.8876e-02.*beta.^4 - 1.3967e-01.*beta.^3.*zeta - 4.6775e-02.*beta.^2.*zeta.^2 - 4.3133e-02.*beta.*zeta.^3 - 1.2205e-01.*zeta.^4 ;
Cz.zeta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 2.1957e-01 - 1.2475e+00.*alpha + 6.6225e-03.*beta + 5.3564e-07.*zeta + 2.4371e+00.*alpha.^2 - 2.9030e-02.*alpha.*beta - 2.0623e-06.*alpha.*zeta - 8.3481e-02.*beta.^2 + 1.4151e-01.*beta.*zeta - 8.0888e-02.*zeta.^2 - 2.0740e+00.*alpha.^3 + 2.3787e-02.*alpha.^2.*beta + 2.3300e-06.*alpha.^2.*zeta + 3.4468e-01.*alpha.*beta.^2 - 1.6264e-01.*alpha.*beta.*zeta + 5.9123e-01.*alpha.*zeta.^2 + 8.7928e-03.*beta.^3 + 3.2763e-08.*beta.^2.*zeta + 2.8833e-08.*beta.*zeta.^2 + 2.4123e-08.*zeta.^3 + 6.5392e-01.*alpha.^4 - 4.7887e-08.*alpha.^3.*beta - 8.0564e-07.*alpha.^3.*zeta - 1.4910e-01.*alpha.^2.*beta.^2 + 1.0676e-01.*alpha.^2.*beta.*zeta - 4.2239e-01.*alpha.^2.*zeta.^2 - 1.2310e-09.*alpha.*beta.^3 - 2.5863e-08.*alpha.*beta.^2.*zeta - 2.2761e-08.*alpha.*beta.*zeta.^2 - 1.9043e-08.*alpha.*zeta.^3 - 9.8876e-02.*beta.^4 - 1.3967e-01.*beta.^3.*zeta - 4.6775e-02.*beta.^2.*zeta.^2 - 4.3133e-02.*beta.*zeta.^3 - 1.2205e-01.*zeta.^4 ;

%% Cm.zeta(alpha,beta,zeta,varargin)
Cm.zeta1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 1.5589e-02 - 3.8899e-02.*alpha + 5.3573e-02.*beta + 1.9266e-07.*zeta + 3.3455e-01.*alpha.^2 - 2.0614e-01.*alpha.*beta + 3.9421e-07.*alpha.*zeta + 3.2331e-02.*beta.^2 - 3.3994e-01.*beta.*zeta + 2.9821e-01.*zeta.^2 - 7.6531e+00.*alpha.^3 + 6.7943e-02.*alpha.^2.*beta - 2.8719e-05.*alpha.^2.*zeta + 7.2346e-01.*alpha.*beta.^2 - 3.5382e-01.*alpha.*beta.*zeta + 2.1168e-02.*alpha.*zeta.^2 - 6.4297e-03.*beta.^3 - 3.9088e-08.*beta.^2.*zeta - 3.4400e-08.*beta.*zeta.^2 - 2.8780e-08.*zeta.^3 + 2.3281e+01.*alpha.^4 + 5.3460e-06.*alpha.^3.*beta + 9.3961e-05.*alpha.^3.*zeta - 1.0373e+00.*alpha.^2.*beta.^2 + 5.8788e-01.*alpha.^2.*beta.*zeta - 1.0127e+00.*alpha.^2.*zeta.^2 + 1.0455e-08.*alpha.*beta.^3 + 2.2949e-07.*alpha.*beta.^2.*zeta + 2.0197e-07.*alpha.*beta.*zeta.^2 + 1.6898e-07.*alpha.*zeta.^3 - 3.1103e-01.*beta.^4 - 3.0418e-01.*beta.^3.*zeta + 2.1273e-01.*beta.^2.*zeta.^2 + 6.7966e-01.*beta.*zeta.^3 - 2.4709e-01.*zeta.^4 ;
Cm.zeta2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 4.2850e-01 - 3.3245e+00.*alpha + 1.9153e-03.*beta + 5.3467e-07.*zeta + 8.0842e+00.*alpha.^2 - 2.6227e-03.*alpha.*beta - 2.0586e-06.*alpha.*zeta + 2.4457e-01.*beta.^2 - 9.7397e-01.*beta.*zeta + 1.9944e-01.*zeta.^2 - 7.8220e+00.*alpha.^3 - 2.5173e-03.*alpha.^2.*beta + 2.3258e-06.*alpha.^2.*zeta - 3.6303e-01.*alpha.*beta.^2 + 2.6167e+00.*alpha.*beta.*zeta + 1.4951e-01.*alpha.*zeta.^2 - 6.4297e-03.*beta.^3 + 3.2704e-08.*beta.^2.*zeta + 2.8781e-08.*beta.*zeta.^2 + 2.4080e-08.*zeta.^3 + 2.6148e+00.*alpha.^4 - 4.5755e-08.*alpha.^3.*beta - 8.0418e-07.*alpha.^3.*zeta + 1.4226e-01.*alpha.^2.*beta.^2 - 1.9574e+00.*alpha.^2.*beta.*zeta - 2.1992e-01.*alpha.^2.*zeta.^2 - 1.1761e-09.*alpha.*beta.^3 - 2.5816e-08.*alpha.*beta.^2.*zeta - 2.2720e-08.*alpha.*beta.*zeta.^2 - 1.9009e-08.*alpha.*zeta.^3 - 3.1103e-01.*beta.^4 - 3.0418e-01.*beta.^3.*zeta + 2.1273e-01.*beta.^2.*zeta.^2 + 6.7966e-01.*beta.*zeta.^3 - 2.4709e-01.*zeta.^4 ;

%% Cx.qhat(alpha,qhat,varargin)
Cx.qhat1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 8.1545e-03 - 2.9881e-02.*alpha + 8.5131e-01.*qhat + 7.5872e-02.*alpha.^2 + 1.2434e+01.*alpha.*qhat + 5.7133e+02.*qhat.^2 + 1.2730e-01.*alpha.^3 + 2.3482e+01.*alpha.^2.*qhat + 2.1967e+03.*alpha.*qhat.^2 + 2.5292e+03.*qhat.^3 ;
Cx.qhat2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 3.3646e-02 - 2.0866e-01.*alpha + 2.2961e+01.*qhat + 2.4615e-01.*alpha.^2 - 7.6364e+01.*alpha.*qhat + 8.2114e+02.*qhat.^2 - 9.7266e-02.*alpha.^3 + 5.9647e+01.*alpha.^2.*qhat + 1.3083e+03.*alpha.*qhat.^2 + 2.5292e+03.*qhat.^3 ;

%% Cz.qhat(alpha,qhat,varargin)
Cz.qhat1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 8.6906e-03 - 9.8091e-04.*alpha - 3.2946e+01.*qhat - 2.0251e-01.*alpha.^2 - 3.2151e+01.*alpha.*qhat + 1.4018e+03.*qhat.^2 - 2.8198e-01.*alpha.^3 - 8.0578e+01.*alpha.^2.*qhat + 1.2392e+03.*alpha.*qhat.^2 + 2.3459e+03.*qhat.^3 ;
Cz.qhat2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 1.4802e-02 - 3.8512e-01.*alpha - 9.1233e+01.*qhat + 9.7764e-01.*alpha.^2 + 1.9905e+02.*alpha.*qhat + 1.3122e+03.*qhat.^2 - 6.7728e-01.*alpha.^3 - 1.6561e+02.*alpha.^2.*qhat + 1.5579e+03.*alpha.*qhat.^2 + 2.3459e+03.*qhat.^3 ;

%% Cm.qhat(alpha,qhat,varargin)
Cm.qhat1 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) - 2.4401e-02 - 3.2817e-02.*alpha - 4.2378e+01.*qhat + 4.9616e-01.*alpha.^2 - 2.9909e+00.*alpha.*qhat + 7.6502e+02.*qhat.^2 + 9.0869e-01.*alpha.^3 + 2.9585e+01.*alpha.^2.*qhat + 2.2148e+03.*alpha.*qhat.^2 + 2.3704e+03.*qhat.^3 ;
Cm.qhat2 = @(alpha,beta,xi,eta,zeta,phat,qhat,rhat,varargin) 1.2767e-01 - 4.7000e-01.*alpha - 1.1370e+01.*qhat + 4.3183e-01.*alpha.^2 - 1.4597e+02.*alpha.*qhat + 1.0293e+03.*qhat.^2 - 1.7286e-01.*alpha.^3 + 1.4591e+02.*alpha.^2.*qhat + 1.2750e+03.*alpha.*qhat.^2 + 2.3704e+03.*qhat.^3 ;

