%% LABELING BLOCK
var c, N1, N2, N, k1, k2, k, i1, i2, i, mtil, y, ytil, d,               % endogenous variables (quantities)
    w1, w2, rk1, rk2, r, lam,                                           % endogenous variables (prices)
    lna, lnatil, lng, lns, lnmu, lnnu, lnptil,                          % exogenous variables
    tb, tbtil, tbagg, gdp, tbaggout, debtout,                           % additional calculations 
    tbout, tbtilout, lngdp, lnc, lni, lngdpfinal, lngdptil, 
    y_growth_obs, c_growth_obs, i_growth_obs, tby_obs, ptil_dev_obs;    % observables in estimation

varexo ea, eatil, eg, es, emu, enu, eptil; 

%% PARAMETER BLOCK
parameters  alphk, alphm, alphktil, del, phi,                   % technology-related
            betta, theta, thetatil, om, omtil, gam,             % preference-related
            psi, xi,                                            % rate-related
            dbar, rstar, s_share, tb_share, nxtil_share,        % steady state
            abar, rhoa, siga,                                   % exogenous processes
            atilbar, rhoatil, sigatil
            gbar, rhog, sigg,
            sbar, rhos, sigs,
            mubar, rhomu, sigmu,
            nubar, rhonu, signu,
            ptilbar, rhoptil1, rhoptil2, sigptil;

alphk = 0.32;
alphm = 0.05;
alphktil = 0.32;
del = 0.1255;
phi = 6;

betta = 0.9224;
theta = 1.6;            % this give total N approx = 1/3 
thetatil = 1.6;
om = 1.6;
omtil = 1.6;
gam = 2;

psi = 2.8;
xi = 0.199;             % Based on regression results

nxtil_share=0.086;     % Calibration target for ptilbar 
                       % Ratio of net exports of commodities relative to GDP
                       % This includes exports of manufactured goods based on commodities
                       % Excluding this gives 0.0417 (possible robustness value)
                       % ptilbar value will be assigned in steady state file

s_share = 0.0938;      % Calibration target for sbar
                       % Ratio of govt spending to final output
                       % sbar value will be assigned in steady state file

tb_share = -0.00041;   % Calibration target for dbar 
                       % Ratio of trade balance to total GDP
                       % dbar value will be assigned in steady state file

gbar = 1.0117204;      % Calibration of the average growth rate of output 

rstar = 1/betta*gbar^gam - 1;   

abar = 1;
atilbar = 1;
mubar = 1;
nubar = 1;

rhoa = 0.9;
rhoatil = 0.9;
rhog = 0.9;
rhos = 0.9;
rhomu = 0.9;
rhonu = 0.9;
rhoptil1 = 0.95; % calibrated based on data (VAR)
rhoptil2 = -0.13; % calibrated based on data (VAR)

siga = 0;
sigatil = 0;
sigg = 0;
sigs = 0;
sigmu = 0;
signu = 0;
sigptil  = 0.1064; % calibrated based on data (VAR)

%% MODEL BLOCK
model; 

% Model equations (normalised model in levels)
(c - theta/om*N1^om - thetatil/omtil*N2^omtil)^(-gam) = lam;
(c - theta/om*N1^om - thetatil/omtil*N2^omtil)^(-gam)*theta*N1^(om-1) = lam * w1;
(c - theta/om*N1^om - thetatil/omtil*N2^omtil)^(-gam)*theta*N2^(omtil-1) = lam * w2;

w1 = exp(lna) * exp(lng)^(1-alphk-alphm) * (1-alphk-alphm) * k1(-1)^alphk * mtil^alphm * N1^(-alphk-alphm);
w2 = exp(lnatil) * exp(lnptil) * exp(lng)^(1-alphktil) * (1-alphktil) * k2(-1)^alphktil * N2^(-alphktil);

i = i1+i2;
k = k1+k2;
N = N1+N2;
i1 = k1*exp(lng) - (1-del)*k1(-1);
i2 = k2*exp(lng) - (1-del)*k2(-1);

lam = betta*(1+r)*exp(lng)^(-gam)*exp(lnnu(+1)-lnnu)*lam(+1);

r=rstar + psi*(exp(d-dbar) - 1) - xi*ptil_dev_obs + (exp(exp(lnmu) - 1) - 1);

%d/(1+r)*exp(lng) = d(-1) - y - exp(lnptil)*ytil + c + i + exp(lns) + exp(lnptil)*mtil + phi/2*(k/k(-1)*exp(lng) - gbar)^2 * k(-1);
d/(1+r)*exp(lng) = d(-1) - w1*N1 - w2*N2 - rk1*k1 - rk2*k2 + c + i + exp(lns) + phi/2*(k/k(-1)*exp(lng) - gbar)^2 * k(-1);

lam * (1+phi*(k/k(-1)*exp(lng)-gbar))  = betta * exp(lng)^-gam * exp(lnnu(+1)-lnnu)* lam(+1) *(rk1(+1) + 1 - del
                            + phi*(k(+1)/k*exp(lng(+1))-gbar)*(k(+1)/k)*exp(lng(+1)) 
                            - phi/2*(k(+1)/k*exp(lng(+1))-gbar)^2);

rk1 = exp(lna) * alphk * exp(lng)^(1-alphk-alphm) * k1(-1)^(alphk-1) * mtil^(alphm) * N1^(1-alphk-alphm);

lam * (1+phi*(k/k(-1)*exp(lng)-gbar))  = betta * exp(lng)^-gam * exp(lnnu(+1)-lnnu)* lam(+1) *(rk2(+1) + 1 - del
                            + phi*(k(+1)/k*exp(lng(+1))-gbar)*(k(+1)/k)*exp(lng(+1)) 
                            - phi/2*(k(+1)/k*exp(lng(+1))-gbar)^2);

rk2 = exp(lnptil) * exp(lnatil) * alphktil * exp(lng)^(1-alphktil) * k2(-1)^(alphktil-1) * N2^(1-alphktil);

exp(lnptil) = exp(lna) * alphm * exp(lng)^(1-alphk-alphm) * k1(-1)^alphk * mtil^(alphm-1) * N1^(1-alphk-alphm);

y = exp(lna) * k1(-1)^alphk * mtil^alphm * (exp(lng)*N1)^(1-alphk-alphm);
ytil = exp(lnatil) * k2(-1)^alphktil * (exp(lng)*N2)^(1-alphktil);

% Additional calculations
tb = y - c - i - exp(lns) - phi/2*(k/k(-1)*exp(lng) - gbar)^2 * k(-1);
tbtil = exp(lnptil)*(ytil-mtil);
tbagg = tb+tbtil;
gdp = y+exp(lnptil)*ytil-exp(lnptil)*mtil;
tbaggout = (tb+tbtil) / gdp;
tbout = tb / gdp;
tbtilout = tbtil / gdp;
debtout = d(-1) / gdp;
lngdp = log(gdp);
lnc = log(c);
lni = log(i);
lngdpfinal = log(y);
lngdptil = log(exp(lnptil)*ytil-exp(lnptil)*mtil);

% Exogenous processes
lna = (1-rhoa)*log(abar) + rhoa*lna(-1) + ea;
lnatil = (1-rhoatil)*log(atilbar) + rhoatil*lnatil(-1) + eatil;
lng = (1-rhog)*log(gbar) + rhog*lng(-1) + eg;
lns = (1-rhos)*log(sbar) + rhos*lns(-1) - es;
lnmu = (1-rhomu)*log(mubar) + rhomu*lnmu(-1) - emu;
lnnu = (1-rhonu)*log(nubar) + rhonu*lnnu(-1) + enu;
lnptil = (1-rhoptil1+rhoptil2)*log(ptilbar) + rhoptil1*lnptil(-1) - rhoptil2*lnptil(-2) + eptil;

% Measurement equations (add M.E. in these equations)
% '_obs' variables are observable in the data
y_growth_obs = log(gdp) - log(gdp(-1)) + lng(-1);   % this is the total empirical growth rate
c_growth_obs = log(c) - log(c(-1)) + lng(-1);       % this is the total empirical growth rate
i_growth_obs = log(i) - log(i(-1)) + lng(-1);       % this is the total empirical growth rate 
tby_obs = tbaggout;                               % this is the empirical ratio
ptil_dev_obs = lnptil - log(ptilbar);               % this is the empirical log-deviation from trend
                                                    % (data will need to be mean zero and in logpoints)

end;


// SEE STEADY STATE FILE !!!

steady;

check;


%% RANDOM SHOCKS BLOCK
shocks;
var ea; stderr siga;
var eatil; stderr sigatil;
var eg; stderr sigg;
var es; stderr sigs;
var emu; stderr sigmu;
var enu; stderr signu;
var eptil; stderr sigptil;
end;

estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF

xi,0.199,0,1,NORMAL_PDF,0.199,0.045;
%psi,0.01,0,3,GAMMA_PDF,0.001,0.025;
psi,2.8,0,6,NORMAL_PDF,2.8,0.5;

rhoa,0.5,0.01,0.99,BETA_PDF,0.5,0.20;
rhoatil,0.5,0.01,0.99,BETA_PDF,0.5,0.20;
rhog,0.5,0.01,0.99,BETA_PDF,0.5,0.20;
rhonu,0.5,0.01,0.99,BETA_PDF,0.5,0.20;
rhos,0.5,0.01,0.99,BETA_PDF,0.5,0.20;
rhomu,0.5,0.01,0.99,BETA_PDF,0.5,0.20;
rhoptil1,0.8,0.01,0.99,BETA_PDF,0.8,0.10;
rhoptil2,0.15,0.01,0.99,BETA_PDF,0.15,0.10;

stderr ea,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;
stderr eatil,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;
stderr eg,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;
stderr enu,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;
stderr es,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;
stderr emu,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;
stderr eptil,0.1,0.001,0.5,INV_GAMMA_PDF,0.1,2;

end;

% COMMODITY PRICE UNOBSERVED VS. OBSERVED
varobs y_growth_obs c_growth_obs i_growth_obs tby_obs;
%varobs y_growth_obs c_growth_obs i_growth_obs tby_obs ptil_dev_obs;

% Estimation
% FULL SAMPLE VS. POST 1950 SAMPLE
estimation(datafile=DTDATA,mode_compute=6,mh_replic=10000000,mh_nblocks=1,mh_jscale=0.50,mh_drop=0.25);
%estimation(datafile=DTDATA,first_obs=50,mode_compute=4,mh_replic=100000,mh_nblocks=1,mh_jscale=0.50,mh_drop=0.5);

% Generate results, at the mode
stoch_simul(conditional_variance_decomposition=[1,2,3,4,5,6,7,8,9,10],irf=10) lngdp lnc lni tby_obs y_growth_obs c_growth_obs i_growth_obs tby_obs ptil_dev_obs;
shock_decomposition y_growth_obs c_growth_obs i_growth_obs tby_obs ptil_dev_obs;
