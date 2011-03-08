% script for testing different algorithms against historic stock
% data
addpath('/home/hans/src/clp/build/');

stock_id_dax    = '^GDAXI'; % DAX
                     
stock_id_basf   = 'BASFY.PK'; % BASF
stock_id_brent  = '^DJUBSCOT'; % BRENT OIL 
stock_id_eurusd = 'EURUSD=X'; % 1 EUR in USD exchange rate


[hist_date_dax, hist_high_dax, hist_low_dax, hist_open_dax, ...
    hist_close_dax, hist_vol_dax] = get_hist_stock_data(stock_id_dax);
%[hist_date_brent, hist_high_brent, hist_low_brent, hist_open_brent, ...
%    hist_close_brent, hist_vol_brent] = get_hist_stock_data(stock_id_brent);
%[hist_date_eurusd, hist_high_eurusd, hist_low_eurusd, hist_open_eurusd, ...
%    hist_close_eurusd, hist_vol_eurusd] = get_hist_stock_data(stock_id_eurusd);

% ---------------------------------------
% SELECT NUMBER OF CHILDREN
nchildren = 3; 
% ---------------------------------------
 
% 4 Quartalsweise scenarien
nstages = 4;
sdat = hist_close_dax;
step = 90;
idx = 0:step:nstages*step;
nscen = length(sdat)-idx(end);
xi = zeros(nscen, nstages+1);
for kk=1:nscen
    xi(kk,:) = sdat(kk+idx)/sdat(kk);
end
xi = xi(:,2:end); % first value is always one

ngenscen = 1000;                 
ntestcen = nscen - ngenscen;

%% randomly select generating scenarios

ri = randi(nscen, ngenscen,1);

genxi   = xi(ri, :);
ngenxi  = size(genxi,1);
genp    = 1/ngenxi*ones(ngenxi,1);
testxi  = xi(setdiff(1:nscen,ri),:);
ntestxi = size(testxi,1);
testp   = 1/ntestxi*ones(ntestxi,1);

%% Calculate Trees

% geneticDE
geneticDE_tr   = geneticDE(genxi, genp, nchildren);
geneticDE_scen = geneticDE_tr.tree2scen;
geneticDE_p    = geneticDE_tr.p;

% backwardkmeans
%[bwkmeans_scen, bwkmeans_p] = backwardtreeKmediods(genxi', genp, ...
%                                                  nchildren);


%% Asses quality of trees

Dk_geneticDE = kantorovich(testxi, testp, geneticDE_scen(2:end,:)', geneticDE_p);
%Dk_bwkmeans  = kantorovich(testxi, testp, bwkmeans_scen, ...
%                           bwkmeans_p);
