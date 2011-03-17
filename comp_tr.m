function [rating, etime] = comp_tr(distribution)
% This script creates nice pictures

n_children = 7;
n_stages   = 3;
n_rating   = 10;
n_rating_scen = 3000;

% distribution information
sigma = 0.5;
mu = 0;
ln_init = 1;

% 1. Datenabhaengigkeit
testscen = 1000;
%testscen = [500 1000 2000];
ntestscen = length(testscen);
stime = zeros(ntestscen,1,'uint64');
etime = zeros(ntestscen,1);
rating = zeros(ntestscen, n_rating);
for kk=1:ntestscen
    % create distribution and tree
    [xi, p] = generate_scenarios(testscen(kk), n_stages, distribution, ln_init, mu, sigma);
    stime(kk) = tic;
    %tr_kmedoid = treeDEkmedoids(xi', p, n_children, 1);
    tr_kmedoid = treeDEsw(xi, p, n_children, 1);
    etime(kk) = toc(stime(kk));
    % test with a lot of distributions
    for jj=1:n_rating
        [nu, q] = generate_scenarios(testscen(kk), n_stages, ...
                                     distribution, ln_init, mu, ...
                                     sigma);
        rating(kk, jj) = tr_kmedoid.kantorovich(nu, q, 1);
    end
end
figure
errorbar(testscen', mean(rating,2), sqrt(var(rating'))')
figure
plot(testscen',etime)
keyboard
