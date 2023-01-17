%Ising Model
format long;
% clear;
%Parameters
N = {10, 10, 10};%Number of martices
Nensemble = 1e5;%Number of systems in an ensemble.
% beta_J = 5;
% J = 1;

dim = numel(N);
Ndot = 1;%Total number of particles
for k = 1:dim
    Ndot = Ndot.*N{k};
end
%Spin array
% Sigma = ones(Ndot, dim) - 2.*(rand(Ndot, dim) <= 0.5);
Sigma = ones(Ndot, 1) - 2.*(rand(Ndot, 1) <= 0.5);
%Ensemble array
% Ensemble = zeros(Nensemble, dim);
Ensemble = zeros(Nensemble, 1);
Ensemble(1, :) = mean(Sigma, 1);
EnsembleH = zeros(Nensemble, 1);
%Cofficient for array searching
indcff = ones(dim, 1);
for k = 2:dim
    indcff(k) = N{k-1}.*indcff(k-1);
end
GenIndVec = cell(dim, 2);
IndArray = cell(1, dim);
for k = 1:dim
    IndArray{k} = 1;
end
for k = 1:dim
    IndArray{k} = N{k};
    GenIndVec{k, 1} = true(IndArray{:});
    GenIndVec{k, 2} = GenIndVec{k, 1};
    GenIndVec{k, 1}(end) = false;
    GenIndVec{k, 2}(1) = false;
    IndArray{k} = 1;
end
IndRepmat = N;
IndH = cell(dim, 2);
for k = 1:dim
    IndRepmat{k} = 1;
    IndH{k, 1} = reshape(repmat(GenIndVec{k, 1}, IndRepmat{:}), 1, Ndot);
    IndH{k, 2} = reshape(repmat(GenIndVec{k, 2}, IndRepmat{:}), 1, Ndot);
    IndRepmat{k} = N{k};
end
%Initialize EnsembelH
for k = 1:dim
    EnsembleH(1) = EnsembleH(1) - sum(Sigma(IndH{k, 1}, :).*Sigma(IndH{k, 2}, :), 'all');
end
%Metropolis
MinInd = zeros(1, dim);
MinInd(1) = 1;
MaxInd = ones(1, dim);
MaxInd(1) = N{1};
for k = 2:dim
    MaxInd(k) = N{k} - 1;
end
NewInd = zeros(1, dim);
for k = 2:Nensemble
    NewInd(1) = randi(N{1});
    for j = 2:dim
        NewInd(j) = randi(N{j}) - 1;
    end
    NewJ = 1;
    OneDNewInd = NewInd*indcff;
    Sigma(OneDNewInd) = - Sigma(OneDNewInd);
    %Calculate Energy Change
    DeltaH = 0;
    VicNewInd = NewInd;
    for j = 1:dim
        if NewInd(j) ~= MinInd(j)
            VicNewInd(j) = NewInd(j) - 1;
            VicOneDNewInd = VicNewInd*indcff;
            DeltaH = DeltaH - Sigma(VicOneDNewInd).*Sigma(OneDNewInd);
            VicNewInd(j) = NewInd(j);
        end
        if NewInd(j) ~= MaxInd(j)
            VicNewInd(j) = NewInd(j) + 1;
            VicOneDNewInd = VicNewInd*indcff;
            DeltaH = DeltaH - Sigma(VicOneDNewInd).*Sigma(OneDNewInd);
            VicNewInd(j) = NewInd(j);
        end
    end
    DeltaH = DeltaH.*2;
    %Acceptable?
    if rand < exp(-beta_J.*DeltaH)
        EnsembleH(k) = EnsembleH(k-1) + DeltaH;
    else
        Sigma(OneDNewInd, :) = - Sigma(OneDNewInd, :);
        EnsembleH(k) = EnsembleH(k-1);
    end
    Ensemble(k, :) = mean(Sigma, 1);
end
EnsembleH_calc = -beta_J.*EnsembleH;
EnsembleH_calc = EnsembleH_calc - max(EnsembleH_calc, [], 'all');
Z_calc = exp(EnsembleH_calc);
EnsAve = mean(Ensemble.*Z_calc, 'all')./mean(Z_calc, 'all');
fprintf('beta*J = %e     Sigma_ave = %e\n', beta_J, EnsAve);
