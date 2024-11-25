function Theta = NM_initialization(data, m, W)
% data: the data to be fitting including complete, right-censored
% m: number of mixing components
% Theta: estimated initial gamma parameters for complete and censored data

if W == "W1"
    V = data(:,2)+data(:,3);
elseif W == "W2"
    V = data(:,2)+data(:,3)+data(:,4);
end

name_cluster = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]; % maximum support 10 clusters/gamma components
% ===== starting values for Nelder-Mead by moments
idx = kmeans(V,m);
dict = struct;          % dict: a struct to store data in k clusters 
for i = [name_cluster(1:m);1:m]  % i(1) the name of cluster, i(2) the value of the cluster
    dict.(i(1)) = V(idx==str2num(i(2)));
end

mixing_proportions = [];   % mixing proportions of the k gamma
for i = name_cluster(1:m)
    mixing_proportions = [mixing_proportions; length(dict.(i))];
end
mixing_proportions = mixing_proportions / sum(mixing_proportions);

Theta = zeros(m, 3);  % parameters of gammas, mixing, shape and scale
for i = [name_cluster(1:m);1:m]
    data_cluster = dict.(i(1));
    fir = sum(data_cluster) / length(data_cluster);
    sec = sum(data_cluster .* data_cluster) / length(data_cluster);
    rate = fir / (sec - fir * fir);   % rate
    shape = rate * fir;           % shape
    Theta(str2num(i(2)),:) = [mixing_proportions(str2num(i(2))), shape, 1/rate];
end