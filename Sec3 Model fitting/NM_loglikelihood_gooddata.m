function [data_good, Lc1, Lc2] = NM_loglikelihood_gooddata(data, parameters_W1, parameters_W2)
    % select good data for log-likelihood calculation
    
    [m1, ~] = size(parameters_W1);
    [m2, ~] = size(parameters_W2);
    IDSS = data(:,1);

    paramW1_log = NM_param2log(parameters_W1, m1);              % transform param to its logrithm form
    paramW2_log = NM_param2log(parameters_W2, m2);              % transform param to its logrithm form

    [~, ~, W1_L] = gammix_in(paramW1_log, data, "W1");
    [~, ~, W2_L] = gammix_in(paramW2_log, data, "W2");

    ind_W1 = NM_Nan_Inf_0(IDSS, W1_L);% IDs of data that with log-likelihood of Nan, inf or 0
    ind_W2 = NM_Nan_Inf_0(IDSS, W2_L);
    ID_bad = [ind_W1; ind_W2]; 

    IDs_good = setdiff(IDSS,ID_bad).';    % case IDs of good instances
    ind_good = [];                        % indices of good instances
    for id = IDs_good
        temp = find(IDSS==id);
        ind_good = [ind_good, temp];
    end

    data_good = data(ind_good,:); 
    [Lc1, ~, ~] = gammix_in(paramW1_log, data_good, "W1");
    [Lc2, ~, ~] = gammix_in(paramW2_log, data_good, "W2");
    [Lc1,Lc2] = deal( roundn(Lc1,-2),roundn(Lc2,-2) );
end
               
    