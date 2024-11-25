clear all; clc; close all;
warning off;
format long g   % not use scientific notation
% ==================================== settings ========================
Path_data = "data/window3years/";

nbins = 50;         % used in stage 6 for plotting histogram
fontSize = 16;
max_hours = 12*31*24*2;        % 2 years
    
% ms = [2,3,4,5,6,7];   % the number of mixing components
ms = [2,3,4,5];   % the number of mixing components

% stage0: prepare original data to 5-columns, ID, y1, y2, y3, tag=1-8 
% stage1: short-runs 
% stage2: select the best model under each sliding window with different m based on BIC
          % use the average likelihood of the three windows (keep all models with the same mixing components) 
% stage3: long-runs of Nelder-Mead based on models from stage 2
% stage4: the number of complete, truncated and censored data within the three windows

n_iters_stage1 = 1000;   % the number of iterations in Nelder-Mead (short run)
n_iters_stage3 = 2000;  % the number of iterations in Nelder-Mead (long run)

Stages = [1,2,3];
stage4_out = [];
nwindow = 3; % the length of the observation window, 12,i.e., 1 year
startT = clock;
for stage = Stages
    if stage == 0
        % ================== transform data to 5 columns, ID, y1,y2,y3 and tag            
        years = [
         "01/01/2013 00:00"...,
         "01/01/2014 00:00"...,
         "01/01/2015 00:00"...,
         "01/01/2016 00:00"...,
        ];

        % windows for concept drift detection
        W1 = [];
        W2 = [];
        for Ts1 = 1:nwindow
            temp_w1 = [years(Ts1), years(Ts1+1)];
            temp_w2 = [years(Ts1+1), years(Ts1+nwindow)];
            W1 = [W1;temp_w1];
            W2 = [W2;temp_w2];
        end

        data_file = strcat('data/CSDM2.csv');
        data = readmatrix(data_file, 'OutputType', 'string');   % data to be fitted
        for window = 1:nwindow
            W1_single = W1(window,:);
            W2_single = W2(window,:);
            new_data = data_extraction(data, W1_single, W2_single);
            new_data2 = sortrows(new_data,5);

            T = array2table(new_data2);
            T.Properties.VariableNames = {'ID', 'y1', 'y2', 'y3', 'tag'};
            output_file = strcat(Path_data, num2str(window), '.csv');
            writetable(T, output_file);
        end
    end
end   
for sliding_window = 1:nwindow
    sliding_window
    
    W_file = strcat(Path_data, num2str(sliding_window), '.csv');
    data = readmatrix(W_file, 'OutputType', 'double');   % data to be fitted
    [n,~] = size(data);
    % ======================== Nelder-Mead initialization ===========
    for stage = Stages 
        if stage == 1
            % create folders for saving model parameters and images for all sliding windows
            Path_slidingwindow = strcat('out3Years/models/', num2str(sliding_window));
            Path_W = strcat(Path_slidingwindow, '/W/');    
            Path_W1 = strcat(Path_slidingwindow, '/W1/');
            Path_W2 = strcat(Path_slidingwindow, '/W2/');
            if exist(Path_slidingwindow, 'dir')==0 
                mkdir (Path_slidingwindow)
                mkdir ( Path_W );
                mkdir ( Path_W1 );
                mkdir ( Path_W2 );
            end
            
            for m = ms
                param_W1 = NM_initialization(data, m, 'W1'); % mixing, shape and scale, the whole window W for initialization
                param_W2 = NM_initialization(data, m, 'W2');
                
                % remove samples if its log-likelohood is nan/inf/0 with initial parameters (ensure optimization can start)
                [data_good, ~, ~] = NM_loglikelihood_gooddata(data, param_W1, param_W2);
                % optimize models (short runs)      
                param_shortW1 = NM_optimization(data_good, param_W1, 'W1', n_iters_stage1);
                param_shortW2 = NM_optimization(data_good, param_W2, 'W2', n_iters_stage1);
                NM_saving(param_shortW1, Path_W1, 'short', m);
                NM_saving(param_shortW2, Path_W2, 'short', m);
            end
        end
        
        if stage == 2
            % create folders for saving best models as short run models
            Path_slidingwindow = strcat('out3Years/short run/', num2str(sliding_window));
            if exist(Path_slidingwindow, 'dir')==0 
                mkdir (Path_slidingwindow)
            end
            
            % the path of models to be estimated   
            Path_W1 = strcat('out3Years/models/', num2str(sliding_window), '/W1/');
            Path_W2 = strcat('out3Years/models/', num2str(sliding_window), '/W2/');  
            target_Path = strcat('out3Years/short run/', num2str(sliding_window), '/' );    
            
            [L1_temp, L2_temp] = deal([],[]);
            for m = ms
                % best model of W under different m
                W_file1 = strcat(Path_W1, num2str(m), '.csv');
                W_file2 = strcat(Path_W2, num2str(m), '.csv');   
                parameters_W1 = readmatrix(W_file1, 'OutputType', 'double');   % data to be fitted
                parameters_W2 = readmatrix(W_file2, 'OutputType', 'double');   % data to be fitted

                param_logW1 = NM_param2log(parameters_W1, m);        % transform param to its logrithm form
                param_logW2 = NM_param2log(parameters_W2, m);        % transform param to its logrithm
                    
                [data_good, ~, ~] = NM_loglikelihood_gooddata(data, parameters_W1, parameters_W2);
                [W1_fg, ~, ~] = gammix_in(param_logW1, data_good, "W1");
                [W2_fg, ~, ~] = gammix_in(param_logW2, data_good, "W2");
                
                BIC1_temp = log(n) * (3*m-1) - 2*W1_fg;
                BIC2_temp = log(n) * (3*m-1) - 2*W2_fg; 
                
                [L1_temp, L2_temp] = deal([L1_temp, BIC1_temp],[L2_temp, BIC2_temp]);
            end
            [LL_W1, ind_W1] = max(L1_temp);
            [LL_W2, ind_W2] = max(L2_temp);
               
            Best_fileW1 = strcat( num2str(ms(ind_W1)), ".csv"); 
            Best_fileW2 = strcat( num2str(ms(ind_W2)), ".csv"); 
            best_files = [Best_fileW1, Best_fileW2];
            
            % copy best models to folder best
            NM_copy(best_files, Path_W1, Path_W2, target_Path)
        end
        
        if stage == 3
            % create folders for saving long-run models 
            Path_longrun = strcat('out3Years/long run/', num2str(sliding_window), '/');
            if exist(Path_longrun, 'dir')==0 
                mkdir (Path_longrun)
            end

            % load models
            parameters_Path = strcat('out3Years/short run/', num2str(sliding_window), '/' );                                          
            parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
            parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double');   
            
            % select good data
            [data_good, ~, ~] = NM_loglikelihood_gooddata(data, parameters_W1, parameters_W2);
            % optimize models            
            param_longW1 = NM_optimization(data_good, parameters_W1, 'W1', n_iters_stage3);
            param_longW2 = NM_optimization(data_good, parameters_W2, 'W2', n_iters_stage3);

            % remove mixing components with extremely small mixing proportions
            param_longW1 = NM_parameters_modification(param_longW1);
            param_longW2 = NM_parameters_modification(param_longW2); 
            
            NM_saving(param_longW1, Path_longrun, 'long', 'W1')
            NM_saving(param_longW2, Path_longrun, 'long', 'W2')
        end
        
        if stage == 4
            Path_out = 'out3Years/summary.csv';
            tags = data(:,5);

            % ===================== window W1            
            tags_W1_complete = tags(tags == 1 |tags == 2);
            tags_W1_RC = tags(tags == 3 | tags == 4 | tags == 5 | tags == 6);

            [n_W1_complete, ~] = size(tags_W1_complete); 
            [n_W1_RC, ~] = size(tags_W1_RC); 
            
            % ===================== window W2           
            tags_W2_complete = tags(tags == 3 | tags == 4 | tags == 7);
            tags_W2_RC = tags(tags == 5 | tags == 6 | tags == 8);
            
            [n_W2_complete, ~] = size(tags_W2_complete); 
            [n_W2_RC, ~] = size(tags_W2_RC);  
            
            % ===================== Table
            n_W1 = n_W1_complete + n_W1_RC;
            n_W2 = n_W2_complete + n_W2_RC;       
            stage4_out = [stage4_out; [sliding_window, ...
                n_W1_complete, round(n_W1_complete/n_W1,3), n_W1_RC, round(n_W1_RC/n_W1,3),  ...
                n_W2_complete, round(n_W2_complete/n_W2,3), n_W2_RC, round(n_W2_RC/n_W2,3)]];
            T = array2table(stage4_out);
            writetable(T, Path_out);            
        end

        if stage == 6
            parameters_Path = strcat('out3Years/long run/', num2str(sliding_window), '/' );  
            parameters_W1 = readmatrix( strcat(parameters_Path, 'W1.csv'), 'OutputType', 'double');   
            parameters_W2 = readmatrix( strcat(parameters_Path, 'W2.csv'), 'OutputType', 'double'); 
            parameters_W1(:,4) = parameters_W1(:,2) .* parameters_W1(:,3);   % calculate means
            parameters_W2(:,4) = parameters_W2(:,2) .* parameters_W2(:,3);   % calculate means
            
            tags = data(:,5);
            data_W1_C  = data(tags == 1 | tags == 2, :);
            data_W1_RC  = data(tags == 3 | tags == 4 | tags == 5 | tags == 6, :);
            data_W2_C = data(tags == 3 | tags == 4 | tags == 7, :);
            data_W2_RC = data(tags == 5 | tags == 6 | tags == 8, :);
            
            duration_W1_C  = data_W1_C(:,2) + data_W1_C(:,3);
            duration_W1_RC  = data_W1_RC(:,2) + data_W1_RC(:,3);
            duration_W2_C  = data_W2_C(:,2) + data_W2_C(:,3) + data_W2_C(:,4);
            duration_W2_RC  = data_W2_RC(:,2) + data_W2_RC(:,3) + data_W2_RC(:,4);
            
            duration_W1_C = duration_W1_C';
            duration_W1_RC = duration_W1_RC';
            duration_W2_C = duration_W2_C';
            duration_W2_RC = duration_W2_RC';
            
            minY = 0; %there are no values lower than this
            maxY = 12000; %there are no values higher than this
            stride = 300;
            binEdges = 0:stride:maxY;
            Nedges = size(binEdges);
            Nedges = Nedges(2);
            
            samples = 1:10:maxY;
            y1 = pdf(parameters_W1(:,[1,2,4]), samples);
            y2 = pdf(parameters_W2(:,[1,2,4]), samples);
            
            % =============== figure 1 for W1
            % Add in extrema placeholders to adjust bins to a common scale 
            h = histogram(duration_W1_C, binEdges,'Normalization','pdf');
            countsA1 = h.Values;            
            h = histogram(duration_W1_RC, binEdges,'Normalization','pdf');
            countsB1 = h.Values;
            
            newX = binEdges(1:Nedges-1)+stride/2;
            newY = [countsA1', countsB1']';
            
            fig1=figure;
            % Plot histograms
            b = bar(newX, newY, 'stacked');
            b(1).FaceColor = 'r';
            b(2).FaceColor = 'b';
            hold on;
            plot(samples, y1, 'g-', 'DisplayName', 'W1', 'lineWidth', 1.5);
            hold off;

            % Labels
            set(gca, 'XLim', [minY maxY], 'FontSize', fontSize)
            xlabel('hours')
            ylabel('pdf')
            legend({'Complete', 'Right-censored'})
            
            output_png = strcat('out3Years/PDF/', num2str(sliding_window), 'W1.png' );
            print(fig1,'-r500','-dpng',output_png);
            close(fig1)
            
            % =============== figure 2 for W2
            % Add in extrema placeholders to adjust bins to a common scale
            h = histogram(duration_W2_C, binEdges,'Normalization','pdf');
            countsA2 = h.Values;            
            h = histogram(duration_W2_RC, binEdges,'Normalization','pdf');
            countsB2 = h.Values;

            newX = binEdges(1:Nedges-1)+stride/2;
            newY = [countsA2', countsB2']';
            
            fig2=figure;
            % Plot histograms
            b = bar(newX, newY, 'stacked');
            b(1).FaceColor = 'r';
            b(2).FaceColor = 'b';
            hold on;
            plot(samples, y2, 'g-', 'DisplayName', 'W2', 'lineWidth', 1.5);
            hold off;

            % Labels
            set(gca, 'XLim', [minY maxY], 'FontSize', fontSize)
            xlabel('hours')
            ylabel('pdf')
            legend({'Complete', 'Right-censored'})
            output_png = strcat('out3Years/PDF/', num2str(sliding_window), 'W2.png' );
            print(fig2,'-r500','-dpng',output_png);  
            close(fig2)
        end
    end
end

endT = clock;
Texcution = etime(endT,startT);