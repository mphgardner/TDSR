function results = sim_TDSR(sim)
    
    % Simulate experimental paradigms.
    %
    % USAGE: results = sim_TDSR(sim)
    %
    % INPUTS:
    %   sim - string specifying simulation:
    %           'sharpe17_opto' - activation experiment from Sharpe et al (2017)
    %           'sharpe17_deval' - devaluation experiment from Sharpe et al (2017)
    %           'sharpe17_inhib' - inhibition experiment from Sharpe et al (2017)
    %           'takahashi17_identity' - identity change experiment from Takahashi et al (2017)
    %           'chang17_identity' - identity unblocking experiment from Chang et al (2017)
    %
    % OUTPUTS:
    %   results - see linearTDSR.m
    %
    % Sam Gershman, Dec 2017, edited by Matt Gardner April 2018
    
    switch sim
        
        case 'sharpe17_opto'
            
            % A = 1
            % C = 2
            % D = 3
            % E = 4
            % F = 5
            % X = 6
            % food = 7
            
            A_X = [1 0 0 0 0 0 0; 1 0 0 0 0 1 0];    %24
            EF_X = [0 0 0 1 1 0 0; 0 0 0 1 1 1 0];   %8
            AD_X = [1 0 1 0 0 0 0; 1 0 1 0 0 1 0];   %8
            AC_X = [1 1 0 0 0 0 0; 1 1 0 0 0 1 0];   %8
            X = [0 0 0 0 0 1 0; 0 0 0 0 0 1 1];      %4*24
            F = [0 0 0 0 1 0 0];
            D = [0 0 1 0 0 0 0];
            C = [0 1 0 0 0 0 0];
            
            %This sets the alphas for the SR and U learning:
            alpha = [0.06 0.03];
            
            %this function provides the logicals of the first state of a stimulus
            st1 = @(x,Y) all(Y == x(1,:),2);
            
            %this function provides the second state of a stimulus
            st2 = @(x) [false; x(1:end - 1)];
            %st2 = @(x) all(M == x(2,:),2);
            
            %stimuli are randomized within each execution of the TDSR model. 
            for i = 1:100
                
                %this randomizes trial order for each stage
                s1 = randTrials(24,A_X);
                s2 = randTrials(8,EF_X,AD_X,AC_X);
                s3 = randTrials(4*24,X);
                s4 = randTrials(6,F,D,C);
                
                %Matrix with all stages included
                M = [s1;s2;s3;s4];
                
                %This keeps track of trials from each stage
                sM = [ones(size(s1,1),1);2*ones(size(s2,1),1);3*ones(size(s3,1),1);4*ones(size(s4,1),1)];
                
                %logicals of the first state of each stimulus. The logicals
                %are conditioned on the stage in which the cue occured. 
                a_x = st1(A_X,M) & sM == 1;
                ef_x = st1(EF_X,M) & sM == 2; ad_x = st1(AD_X,M) & sM == 2; ac_x = st1(AC_X,M) & sM == 2;
                x = st1(X,M) & sM == 3;
                f = st1(F,M); d = st1(D,M); c = st1(C,M);

                %rewarded states
                r = double(st2(x));

                opto = 1*(st2(ac_x) | st2(st2(ad_x)));
                               
                %This sets the initial value of the US
                ui = [0 0 0 0 0 0 0.5]';
                
                %this runs the model for the current iteration and saves
                %the results in the following cell arrays
                results{1} = linearTDSR(M,r,ui,alpha,opto,'phasic');
                f_m{1}(i,:) = [results{1}(f).V];
                d_m{1}(i,:) = [results{1}(d).V];
                c_m{1}(i,:) = [results{1}(c).V];
                
                results{2} = linearTDSR(M,r,ui,alpha);
                f_m{2}(i,:) = [results{2}(f).V];
                d_m{2}(i,:) = [results{2}(d).V];
                c_m{2}(i,:) = [results{2}(c).V];

            end
            
            V(1,:) = mean([f_m{1}(:,1) d_m{1}(:,1) c_m{1}(:,1)]);
            V(2,:) = mean([f_m{2}(:,1) d_m{2}(:,1) c_m{2}(:,1)]);
            
            figure;
            subplot(1,2,2);
            b = bar(V, 'FaceColor', 'flat');
            CData = lines(3);
            for c = 1:3
                b(c).FaceColor = CData(c,:);
            end   
            set(gca,'FontSize',25,'XLim', [.5 2.5],'XTickLabel',{'ChR2' 'eYFP'});
            ylabel('V','FontSize',25);
            title('Model','FontSize',25,'FontWeight','Bold')
            
            subplot(1,2,1);
            V = [2.05 1.66 2.66; 2.15 1.47 1.52];
            b = bar(V);
            for c = 1:3
                b(c).FaceColor = CData(c,:);
            end 
            legend({'F' 'D' 'C'},'FontSize',10, 'Location', 'North');
            set(gca,'FontSize',25,'XLim', [.5 2.5],'XTickLabel',{'ChR2' 'eYFP'});
            ylabel('CR','FontSize',25);
            title('Data','FontSize',25,'FontWeight','Bold')
            set(gcf,'Position',[200 200 1000 400])
            
            
        case 'sharpe17_deval'
            
            % A = 1
            % C = 2
            % D = 3
            % E = 4
            % F = 5
            % X = 6
            % food = 7
            
            A_X = [1 0 0 0 0 0 0; 1 0 0 0 0 1 0];    %24
            EF_X = [0 0 0 1 1 0 0; 0 0 0 1 1 1 0];   %8
            AD_X = [1 0 1 0 0 0 0; 1 0 1 0 0 1 0];   %8
            AC_X = [1 1 0 0 0 0 0; 1 1 0 0 0 1 0];   %8
            X = [0 0 0 0 0 1 0; 0 0 0 0 0 1 1];      %4*24
            C = [0 1 0 0 0 0 0];
            FOOD = [0 0 0 0 0 0 1];

            %This sets the alphas for the SR and U learning:
            alpha = [0.06 0.03];
            
            %this function provides the logicals of the first state of a stimulus
            st1 = @(x,Y) all(Y == x(1,:),2);
            
            %this function provides the second state of a stimulus
            st2 = @(x) [false; x(1:end - 1)];
            %st2 = @(x) all(M == x(2,:),2);
            
            %stimuli are randomized within each execution of the TDSR model. 
            for i = 1:100
                
                %this randomizes trial order for each stage
                s1 = randTrials(24,A_X);
                s2 = randTrials(8,EF_X,AD_X,AC_X);
                s3 = randTrials(4*24,X);
                s4 = repmat(FOOD,1,1);
                s5 = repmat(C,6,1);
                
                M = [s1;s2;s3;s4;s5];
                
                %This keeps track of trials from each stage
                sM = [ones(size(s1,1),1);2*ones(size(s2,1),1);3*ones(size(s3,1),1);4*ones(size(s4,1),1);5*ones(size(s5,1),1)];
                
                %logicals of the first state of each stimulus. The logicals
                %are conditioned on the stage in which the cue occured. (X
                %occurs in stage 1 and 2 as part of the other cues)
                a_x = st1(A_X,M) & sM == 1;
                ef_x = st1(EF_X,M) & sM == 2; ad_x = st1(AD_X,M) & sM == 2; ac_x = st1(AC_X,M) & sM == 2;
                x = st1(X,M) & sM == 3;
                food = st1(FOOD,M) & sM == 4;
                c = st1(C,M) & sM == 5;
                
                %rewarded states
                r1 = double(st2(x)) + double(food);
                
                %The -10 devalues the outcome
                r2 = double(st2(x)) + -10*double(food);
                
                opto = 1*(st2(ac_x) | st2(st2(ad_x)));
                
                %This sets the initial value of the US to 1
                ui = [0 0 0 0 0 0 0.5]';
                
                %this runs the model for the current iteration and saves
                %the results in the following cell arrays
                results{1} = linearTDSR(M,r1,ui,alpha,opto,'phasic');
                c_m{1}(i,:) = [results{1}(c).V];
                
                results{2} = linearTDSR(M,r2,ui,alpha,opto,'phasic');
                c_m{2}(i,:) = [results{2}(c).V];
                
            end

            V(1,:) = mean(c_m{1}(:,1));
            V(2,:) = mean(c_m{2}(:,1));
            
            figure;
            subplot(1,2,2);
            b = bar(V);%, 'FaceColor', 'flat');
            CData = lines(2);
            for c = 1:1
                b(c).FaceColor = CData(c,:);
            end 
            set(gca,'FontSize',25,'XTickLabel',{'Nondevalued' 'Devalued'},'XLim',[0.5 2.5]);
            ylabel('V','FontSize',25);
            title('Model','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,1)
            V = [3.22 1.23];
            b = bar(V);
            for c = 1:1
                b(c).FaceColor = CData(c,:);
            end 
            set(gca,'FontSize',25,'XTickLabel',{'Nondevalued' 'Devalued'},'XLim',[0.5 2.5]);
            ylabel('CR','FontSize',25);
            title('Data','FontSize',25,'FontWeight','Bold');
            set(gcf,'Position',[200 200 900 350])

            case 'sharpe17_deval_rehearsal'
            
            % A = 1
            % C = 2
            % D = 3
            % E = 4
            % F = 5
            % X = 6
            % food = 7
            
            A_X = [1 0 0 0 0 0 0; 1 0 0 0 0 1 0];    %24
            EF_X = [0 0 0 1 1 0 0; 0 0 0 1 1 1 0];   %8
            AD_X = [1 0 1 0 0 0 0; 1 0 1 0 0 1 0];   %8
            AC_X = [1 1 0 0 0 0 0; 1 1 0 0 0 1 0];   %8
            X = [0 0 0 0 0 1 0; 0 0 0 0 0 1 1];      %4*24
            C = [0 1 0 0 0 0 0];
            FOOD = [0 0 0 0 0 1 1];

            %This sets the alphas for the SR and U learning:
            alpha = [0.06 0.03];
            
            %this function provides the logicals of the first state of a stimulus
            st1 = @(x,Y) all(Y == x(1,:),2);
            
            %this function provides the second state of a stimulus
            st2 = @(x) [false; x(1:end - 1)];
            %st2 = @(x) all(M == x(2,:),2);
            
            %stimuli are randomized within each execution of the TDSR model. 
            for i = 1:100
                
                %this randomizes trial order for each stage. %Note that
                %randTrials adds a zeros vector between trials (models the
                %ITI)
                s1 = randTrials(24,A_X);
                s2 = randTrials(8,EF_X,AD_X,AC_X);
                s3 = randTrials(4*24,X);
                s4 = repmat(FOOD,1,1);
                s5 = repmat(C,6,1);
                
                M = [s1;s2;s3;s4;s5];
                
                %This keeps track of trials from each stage
                sM = [ones(size(s1,1),1);2*ones(size(s2,1),1);3*ones(size(s3,1),1);4*ones(size(s4,1),1);5*ones(size(s5,1),1)];
                
                %logicals of the first state of each stimulus. The logicals
                %are conditioned on the stage in which the cue occured. (X
                %occurs in stage 1 and 2 as part of the other cues)
                a_x = st1(A_X,M) & sM == 1;
                ef_x = st1(EF_X,M) & sM == 2; ad_x = st1(AD_X,M) & sM == 2; ac_x = st1(AC_X,M) & sM == 2;
                x = st1(X,M) & sM == 3;
                food = st1(FOOD,M) & sM == 4;
                c = st1(C,M) & sM == 5;
                
                %rewarded states
                r1 = double(st2(x)) + double(food);
                
                %-10 devalues the outcome
                r2 = double(st2(x)) + -10*double(food);
                
                opto = 1*(st2(ac_x) | st2(st2(ad_x)));
                
                %This sets the initial value of the US to 1
                ui = [0 0 0 0 0 0 0.5]';
                
                %this runs the model for the current iteration and saves
                %the results in the following cell arrays
                results{1} = linearTDSR(M,r1,ui,alpha,opto,'phasic');
                c_m{1}(i,:) = [results{1}(c).V];
                All{1}(:,:,i) = [results{1}.AllV];
                
                results{2} = linearTDSR(M,r2,ui,alpha,opto,'phasic');
                c_m{2}(i,:) = [results{2}(c).V];
                All{2}(:,:,i) = [results{2}.AllV];
                
            end
            
            V(1,:) = mean(c_m{1}(:,1));
            V(2,:) = mean(c_m{2}(:,1));
            
            figure;
            subplot(1,2,2);
            b = bar(V);%, 'FaceColor', 'flat');
            CData = lines(2);
            for c = 1:1
                b(c).FaceColor = CData(c,:);
            end 
            set(gca,'FontSize',25,'XTickLabel',{'Nondevalued' 'Devalued'},'XLim',[0.5 2.5]);
            ylabel('Value','FontSize',25);
            title('Model','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,1)
            V = [3.22 1.23];
            b = bar(V);
            for c = 1:1
                b(c).FaceColor = CData(c,:);
            end 
            set(gca,'FontSize',25,'XTickLabel',{'Nondevalued' 'Devalued'},'XLim',[0.5 2.5]);
            ylabel('Food cup entries','FontSize',25);
            title('Data','FontSize',25,'FontWeight','Bold');
            set(gcf,'Position',[200 200 900 350])

        case 'sharpe17_inhib'
            
            % A = 1
            % B = 2
            % X = 3
            % Y = 4
            % flavor1 = 5
            % flavor2 = 6
            
            A_X = [1 0 0 0 0 0; 1 0 1 0 0 0];   %12
            B_Y = [0 1 0 0 0 0; 0 1 0 1 0 0];   %12
            X   = [0 0 1 0 0 0; 0 0 1 0 1 0];   %4*24
            Y   = [0 0 0 1 0 0; 0 0 0 1 0 1];   %4*24
            A   = [1 0 0 0 0 0];
            B   = [0 1 0 0 0 0];
            
            %This sets the alphas for the SR and U learning:
            alpha = [0.06 0.03];
            
            %this function provides the logicals of the first state of a stimulus
            st1 = @(x,Y) all(Y == x(1,:),2);
            
            %this function provides the second state of a stimulus
            st2 = @(x) [false; x(1:end - 1)];
            %st2 = @(x) all(M == x(2,:),2);        
            
            %This gets the prior state
            stpre = @(x) [x(2:end); false];
            
            %stimuli are randomized within each execution of the TDSR model. 
            for i = 1:100
                
                %this randomizes trial order for each stage
                s1 = randTrials(12,A_X,B_Y);
                s2 = randTrials(24*4,X,Y);
                s3 = randTrials(6,A,B);
                
                M = [s1;s2;s3];

                %This keeps track of the stage
                sM = [ones(size(s1,1),1);2*ones(size(s2,1),1);3*ones(size(s3,1),1)];
                
                
                %logicals of each stimulus:
                b_y = st1(B_Y,M) & sM == 1;
                a = st1(A,M) & sM == 3;
                b = st1(B,M) & sM == 3;
                x = st1(X,M) & sM == 2;
                y = st1(Y,M) & sM == 2;
                
                
                %rewarded states
                r = double(st2(x | y));

                %tonic inhibition is on through the full cue
                opto = -0.8*(stpre(b_y) | b_y | st2(b_y));
                
                ui = [0 0 0 0 0.5 0.5]';
                
                results{1} = linearTDSR(M,r,ui,alpha,opto,'tonic');
                a_m{1}(i,:) = [results{1}(a).V];
                b_m{1}(i,:) = [results{1}(b).V];
                All{1}(:,:,i) = [results{1}.AllV];
                
                results{2} = linearTDSR(M,r,ui,alpha);
                a_m{2}(i,:) = [results{2}(a).V];
                b_m{2}(i,:) = [results{2}(b).V];
                All{2}(:,:,i) = [results{2}.AllV];
                
            end           
            
            V(1,:) = mean([a_m{1}(:,1) b_m{1}(:,1)]);
            V(2,:) = mean([a_m{2}(:,1) b_m{2}(:,1)]);
            
            figure;
            subplot(1,2,2);
            b = bar(V);
            CData = lines(2);
            for c = 1:2
                b(c).FaceColor = CData(c,:);
            end 
            set(gca,'FontSize',25,'XTickLabel',{'NpHR' 'eYFP'},'XLim', [.5 2.5]);
            ylabel('V','FontSize',25);
            title('Model','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,1);
            V = [51.4 29.3; 38.5 36.9];
            b = bar(V);
            for c = 1:2
                b(c).FaceColor = CData(c,:);
            end 
            legend({'A' 'B'},'FontSize',25);
            set(gca,'FontSize',25,'XTickLabel',{'NpHR' 'eYFP'},'XLim', [.5 2.5]);
            ylabel('CR','FontSize',25);
            title('Data','FontSize',25,'FontWeight','Bold');
            set(gcf,'Position',[200 200 1000 400])           
            
      
        case 'takahashi17_identity'
            
            % A (odor 1) = 1
            % B (odor 2) = 2
            % flavor 1   = 3
            % flavor 2   = 4
            %i: identity block shift
            
            n = 30;
            A  = [1 0 0 0; 1 0 1 0];
            B  = [0 1 0 0; 0 1 0 1];
            Ai = [1 0 0 0; 1 0 0 1];
            Bi = [0 1 0 0; 0 1 1 0];
            
            %This sets the alphas for the SR and U learning:
            alpha = [0.06 .03];
            
            %this function provides the logicals of the first state of a stimulus
            st1 = @(x,Y) all(diff([0 0 0 0;Y]) == x(1,:),2);
            
            %This runs identity then the upshift/downshift
            M = repmat([randTrials(n,A,B);randTrials(n,Ai,Bi)], 4,1);
            
            r = repmat([0; 1; 0],n*2*2*4,1);
            
            ai = st1(A,M);
            bi = st1(B,M); 
            
            ui = [0 0 .5 .5]';
            
            %First get the identity errors:
            results = linearTDSR(M,r,ui,alpha);
            
            %This gets the last two blocks
            ni = 0;
            for i=length(results) - n*2*4*3 + 1:length(results); ni = ni+1; dt(ni,:) = results(i).dt; end
            
            %use A to assess the identity error over the final two blocks:
            A_dt = dt(ai(length(results) - n*2*4*3 + 1:length(results)),:);
            
            %use 5 trials
            k = 5;
            
            D = [sum(mean(A_dt(1:k,:))); sum(mean(A_dt(n - k + 1:n,:)))];
            
            clear dt
            
            %This is the actual data:
            T_Data = [6.82 2.34]; %Identity

            figure;
            hold on
                subplot(1,2,1)
                bar(D);
                set(gca,'FontSize',25,'XTickLabel',{'Early' 'Late'},'XLim',[0.5 2.5],'YLim',[-.5 1]);
                ylabel('Prediction error','FontSize',25);
                title('Model','FontSize',25,'FontWeight','Bold');
                
                subplot(1,2,2)
                
                bar(T_Data);
                set(gca,'FontSize',25,'XTickLabel',{'Early' 'Late'},'XLim',[0.5 2.5],'YLim',[-5 15]);
                ylabel('Spikes/sec','FontSize',25);
                title('Data','FontSize',25,'FontWeight','Bold');
                set(gcf,'Position',[200 200 900 350])
        
      case 'chang17_identity'
            
            %randTrials adds the ITI between trials
            A =  [1 0 0 0 0 0; 1 0 0 0 1 0];
            B =  [0 1 0 0 0 0; 0 1 0 0 0 1];
            AX = [1 0 1 0 0 0; 1 0 1 0 1 0];
            BY = [0 1 0 1 0 0; 0 1 0 1 1 0];     
            X =  [0 0 1 0 0 0];
            Y =  [0 0 0 1 0 0];
            
            %This sets the alphas for the SR and U learning:
            alpha = [0.06 0.03];
            
            %this function provides the logicals of the first state of a stimulus
            st1 = @(x,Y) all(Y == x(1,:),2);
            
            %this function provides the second state of a stimulus
            st2 = @(x) [false; x(1:end - 1)];
            %st2 = @(x) all(M == x(2,:),2);        
            
            %This gets the prior state
            stpre = @(x) [x(2:end); false];
            
            %stimuli are randomized within each execution of the TDSR model. 
            for i = 1:10
                
                %this randomizes trial order for each stage
                s1 = randTrials(48,A, B);
                s2 = randTrials(16,AX,BY);
                s3 = randTrials(6,X,Y);
                
                M = [s1;s2;s3];
                
                %This keeps track of the stage
                sM = [ones(size(s1,1),1);2*ones(size(s2,1),1);3*ones(size(s3,1),1)];
                
                %logicals of the first state of each stimulus
                a = st1(A,M) & sM == 1;
                b = st1(B,M) & sM == 1;
                ax = st1(AX,M) & sM == 2;
                by = st1(BY,M) & sM == 2;
                x = st1(X,M) & sM == 3;
                y = st1(Y,M) & sM == 3;
                                
                %rewarded states
                r = double(st2(a) | st2(b) | st2(ax) | st2(by));
                
                %This implements tonic inhibition which occurs through the
                %full trial: ITI, CS, CS + US
                opto{1} = -.8*(stpre(by) | by | st2(by));
                opto{2} = -.8*([by(2:end);0]);
                
                ui = [0 0 0 0 .5 .5]';
                
                results{1} = linearTDSR(M,r,ui,alpha,opto{1},'tonic');
                x_m{1}(i,:) = [results{1}(x).V];
                y_m{1}(i,:) = [results{1}(y).V];
                All{1}(:,:,i) = [results{1}.AllV];
                
                results{2} = linearTDSR(M,r,ui,alpha,opto{2},'tonic');
                x_m{2}(i,:) = [results{2}(x).V];
                y_m{2}(i,:) = [results{2}(y).V];
                All{2}(:,:,i) = [results{2}.AllV];
                
            end
            
            V(1,:) = mean([x_m{1}(:,1) y_m{1}(:,1)]);
            V(2,:) = mean([x_m{2}(:,1) y_m{2}(:,1)]);
            
            figure
            subplot(1,2,2);
            b = bar(V);
            CData = lines(2);
            for c = 1:2
                b(c).FaceColor = CData(c,:);
            end 
            set(gca,'FontSize',25,'XLim', [.5 2.5],'XTickLabel',{'Exp' 'ITI'});
            ylabel('V','FontSize',25);
            title('Model','FontSize',25,'FontWeight','Bold');
            
            subplot(1,2,1);
            V = [27.6 17.1; 26.7 47.4];
            b = bar(V);
            for c = 1:2
                b(c).FaceColor = CData(c,:);
            end 
            legend({'A_B' 'A_{UB}'},'FontSize',15,'Location','North');
            set(gca,'FontSize',25,'XLim', [.5 2.5],'XTickLabel',{'Exp' 'ITI'});
            ylabel('CR','FontSize',25);
            title('Data','FontSize',25,'FontWeight','Bold');
            set(gcf,'Position',[200 200 900 350])      
            
    end