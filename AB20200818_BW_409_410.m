%% Script to analyse plate-reader data from 18/08/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strain: BW25113
% Constructs: BW25113 + AB409
%             BW25113 + AB410
% Media: M9 0.8% fructose with Casaminoacids
% Level of ara & rha induction: 0%, 0.0002%, 0.002%, 0.005%, 0.01%, 0.02%,...
%                               0.05%, 0.1%, 0.2%, 0.5%, 1%, 2%
% GFP Gain: 50
% Plate-reader: Tecan Spark
% Temperature: 37�C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define indexes for each construct (in the rates matrixes)

idx_m9 = [1, 2, 3];
idx_bw = [4, 5, 6];

AB409_0 = [4, 5, 6];
AB409_1 = [1+12*1, 1+12*2, 1+12*3]; % AB409 + 0% ara
AB409_2 = [2+12*1, 2+12*2, 2+12*3]; % AB409 +0.0002% ara
AB409_3 = [3+12*1, 3+12*2, 3+12*3]; % AB409 + 0.002% ara
AB409_4 = [4+12*1, 4+12*2, 4+12*3]; % AB409 + 0.005% ara
AB409_5 = [5+12*1, 5+12*2, 5+12*3]; % AB409 + 0.01% ara
AB409_6 = [6+12*1, 6+12*2, 6+12*3]; % AB409 + 0.02% ara
AB409_7 = [7+12*1, 7+12*2, 7+12*3]; % AB409 + 0.05% ara
AB409_8 = [8+12*1, 8+12*2, 8+12*3]; % AB409 + 0.1% ara
AB409_9 = [9+12*1, 9+12*2, 9+12*3]; % AB409 + 0.2% ara
AB409_10 = [10+12*1, 10+12*2, 10+12*3]; % AB409 + 0.5% ara
AB409_11 = [11+12*1, 11+12*2, 11+12*3]; % AB409 + 1% ara
AB409_12 = [12+12*1, 12+12*2, 12+12*3]; % AB409 + 2% ara

AB410_0 = [4, 5, 6];
AB410_1 = [1+12*4, 1+12*5, 1+12*6]; % AB410 + 0% rha
AB410_2 = [2+12*4, 2+12*5, 2+12*6]; % AB410 + 0.0002% rha
AB410_3 = [3+12*4, 3+12*5, 3+12*6]; % AB410 + 0.002% rha
AB410_4 = [4+12*4, 4+12*5, 4+12*6]; % AB410 + 0.005% rha
AB410_5 = [5+12*4, 5+12*5, 5+12*6]; % AB410 + 0.01% rha
AB410_6 = [6+12*4, 6+12*5, 6+12*6]; % AB410 + 0.02% rha
AB410_7 = [7+12*4, 7+12*5, 7+12*6]; % AB410 + 0.05% rha
AB410_8 = [8+12*4, 8+12*5, 8+12*6]; % AB410 + 0.1% rha
AB410_9 = [9+12*4, 9+12*5, 9+12*6]; % AB410 + 0.2% rha
AB410_10 = [10+12*4, 10+12*5, 10+12*6]; % AB410 + 0.5% rha
AB410_11 = [11+12*4, 11+12*5, 11+12*6]; % AB410 + 1% rha
AB410_12 = [12+12*4, 12+12*5, 12+12*6]; % AB410 + 2% rha

% End of index definition section

%% Import Data

load('AB20200818_BW_409_410.mat');

% Import the OD data and substract the blank
OD = OD700 - mean(OD700(:,idx_m9),2);

% Substract GFP background
for j = 1:2
    for i = 1:13
        if j == 1
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
            GFP(:,idx) = GFP(:,idx) - mean(GFP(:,AB409_0),2);
        elseif j == 2
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
            GFP(:,idx) = GFP(:,idx) - mean(GFP(:,AB410_0),2); %#ok<*SAGROW>
        end
    end
end

% Add time to OD and GFP data
timeLength = size(OD700,1);
t = 0;
inc = 0;
delaT = 15/60;

for n = 2:timeLength
    t = [t inc+delaT]; %#ok<AGROW>
    inc = inc + delaT;
end

OD = horzcat(t',OD(1:length(t'),:));
GFP = horzcat(t',GFP(1:length(t'),:));

% Replace NaNs with 0s
OD(isnan(OD)) = 0;
GFP(isnan(GFP)) = 0;

% end of Import data section

%% Calculate the growth rate, and the CFP production rate per cell

timeLength = length(delaT:delaT:t(end-1));
timeIDX = 1;
timeIdx = 0;

% growth_rate = zeros(1,96);
growth_rate = zeros(timeLength,96);
GFP_rate = zeros(timeLength,96);
G_rate = cell(2,13);
MAXG = cell(2,13);
gfp_rate = cell(2,13);

%% Calculate Smoothed Growth Rate

for n = 1:96
    
    f = fit(t',OD(:,n+1),'smoothingspline','SmoothingParam',0.8648426188005848);
    fOD = f(t);
    growth_rate(:,n) = (log(fOD(3:end))-log(fOD(1:end-2)))/0.5;
    
    f1 = fit(t',GFP(:,n+1),'smoothingspline','SmoothingParam',0.8648426188005848);
    fGFP = f1(t);
    GFP_rate(:,n) = (fGFP(3:end)-fGFP(1:end-2))./fOD(2:end-1)/0.5;
  
end

for j = 1:2
    
    for i = 1:13
       
        if j == 1
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
        else
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
        end
        
        G_rate{j,i} = growth_rate(:,idx); % record growth rate
        gfp_rate{j,i} = GFP_rate(:,idx); % record GFP production rate per cell
        
        [~ , maxGIdx] = max(mean(growth_rate(:,idx),2));
        MAXG{j,i} = growth_rate(maxGIdx,idx);
    end
end
    
%% Plot Bar Graphs
    
for time = [1 2]
    
    timeIdx = find(round(t,2) == time);
    
%%%%%%%%%%%%%%%%%%%%%%%% Plot Growth Rate Figure %%%%%%%%%%%%%%%%%%%%%%%% 

    figure;
    names = {'pBAD-GFP','pRha-GFP'};
    legend_barweb = {'0% ','0.0002% ','0.002% ','0.005% ','0.01% ','0.02% ','0.05% ','0.1% ','0.2% ','0.5% ','1% ','2% '};
    subplot(2,1,1);
    
    growth = cell2mat(cellfun(@(x) mean(x(timeIdx,:),2),G_rate(:,:),'UniformOutput',0));
    growth_std = cell2mat(cellfun(@(x) mean(std(x(timeIdx,:),0,2),2),G_rate(:,:),'UniformOutput',0));
    
    barweb(growth(:,2:end), growth_std(:,2:end), 0.8, names, ['Growth Rate ' int2str(time) 'h Post-Induction'], [], ...
        'Growth Rate [h^{-1}]', [0.984, 0.639, 0.274], 'xy', legend_barweb, 2,'axis');

    set(gca,'FontSize',16,'FontWeight','bold','linewidth',1.5); ylim([0 1.2]);
    h = get(gca,'position'); set(gca,'position',[h(1), h(2)+0.02, 0.6, 0.2]);
    
    %%%%%%%%%%%%%%%%%%%%%%%% Plot GFP Rate Figure %%%%%%%%%%%%%%%%%%%%%%%% 

    subplot(2,1,2);
    
    gfp = cell2mat(cellfun(@(x) mean(x(timeIdx,:),2),gfp_rate(:,:),'UniformOutput',0));
    GFP_std = cell2mat(cellfun(@(x) mean(std(x(timeIdx,:),0,2),2),gfp_rate(:,:),'UniformOutput',0));
    
    barweb(gfp(:,2:end), GFP_std(:,2:end), 0.8, names, ['GFP Production Rate per Cell ' int2str(time) 'h Post-Induction'], [], ...
        {'GFP Production Rate','per Cell [F.U.ABS_{700}^{-1}h^{-1}]'}, [0.435, 0.980, 0.380], 'xy', legend_barweb, 2,'axis');

    set(gca,'FontSize',16,'FontWeight','bold','linewidth',1.5); ylim([0 6e4]);
    h = get(gca,'position'); set(gca,'position',[h(1), h(2)+0.1, 0.6, 0.2]);

end

%% Plot production rates

colour = jet(13);
% colour(9,:) = [1, 0.937, 0.141];
tEND = t(end);
    
LEG ={'WT ','0% ','0.0002% ','0.002% ','0.005% ','0.01% ','0.02% ','0.05% ','0.1% ','0.2% ','0.5% ','1% ','2% '};

%% OD & Growth Rate Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OD Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(11);

for j = 1:2
    
    subplot(2,2,j);
    
    for i = 1:13
  
        if j == 1
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
            Title = {'OD of pBAD-GFP'};
        elseif j == 2
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
            Title = {'OD of pRha-GFP'};
        end

%         p(i) = stdshade(OD(:,idx+1)',0.2,colour(i,:),OD(:,1),1); %#ok<*SAGROW>
%         hold on
        
        for m = 1:3
            f = fit(t',OD(:,idx(m)+1),'smoothingspline','SmoothingParam',0.8648426188005848);
            fOD(:,m) = f(t);
        end
        
        p(i) = stdshade(fOD',0.1,colour(i,:),OD(:,1),1);
        hold on
    end
    
    grid on; axis square; xlim([delaT tEND-1]); ylim([0 1]);
    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2,'TickLength',[0.02,0.02]);
    xlabel('Time [h]'); legend(p(:),LEG(:),'location','NorthEastOutside');
    ylabel({'ABS_{700}'}); title(Title);hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Growth Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11);

for j = 1:2
  
    subplot(2,2,j+2);
    
    for i = 1:13
        
        if j == 1
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
            Title = {'Growth Rate of pBAD-GFP'};
        elseif j == 2
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
            Title = {'Growth Rate of pRha-GFP'};
        end
        
        for m = 1:3
            f = fit(t',OD(:,idx(m)+1),'smoothingspline','SmoothingParam',0.8648426188005848);
            fOD(:,m) = f(t);
            fG(:,m) = (log(fOD(3:end,m))-log(fOD(1:end-2,m)))/0.5;
        end
        
        p(i) = stdshade(fG',0.1,colour(i,:),OD(2:end-1,1),1);
        hold on
        
%         p(i) = stdshade(G_rate{j,i}',0.2,colour(i,:),OD(2:end-1,1),1);
%         hold on
    end
    
    grid on; axis square; xlim([delaT tEND-1]); ylim([-0.2 1.5]);
    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2,'TickLength',[0.02,0.02]);
    xlabel('Time [h]'); legend(p(:),LEG(:),'location','NorthEastOutside');
    ylabel({'Growth Rate [h^{-1}]'}); title(Title);hold off;
end

%% GFP/cell & GFP Capacity Figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GFP/Cell Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(12);

for j = 1:2
    
    subplot(2,2,j);
    
    for i = 1:13
        
        if j == 1
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
            Title = {'Total GFP of pBAD-GFP'};
        elseif j == 2
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
            Title = {'Total GFP of pRha-GFP'};
        end
        
        for m = 1:3
            f = fit(t',GFP(:,idx(m)+1),'smoothingspline','SmoothingParam',0.8648426188005848);
            fGFP(:,m) = f(t);
        end
        
        p(i) = stdshade(fGFP',0.1,colour(i,:),OD(:,1),1);

        %p(i) = stdshade((GFP(:,idx+1))',0.2,colour(i,:),OD(:,1),1); %#ok<*SAGROW>
        hold on
    end
    
    grid on; axis square; xlim([delaT tEND-1]); ylim([-1e4 5e4]);
    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2,'TickLength',[0.02,0.02]);
    xlabel('Time [h]'); legend(p(:),LEG(:),'location','NorthEastOutside');
    ylabel({'F.U.'}); title(Title);hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%% GFP Production Rate per Cell Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:2
  
    subplot(2,2,j+2);
    
    for i = 1:13
        
        if j == 1
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
            Title = {'GFP Production Rate', 'per Cell of pBAD-GFP'};
        elseif j == 2
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
            Title = {'GFP Production Rate', 'per Cell of pRha-GFP'};
        end
        
        for m = 1:3
            f1 = fit(t',OD(:,idx(m)+1),'smoothingspline','SmoothingParam',0.8648426188005848);
            fOD(:,m) = f1(t);
            f2 = fit(t',GFP(:,idx(m)+1),'smoothingspline','SmoothingParam',0.8648426188005848);
            fGFP(:,m) = f2(t);
            frGFP(:,m) = ((fGFP(3:end,m)-fGFP(1:end-2,m))./fOD(1:end-2,m))/0.5;
        end
        
        p(i) = stdshade(frGFP',0.1,colour(i,:),OD(2:end-1,1),1);

%         p(i) = stdshade(GFP_rate{j,i}',0.2,colour(i,:),OD(2:end-1,1),1);
        hold on
    end
    
    grid on; axis square; xlim([delaT tEND-1]); ylim([-1e4 7e4]);
    set(gca,'FontSize',14,'FontWeight','bold','linewidth',2,'TickLength',[0.02,0.02]);
    xlabel('Time [h]'); legend(p(1:8),LEG(1:8),'location','NorthEastOutside');
    ylabel({'GFP Production Rate per Cell', '[F.U.ABS_{700}^{-1}h^{-1}]'}); title(Title);hold off;
end


%% Static Input-Output Curves

colour = [1 0 0; 0 0 1];
timePoint = find(t==3);
concentrations = [0,0.0002,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2];

figure
hAx = axes;
hAx.XScale = 'log';
hold all

for j = 1:2

    for i = 2:13
        
        if j == 1 
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
        else
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
        end

        TF{j}(i-1,:) = [concentrations(i-1),gfp_rate{j,i}(timePoint,:)];
    end
    
    f = @(F,x) F(1)./(F(2)+exp(-F(3).*x));
    F_fitted = nlinfit(TF{j}(:,1),mean(TF{j}(:,2:end),2),f,[0 0 0]);
    disp(['F = ',num2str(F_fitted)])
    
    errorbar(TF{j}(:,1),mean(TF{j}(:,2:end),2),std(TF{j}(:,2:end),0,2),'o','linewidth',2,'Color',colour(j,:));
    hold on;
    p(j) = plot(linspace(0.0001,3,1e6),f(F_fitted,linspace(0.0001,3,1e6)),'linewidth',2,'Color',colour(j,:));
end

title({'Static Input-Output Curve', ['at ' num2str(t(timePoint)) 'h post-induction']});
xlabel('3OC6-HSL/pC-HSL Concentration [M]');ylabel({'GFP Production Rate per Cell','[F.U.ABS_{700}^{-1}h^{-1}]'});
set(gca,'FontSize',16,'FontWeight','bold','linewidth',2,'TickLength',[0.02,0.02]); ylim([-5e3 4e4]);xlim([0.0001 3]);
set(gcf,'position',[10,10,550,350]); legend(p(1:2),{'pBAD-GFP','pRha-GFP'},'Location','NorthWest');

%% Static Input-Output Curves Total GFP

colour = [1 0 0; 0 0 1];
timePoint = find(t==6);
concentrations = [0,0.0002,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2];

figure
hAx = axes;
hAx.XScale = 'log';
hold all

for j = 1:2

    for i = 2:13
        
        if j == 1 
            idx =  eval(matlab.lang.makeValidName(strcat('AB409_',int2str(i-1))));
        else
            idx =  eval(matlab.lang.makeValidName(strcat('AB410_',int2str(i-1))));
        end

        TF{j}(i-1,:) = [concentrations(i-1),GFP(timePoint,idx)];
    end
    
    f = @(F,x) F(1)./(F(2)+exp(-F(3).*x)); %F(1) + (F(2)./(1 + exp(F(3).*x+F(4))));
    F_fitted = nlinfit(TF{j}(:,1),mean(TF{j}(:,2:end),2),f,[0 0 0 0]);
    disp(['F = ',num2str(F_fitted)])
    
    errorbar(TF{j}(:,1),mean(TF{j}(:,2:end),2),std(TF{j}(:,2:end),0,2),'o','linewidth',2,'Color',colour(j,:));
    hold on;
    p(j) = plot(linspace(0.0001,3,1e7),f(F_fitted,linspace(0.0001,3,1e7)),'linewidth',2,'Color',colour(j,:));
end

title({'Static Input-Output Curve', ['at ' num2str(t(timePoint)) 'h post-induction']});
xlabel('3OC6-HSL/pC-HSL Concentration [M]');ylabel({'Total End-Point GFP','[F.U.]'});
set(gca,'FontSize',16,'FontWeight','bold','linewidth',2,'TickLength',[0.02,0.02]); ylim([-5e3 3e4]);xlim([0.0001 3]);
set(gcf,'position',[10,10,550,350]); legend(p(1:2),{'pBAD-GFP','pRha-GFP'},'Location','NorthWest');


%% Maximum Growth Rates Figure

for j = 1:2 
    
    for i = 1:15
       
        growthRate = G_rate{j,i};
        [~ , maxGIdx] = max(mean(growthRate(2:end,:),2));
        MAXG{j,i} = growthRate(maxGIdx,:);
        
    end
end

figure;
names = {'pBAD-GFP','pRha-GFP'};
legend_barweb = {'WT ','0M ','10^{-14}M ','10^{-15}M ','10^{-12}M ','10^{-11}M ',...
        '10^{-10}M ','10^{-9}M ','2.5*10^{-9}M ','5*10^{-9}M ','7.5*10^{-9}M ',...
        '10^{-8}M ','10^{-7}M ','10^{-6}M ','10^{-5}M ','10^{-4}M '};
subplot(2,1,1);

growth = cell2mat(cellfun(@(x) mean(x(1,:),2),MAXG(:,:),'UniformOutput',0));
growth_std = cell2mat(cellfun(@(x) mean(std(x(1,:),0,2),2),MAXG(:,:),'UniformOutput',0));

barweb(growth(:,2:end), growth_std(:,2:end), 0.8, names, 'Maximum Growth Rate', [], ...
    'Growth Rate [h^{-1}]', [0.992, 0.894, 0.305], 'xy', legend_barweb, 2,'axis');

set(gca,'FontSize',16,'FontWeight','bold','linewidth',1.5); ylim([0 1.2]);
h = get(gca,'position'); set(gca,'position',[h(1), h(2)+0.07, 0.6, 0.2]);

