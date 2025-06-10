clear; clc;
close all
clear all
% Random seed
% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));
% addpath('/Volumes/MyBook2/zebrafish vs drosophila/ModelCompare')
addpath('./Function')
addpath('./ExpData')
%datapath =  'front_041223/';
datapath =  'main_code_linlin_P_front/front_100422_backup/';
% datapath = 'front_041223/';

znames = {'k2','decB','decS','decBS','DB','DS','DBS','lambdaS',...
           'lambdaBS','j1','j2'};
f1 = figure;
hold on
f2 = figure;
hold on
f3 = figure;
hold on
f4 = figure;
hold on


% %figure for pramters
% f1_1 = figure;
% hold on
% f1_2 = figure;
% hold on
% f1_3 = figure;
% hold on
% f1_4 = figure;
% hold on
% f1_5 = figure;
% hold on
% f1_6 = figure;
% hold on
% f1_7 = figure;
% hold on
% f1_8 = figure;
% hold on
% f1_9 = figure;
% hold on
% f1_10 = figure;
% hold on
% f1_11 = figure;
% hold on



plot_case = 1; % check if plot anchor cases ==1(plot) else(no_plot)

% cmap = colormap(jet(d+1));
% cmap = cmap(randperm(length(cmap)),:);

for m = 0%:3 %:1%:3 %'0 == Origianl screen'; '1 == 1D cases'; '2 == 2D cases'; '3 == 3D cases'
% datapath = './Result_120321_simple/SSE_scenario_Simple_';
    
    if m == 0
         
       folders_names{1} = '';
       folder_n = 1;
       plot_name = 'Degree0';
       plot_color = 'k';
                                   
    elseif m == 1
        
   
       folder_n = 11;
       for i = 1:size(znames,2)
           folders_names{i} = [znames{i}];
       end
       plot_name = 'Degree1';
        plot_color = 'b';
       
    elseif m == 2
        
       folder_n = 11*10/2
       p_count=1;
      for i = 1:size(znames,2)
        for j = (i+1):size(znames,2)
            if i~=j
                folders_names{p_count} = [znames{i} znames{j}];
                p_count = p_count+1;
            end
        end
      end
       plot_name = 'Degree2';
        plot_color = 'g';
       
    elseif m == 3
        
       folder_n = 11*10*9/3/2
       p_count=1;
       for i = 1:size(znames,2)
            for j = (i+1):size(znames,2)
                for k = (j+1):size(znames,2)
                    if (k~=i) && (k~=j)
                        folders_names{p_count} = [znames{i} znames{j} znames{k}];
                        p_count = p_count+1;
      
                    end
                end
            end
       end
       plot_name = 'Degree3';
       plot_color = 'r';
    end
    
    total_ia = [];
    total_ia1 = [];
    total_ia2 = [];
    
    total_ndx = [];
    total_ndx1 = [];
    total_ndx2 = [];
    
    total_ndy = [];
    total_ndy1 = [];
    total_ndy2 = [];

   total_nond_drosS = [];
   total_nond_drosWT = [];
   total_nond_zebraWT = [];
   total_nond_zebraC = [];

   total_nond1_drosS = [];
   total_nond1_drosWT = [];
   total_nond1_zebraWT = [];
   total_nond1_zebraC = [];

   total_nond2_drosS = [];
   total_nond2_drosWT = [];
   total_nond2_zebraWT = [];
   total_nond2_zebraC = [];
    
    for f = 1:folder_n
        
        
       
       dataname =  ['Front_' folders_names{f}];
       
       load(strcat(datapath,dataname,'.mat'))
       
%        if plot_case ==1
       % get the solution for the WT front
%        drosWT = ones(length(ia),72);
%        zebraWT = ones(length(ia),72);
%        drosS = ones(length(ia),72);
%        zebraC = ones(length(ia),72);
       
%        for i0 = 1:length(ia)
%            ia_file = ceil(ia(i0)/2000);
%            ia_case = rem(ia(i0),2000);
%            filename1 = ['results/Result_Simple_' folders_names{f} '/SSE_scenario_Simple_' folders_names{f} num2str(ia_file) '.mat'];
%            if isfile(filename1)
%                load(filename1);
%                drosWT(i0,:)= Stored_Soln.dwt(ia_case,:);
%                drosS(i0,:)= Stored_Soln.dshomo(ia_case,:);
%                zebraWT(i0,:)= Stored_Soln.zwt57(ia_case,:);
%                zebraC(i0,:)= Stored_Soln.zclf57(ia_case,:);
% %            else
% %                continue
%            end
%            
%        end
       
       % get the solution for the CLF front
       
%        drosWT1 = ones(length(ia1),72);
%        zebraWT1 = ones(length(ia1),72);
%        drosS1 = ones(length(ia1),72);
%        zebraC1 = ones(length(ia1),72);
%        
%        for i1 = 1:length(ia1)
%            ia_file = ceil(ia1(i1)/2000);
%            ia_case = rem(ia1(i1),2000);
%            filename1 = ['results/Result_Simple_' folders_names{f} '/SSE_scenario_Simple_' folders_names{f} num2str(ia_file) '.mat'];
%            if isfile(filename1)
%                load(filename1);
%                drosWT1(i1,:)= Stored_Soln.dwt(ia_case,:);
%                drosS1(i1,:)= Stored_Soln.dshomo(ia_case,:);
%                zebraWT1(i1,:)= Stored_Soln.zwt57(ia_case,:);
%                zebraC1(i1,:)= Stored_Soln.zclf57(ia_case,:);
% %            else
% %                continue
%            end
%            
%        end
%        
%        % get the solution for the CLF front
%        drosWT2 = ones(length(ia2),72);
%        zebraWT2 = ones(length(ia2),72);
%        drosS2 = ones(length(ia2),72);
%        zebraC2 = ones(length(ia2),72);
%        
%        for i2 = 1:length(ia2)
%            ia_file = ceil(ia2(i2)/2000);
%            ia_case = rem(ia2(i2),2000);
%            filename1 = ['results/Result_Simple_' folders_names{f} '/SSE_scenario_Simple_' folders_names{f} num2str(ia_file) '.mat'];
%            if isfile(filename1)
%                load(filename1);
%                drosWT2(i2,:)= Stored_Soln.dwt(ia_case,:);
%                drosS2(i2,:)= Stored_Soln.dshomo(ia_case,:);
%                zebraWT2(i2,:)= Stored_Soln.zwt57(ia_case,:);
%                zebraC2(i2,:)= Stored_Soln.zclf57(ia_case,:);
% %            else
% %                continue
%            end
%        end
%        end
       
       total_ndx = [total_ndx;ndx];
       total_ndx1 = [total_ndx1;ndx1];
       total_ndx2 = [total_ndx2;ndx2];
       
       total_ndy = [total_ndy;ndy];
       total_ndy1 = [total_ndy1;ndy1];
       total_ndy2 = [total_ndy2;ndy2];  
       
       total_nond_drosS = [total_nond_drosS;nond_drosS];
       total_nond_drosWT = [total_nond_drosWT;nond_drosWT];
       total_nond_zebraWT = [total_nond_zebraWT;nond_zebraWT];
       total_nond_zebraC = [total_nond_zebraC;nond_zebraC];
       
       total_nond1_drosS = [total_nond1_drosS;nond1_drosS];
       total_nond1_drosWT = [total_nond1_drosWT;nond1_drosWT];
       total_nond1_zebraWT = [total_nond1_zebraWT;nond1_zebraWT];
       total_nond1_zebraC = [total_nond1_zebraC;nond1_zebraC];
       
       total_nond2_drosS = [total_nond2_drosS;nond2_drosS];
       total_nond2_drosWT = [total_nond2_drosWT;nond2_drosWT];
       total_nond2_zebraWT = [total_nond2_zebraWT;nond2_zebraWT];
       total_nond2_zebraC = [total_nond2_zebraC;nond2_zebraC];
       
    end
    
    if folder_n>1
       [fornt_x_wt, fornt_y_wt, iia] = find_nondominated_lim(total_ndx , total_ndy,0.2);
       [fornt_x_chd, fornt_y_chd, iia1] = find_nondominated_lim(total_ndx1 , total_ndy1,0.2);
       [fornt_x_wtchd, fornt_y_wtchd, iia2] = find_nondominated_lim(total_ndx2 , total_ndy2,0.2);
       
       % leave the data only at front
        total_nond_drosS = total_nond_drosS(iia,:);
        total_nond_drosWT = total_nond_drosWT(iia,:);
        total_nond_zebraWT = total_nond_zebraWT(iia,:);
        total_nond_zebraC = total_nond_zebraC(iia,:);
        
        total_nond1_drosS = total_nond1_drosS(iia1,:);
        total_nond1_drosWT = total_nond1_drosWT(iia1,:);
        total_nond1_zebraWT = total_nond1_zebraWT(iia1,:);
        total_nond1_zebraC = total_nond1_zebraC(iia1,:);
        
        total_nond2_drosS = total_nond2_drosS(iia2,:);
        total_nond2_drosWT = total_nond2_drosWT(iia2,:);
        total_nond2_zebraWT = total_nond2_zebraWT(iia2,:);
        total_nond2_zebraC = total_nond2_zebraC(iia2,:);
       
        
    else %for the initial sample only have one folder
        fornt_x_wt = total_ndx;
        fornt_y_wt = total_ndy;
        fornt_x_chd = total_ndx1;
        fornt_y_chd = total_ndy1;
        fornt_x_wtchd = total_ndx2;
        fornt_y_wtchd = total_ndy2;
        
        iia = ia;
        iia1= ia1;
        iia2 = ia2;
        
    end
    
    figure(f1)
    lineH1 = plot(fornt_x_wt, fornt_y_wt, '-o','Color',plot_color,'LineWidth',2.5,'DisplayName',plot_name);
    set(gca,'FontSize',15)
    title('WT front')
    xlim([0,0.2])
    ylim([0,0.2])
    xlabel('Drosophila RMSE')
    ylabel('Zebrafish RMSE')
    legend()
    set(gcf,'color','w');
    color1 = get(lineH1, 'Color');
    
    figure(f2)
    lineH2 = plot(fornt_x_chd, fornt_y_chd, '-o','Color',plot_color,'LineWidth',2.5,'DisplayName',plot_name);
    set(gca,'FontSize',15)
    title('Mutant front')
    xlim([0,0.2])
    ylim([0,0.2])
    xlabel('Drosophila RMSE')
    ylabel('Zebrafish RMSE')
    legend()
    set(gcf,'color','w');
    color2 = get(lineH2, 'Color');
    
    figure(f3)
    lineH3 = plot(fornt_x_wtchd*2, fornt_y_wtchd*2, '-o','Color',plot_color,'LineWidth',2.5,'DisplayName',plot_name);
    set(gca,'FontSize',15)
    title('Mutant+WT front')
    xlim([0,0.4])
    ylim([0,0.4])
    xlabel('Drosophila RMSE')
    ylabel('Zebrafish RMSE')
    legend()
    set(gcf,'color','w');
    color3 = get(lineH3, 'Color');
    

    
    figure(f4)
    lineH1 = plot(fornt_x_wt, fornt_y_wt, '-*','Color',plot_color,'LineWidth',2.5,'DisplayName',['WT' plot_name]);
    set(gca,'FontSize',15)
    title('WT front')
    xlim([0,0.2])
    ylim([0,0.2])
    xlabel('Drosophila RMSE')
    ylabel('Zebrafish RMSE')
    legend()
    set(gcf,'color','w');
    color1 = get(lineH1, 'Color');
    lineH3 = plot(fornt_x_wtchd*2, fornt_y_wtchd*2, '--o','Color',plot_color,'LineWidth',2.5,'DisplayName',['WT+mutant' plot_name]);
    set(gca,'FontSize',15)
    title('Mutant+WT front')
    xlim([0,0.4])
    ylim([0,0.4])
    xlabel('Drosophila RMSE')
    ylabel('Zebrafish RMSE')
    legend()
    set(gcf,'color','w');
    color3 = get(lineH3, 'Color');
    

    
    if plot_case ==1 
    load('experData'); % experimental data
    % WT anchor points-----------------------------------------------------
    nond_drosWT= total_nond_drosWT(:,:);%./repmat(max(total_nond_drosS(:,:),[],2),1,72);  
    nond_zebraWT = total_nond_zebraWT(:,:);%./repmat(max(total_nond_zebraWT(:,:),[],2),1,72);
    nond_drosS= total_nond_drosS(:,:);%./repmat(2*max(total_nond_drosS(:,:),[],2),1,72);  
    nond_zebraC= total_nond_zebraC(:,:);%./repmat(max(total_nond_zebraC(:,:),[],2),1,72);
    
    nodom = [fornt_x_wt,fornt_y_wt];
    cAnch = [nodom(1,:);nodom(end,:)]; 
    pUtopia = min(cAnch);
         
    dist = sqrt((fornt_x_wt-pUtopia(1)).^2 + (fornt_y_wt-pUtopia(2)).^2);
    [~, ind] = sort(dist);
    figure(f1)
    scatter(nodom(ind(1),1),nodom(ind(1),2),50,color1,"^",'DisplayName',['UP' plot_name])
    
%     % get_utopia parameter
%     Front_parameter_WTMU= Par(ia,:);
%     Utopia_parameter_WTMU = Front_parameter_WTMU(ind(1),:);
    
    selectd = [3,4;1,2;3,4;1,2;3,4;1,2;];
    ysim = [(nond_drosWT(1,1:36)), (nond_zebraWT(1,1:36)),...
        (nond_drosWT(ind(1),1:36)), (nond_zebraWT(ind(1),1:36)),...
        (nond_drosWT(end,1:36)), (nond_zebraWT(end,1:36));...
        (nond_drosS(1,1:36)), (nond_zebraC(1,1:36)),...
        (nond_drosS(ind(1),1:36)), (nond_zebraC(ind(1),1:36)),...
        (nond_drosS(end,1:36)), (nond_zebraC(end,1:36))];

    % selectd = [3,1]';
    % ysim = [(nond_drosWT(ind(1),1:36))', (nond_zebraWT(ind(1),1:36))'];
     xsim = [linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36);...
        linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36)];
    plotParetoFigures(experData, xsim,ysim,selectd)

    [~,h1]=suplabel('DV axis', 'x'); 
    set(h1,'FontSize',16) 
    [~,h2]=suplabel('Normalized BMP intensity','y'); 
    set(h2,'FontSize',16)
    [~,h3]=suplabel([plot_name ' WT fitting'] ,'t');  
    set(h3,'FontSize',20)
    set(gcf, 'Position',  [100, 100, 700, 800])
    set(gcf,'color','w');
    
    
    % CLF anchor points----------------------------------------------------
    nond1_drosWT= total_nond1_drosWT(:,:);%./repmat(max(total_nond1_drosWT(:,:),[],2),1,72);  
    nond1_zebraWT = total_nond1_zebraWT(:,:);%./repmat(max(total_nond1_zebraWT(:,:),[],2),1,72);
    nond1_drosS= total_nond1_drosS(:,:);%./repmat(2*max(total_nond1_drosS(:,:),[],2),1,72);  
    nond1_zebraC= total_nond1_zebraC(:,:);%./repmat(max(total_nond1_zebraC(:,:),[],2),1,72);
    
    nodom1 = [fornt_x_chd,fornt_y_chd];
    cAnch1 = [nodom1(1,:);nodom1(end,:)]; 
    pUtopia1 = min(cAnch1);
         
    dist = sqrt((fornt_x_chd-pUtopia1(1)).^2 + (fornt_y_chd-pUtopia1(2)).^2);
    [~, ind] = sort(dist);
    figure(f2)
    scatter(nodom1(ind(1),1),nodom1(ind(1),2),50,color2,"^",'DisplayName',['UP' plot_name])
    
    
    
%     % get_utopia parameter
%     Front_parameter_WTMU= Par(ia,:);
%     Utopia_parameter_WTMU = Front_parameter_WTMU(ind(1),:);
    
    selectd = [3,4;1,2;3,4;1,2;3,4;1,2;];
    ysim = [(nond1_drosWT(1,1:36)), (nond1_zebraWT(1,1:36)),...
        (nond1_drosWT(ind(1),1:36)), (nond1_zebraWT(ind(1),1:36)),...
        (nond1_drosWT(end,1:36)), (nond1_zebraWT(end,1:36));...
        (nond1_drosS(1,1:36)), (nond1_zebraC(1,1:36)),...
        (nond1_drosS(ind(1),1:36)), (nond1_zebraC(ind(1),1:36)),...
        (nond1_drosS(end,1:36)), (nond1_zebraC(end,1:36))];

    % selectd = [3,1]';
    % ysim = [(nond_drosWT(ind(1),1:36))', (nond_zebraWT(ind(1),1:36))'];
    xsim = [linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36);...
        linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36),...
        linspace(-275,0,36), linspace(-700,0,36)];
    plotParetoFigures(experData,xsim,ysim,selectd)

    [~,h1]=suplabel('DV axis', 'x'); 
    set(h1,'FontSize',16) 
    [~,h2]=suplabel('Normalized BMP intensity','y'); 
    set(h2,'FontSize',16)
    [~,h3]=suplabel([plot_name ' Mutant fitting'] ,'t');  
    set(h3,'FontSize',20)
    set(gcf, 'Position',  [100, 100, 700, 800])
    set(gcf,'color','w');

    % WT+CLF anchor points----------------------------------------------------
    nond2_drosWT= total_nond2_drosWT(:,:);%./repmat(max(total_nond2_drosWT(:,:),[],2),1,72);  
    nond2_zebraWT = total_nond2_zebraWT(:,:);%./repmat(max(total_nond2_zebraWT(:,:),[],2),1,72);
    nond2_drosS= total_nond2_drosS(:,:);%./repmat(2*max(total_nond2_drosS(:,:),[],2),1,72);  
    nond2_zebraC= total_nond2_zebraC(:,:);%./repmat(max(total_nond2_zebraC(:,:),[],2),1,72);
    
    nodom2 = [fornt_x_wtchd,fornt_y_wtchd];
    cAnch2 = [nodom2(1,:);nodom2(end,:)]; 
    pUtopia2 = min(cAnch2);
         
    dist = sqrt((fornt_x_wtchd-pUtopia2(1)).^2 + (fornt_y_wtchd-pUtopia2(2)).^2);
    [~, ind] = sort(dist);
    figure(f3)
    scatter(nodom2(ind(1),1),nodom2(ind(1),2),50,color3,"^",'DisplayName',['UP' plot_name])
    
    
    
%     % get_utopia parameter
%     Front_parameter_WTMU= Par(ia,:);
%     Utopia_parameter_WTMU = Front_parameter_WTMU(ind(1),:);
    
    selectd = [3,4;1,2;3,4;1,2;3,4;1,2;];
    ysim = [(nond2_drosWT(1,1:36)), (nond2_zebraWT(1,1:36)),...
        (nond2_drosWT(ind(1),1:36)), (nond2_zebraWT(ind(1),1:36)),...
        (nond2_drosWT(end,1:36)), (nond2_zebraWT(end,1:36));...
        (nond2_drosS(1,1:36)), (nond2_zebraC(1,1:36)),...
        (nond2_drosS(ind(1),1:36)), (nond2_zebraC(ind(1),1:36)),...
        (nond2_drosS(end,1:36)), (nond2_zebraC(end,1:36))];

    % selectd = [3,1]';
    % ysim = [(nond_drosWT(ind(1),1:36))', (nond_zebraWT(ind(1),1:36))'];
%     xsim = linspace(-180,0,36);
    plotParetoFigures(experData, xsim,ysim,selectd)

    [~,h1]=suplabel('DV axis', 'x'); 
    set(h1,'FontSize',16) 
    [~,h2]=suplabel('Normalized BMP intensity','y'); 
    set(h2,'FontSize',16)
    [~,h3]=suplabel([plot_name ' WT+Mutant fitting'] ,'t');  
    set(h3,'FontSize',20)
    set(gcf, 'Position',  [100, 100, 700, 800])
    set(gcf,'color','w');
    
    end
    
end