clear; clc;
close all
clear all
% Random seed
% RandStream.setGlobalStream(RandStream('mt19937ar','Seed',0));
% addpath('/Volumes/MyBook2/zebrafish vs drosophila/ModelCompare')
addpath('./Function')
addpath('./ExpData')
%%
%% Load data
folder_path = 'Path to the simulation result';
save_folder = 'front save folder path';    


if ~exist(save_folder, 'dir')
       mkdir(save_folder)
end


save_output = 0 ; % check =1 if save all the SSE and WT and Chd/S mutant result 

nres =10;   %Screen file number in simple/Extended
n = 1000;    %Screen cases number in each file simple/Extended
d = 11;      %Screen parameter number
znames = {'k2','decB','decS','decBS','DB','DS','DBS','lambdaS',...
           'lambdaBS','j1','j2'};

f1 = figure;
hold on
f2 = figure;
hold on
f3 = figure;
hold on

cmap = colormap(jet(d+1));
cmap = cmap(randperm(length(cmap)),:);

modeltype ='Simple';
folders_names={};

% total_folder = 1+11+11*10/2+11*10*9/3/2;
folder_count_tot = 0;

for m = 0%'0 == Origianl screen'; '1 == 1D cases'; '2 == 2D cases'; '3 == 3D cases'
% datapath = './Result_120321_simple/SSE_scenario_Simple_';
    
    if m == 0
         
       folders_names{1} = '';
       folder_n = 1;
       plot_name = 'Degree0';
                                   
    elseif m == 1
        
   
       folder_n = 11;
       for i = 1:size(znames,2)
           folders_names{i} = [znames{i}];
       end
       plot_name = 'Degree1';
       
    elseif m == 2
        
       folder_n = 11*10/2;
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
       
    elseif m == 3
        
       folder_n = 11*10*9/3/2;
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
    end
       
                
                
       
   

    cnt = 1;     % zebrafish file number
    temp_i = 1;
    folder_count = 0;

    for f = 1:folder_n
        drosSSE = zeros(nres*n,5);
        zebraSSE = zeros(nres*n,12);
        Par = zeros(nres*n,d);
        drosWT = zeros(nres*n,72);
        drosS = zeros(nres*n,72);
        drosSH = zeros(nres*n,72);
        zebraWT = zeros(nres*n,72);
        zebraWT1 = zeros(nres*n,72);
        zebraWT2 = zeros(nres*n,72);
        zebraWT3 = zeros(nres*n,72);
        zebraWT4 = zeros(nres*n,72);
        zebraWT5 = zeros(nres*n,72);
        zebraC = zeros(nres*n,72);
        zebraCh = zeros(nres*n,72);
        zebraN = zeros(nres*n,72);
       
       datapath =  strcat(folder_path, ['Result_Simple_' folders_names{f}], '/');
       dataname =  ['SSE_scenario_Simple_' folders_names{f}];
%            fw = waitbar(0,'Total Pragrass');
       disp(['Processing folder' datapath])
       folder_start_n = 0;
       parcnt = 1;     % Results counter
        for i = 1:nres 
%                 waitbar((folder_count*nres+i)/(nres*folder_n),fw,'Total Pragrass');
             %filename = strcat('../ModelCompare/MOO_simplemodel_10_30_2017/SSE_scenario_SIMPLEMODELS_',num2str(i),'.mat');
        %     filename = strcat('../ModelCompare/SSE_scenario_EXTENDEDMODELS_H_',num2str(i),'.mat');
            filename = strcat(datapath,dataname,num2str(i),'.mat');

            % check if file exist
            if isfile(filename)
                load(filename);
%             elseif isfile(strcat(datapath,'SSE_scenario_Extanded_',num2str(p),'.mat'))
%                 load(strcat(datapath,'SSE_scenario_Extanded_',num2str(p),'.mat'));
            else
                disp(['File not exist :' filename]);
            end

        %      load(strcat('../SSE_scenario_EXTENDEDMODELS_H_',num2str(cnt),'.mat'));
            drosSSE(parcnt:parcnt+n-1,:) = SSE(:,13:17);
%                 Par(parcnt+folder_start_n:parcnt+n-1+folder_start_n,1:11)  = parameters_store;
            drosWT(parcnt:parcnt+n-1,:)  = Stored_Soln.dwt;     % BMP signal drosophila
            drosS(parcnt:parcnt+n-1,:)  = Stored_Soln.dshomo;     % BMP signal drosophila
        %     drosSH(parcnt:parcnt+n-1,:)  = Stored_Soln.dshet;     % BMP signal drosophila

            % zebrafish    
            zebraSSE(parcnt:parcnt+n-1,:)  = SSE(:,1:12);
            Vard(parcnt+folder_start_n:parcnt+n-1+folder_start_n,:) = Bd;

% 
             zebraWT(parcnt:parcnt+n-1,:)  = Stored_Soln.zwt57;     % BMP signal drosophila
%             zebraWT1(parcnt+folder_start_n:parcnt+n-1+folder_start_n,:)  = Stored_Soln.zwt47;     % BMP signal drosophila
%             zebraWT2(parcnt+folder_start_n:parcnt+n-1+folder_start_n,:)  = Stored_Soln.zwt53;     % BMP signal drosophila
%             zebraWT3(parcnt+folder_start_n:parcnt+n-1+folder_start_n,:)  = Stored_Soln.zwt57;     % BMP signal drosophila
%             zebraWT4(parcnt+folder_start_n:parcnt+n-1+folder_start_n,:)  = Stored_Soln.zwt63;     % BMP signal drosophila
%             zebraWT5(parcnt+folder_start_n:parcnt+n-1+folder_start_n,:)  = Stored_Soln.zwt67;     % BMP signal drosophila
             zebraC(parcnt:parcnt+n-1,:)  = Stored_Soln.zclf57;     % BMP signal drosophila
        %     zebraN(parcnt:parcnt+n-1,:)  = Stored_Soln.zhclf57;     % BMP signal drosophila
        %     zebraCh(parcnt:parcnt+n-1,:)  = Stored_Soln.znlf57;     % BMP signal drosophila

            cnt = cnt + 1;
            parcnt = parcnt+n;

        end
        
         folder_count = folder_count+1;
         folder_count_tot = folder_count_tot+1;
         % disp(['progess ' num2str(folder_count_tot/total_folder*100) '%'])
    
    
%     figure(f1)
    % plot WT
        WTzebr = zebraSSE(:,3);
        WTfly = drosSSE(:,1);
        [ndx, ndy, ia] = find_nondominated(WTfly, WTzebr);
        nond_drosWT= drosWT(ia,:)./repmat(max(drosWT(ia,:),[],2),1,72);  
        nond_zebraWT = zebraWT(ia,:)./repmat(max(zebraWT(ia,:),[],2),1,72);
        nond_drosS= drosS(ia,:)./repmat(max(drosWT(ia,:),[],2),1,72);  
        nond_zebraC= zebraC(ia,:)./repmat(max(zebraWT(ia,:),[],2),1,72);
        bd_dros = drosBd(ia,:);
        bd_zebra = zebraBd(ia,:);
     
        
        % scatter plot for all solutions
        figure(f1)
        scatter(WTfly,WTzebr,'k','filled')
        xlim([0,1])
        ylim([0,1])
    
%     plot(ndx, ndy, '-o','LineWidth',2.5,'DisplayName',plot_name)
%     set(gca,'FontSize',15)
%     title('WT front')
%     xlim([0,0.6])
%     ylim([0,0.2])
%     xlabel('Drosophila RMSE')
%     ylabel('Zebrafish RMSE')
%     legend()
%     nodom = [ndx,ndy];
% %     
%     cAnch = [nodom(1,:);nodom(end,:)]; 
%     pUtopia = min(cAnch);
%     
%     
%      
%     dist = sqrt((ndx-pUtopia(1)).^2 + (ndy-pUtopia(2)).^2);
%     [~, ind] = sort(dist);
%     load('experData'); % experimental data
%     
%     % get_utopia parameter
%     Front_parameter_WTMU= Par(ia,:);
%     Utopia_parameter_WTMU = Front_parameter_WTMU(ind(1),:);
%     
%     selectd = [3,4;1,2;3,4;1,2;3,4;1,2;];
%     ysim = [(nond_drosWT(1,1:36)), (nond_zebraWT(1,1:36)),...
%         (nond_drosWT(ind(1),1:36)), (nond_zebraWT(ind(1),1:36)),...
%         (nond_drosWT(end,1:36)), (nond_zebraWT(end,1:36));...
%         (nond_drosS(1,1:36)), (nond_zebraC(1,1:36)),...
%         (nond_drosS(ind(1),1:36)), (nond_zebraC(ind(1),1:36)),...
%         (nond_drosS(end,1:36)), (nond_zebraC(end,1:36))];
% 
%     selectd = [3,1]';
%     % ysim = [(nond_drosWT(ind(1),1:36))', (nond_zebraWT(ind(1),1:36))'];
%     xsim = linspace(-180,0,36);
%     plotParetoFigures(experData, xsim,ysim,selectd)
% 
%     [~,h1]=suplabel('DV axis (angle)', 'x'); 
%     set(h1,'FontSize',16) 
%     [~,h2]=suplabel('Normalized BMP intensity','y'); 
%     set(h2,'FontSize',16)
%     [~,h3]=suplabel([modeltype ' model : WT fitting'] ,'t');  
%     set(h3,'FontSize',20)
%     set(gcf, 'Position',  [100, 100, 700, 800])
%     set(gcf,'color','w');

    % plot S/C mutant
    Muzebr = (zebraSSE(:,8));
    Mufly = (drosSSE(:,2)); 
    [ndx1, ndy1, ia1] = find_nondominated(Mufly, Muzebr);
    nond1_drosWT= drosWT(ia1,:)./repmat(max(drosWT(ia1,:),[],2),1,72);  
    nond1_zebraWT = zebraWT(ia1,:)./repmat(max(zebraWT(ia1,:),[],2),1,72);
    nond1_drosS= drosS(ia1,:)./repmat(max(drosWT(ia1,:),[],2),1,72);  
    nond1_zebraC= zebraC(ia1,:)./repmat(max(zebraWT(ia1,:),[],2),1,72);
    bd_dros1 = drosBd(ia1,:);
    bd_zebra1 = zebraBd(ia1,:);
    
    figure(f2)
    scatter(Mufly,Muzebr,'k','filled');
    xlim([0,1])
    ylim([0,1])
%     plot(ndx1, ndy1, '-*','LineWidth',2.5,'DisplayName',plot_name)
%     set(gca,'FontSize',15)
%     title('CLF front')
%     xlim([0,0.6])
%     ylim([0,0.2])
%     xlabel('Drosophila RMSE')
%     ylabel('Zebrafish RMSE')
%     legend()

% 
%   
    WTmuzebra = (WTzebr+ Muzebr);
    WTmufly = (WTfly +Mufly );  
    [ndx2, ndy2, ia2] = find_nondominated(WTmufly, WTmuzebra);
    nond2_drosWT= drosWT(ia2,:)./repmat(max(drosWT(ia2,:),[],2),1,72);  
    nond2_zebraWT = zebraWT(ia2,:)./repmat(max(zebraWT(ia2,:),[],2),1,72);
    nond2_drosS= drosS(ia2,:)./repmat(max(drosWT(ia2,:),[],2),1,72);  
    nond2_zebraC= zebraC(ia2,:)./repmat(max(zebraWT(ia2,:),[],2),1,72);
    bd_dros2 = drosBd(ia2,:);
    bd_zebra2 = zebraBd(ia2,:);
    
    
    figure(f3)
    scatter(WTmufly,WTmuzebra,'k','filled');
    xlim([0,1])
    ylim([0,1])
%     plot(ndx2, ndy2, '-*','LineWidth',2.5,'DisplayName',plot_name)
%     set(gca,'FontSize',15)
%     title('WT+CLF front')
%     xlim([0,0.6])
%     ylim([0,0.5])
%     xlabel('Drosophila RMSE')
%     ylabel('Zebrafish RMSE')
%     legend()
   
    
    save([save_folder '/Front_' folders_names{f} '.mat'],'ndx','ndy','ia','ndx1','ndy1','ia1','ndx2','ndy2','ia2','nond_drosWT','nond1_drosWT','nond2_drosWT',...
        'nond_zebraWT','nond1_zebraWT','nond2_zebraWT','nond_drosS','nond1_drosS','nond2_drosS','nond_zebraC','nond1_zebraC','nond2_zebraC')
    end

end
%  
%     savefig(f1,'WT.fig')
%     savefig(f2,'CLF.fig')
%     savefig(f3,'WTCLF.fig')
                % errorbarxy(ndx,ndy,vX,vY)

    
 