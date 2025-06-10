function plotParetoFigures(edata, xsim,ysim,selectd)
% Plots figures for data ve simulation:
%
%    edata - structure containing experimental data
%          - x: positional information
%          - y: yvalues 
%          - v: uncertainty information
%    xsim  - positional information of simulation
%    ysim  - yvalues from simulation
%parameters for figure and panel size
nd = size(selectd,1); 
np = ceil(nd/2);         % number of subplots
f=figure;
% subplot1(np,2,'Gap',[0.02 0.02],'XTickL','None','YTickL','All','FontS',16);
exp_name = [" Zeb WT"," Zeb CLF", " Dro WT"," Dro SLF"," Dro Chd like"," Dro Chd like CLF"];
for i = 1:nd
%     subplot1(i)
    colors = distinguishable_colors(5);
    subplot(np,2,i)
%     if size(selectd,2) == 2
%         errorbar(-edata{selectd(i,1)}.x(1,2:2:36),edata{selectd(i,1)}.y(1,2:2:36),...
%             edata{selectd(i,1)}.v(1,2:2:36),'linewidth',2,'color',colors(1,:));        
%         hold on
%         errorbar(-edata{selectd(i,2)}.x(1,2:2:36),edata{selectd(i,2)}.y(1,2:2:36),...
%             edata{selectd(i,2)}.v(1,2:2:36),'linewidth',2,'color',colors(2,:));
%         plot(xsim,ysim(:,i),'linewidth',2,'color',colors(3,:));
%         plot(xsim,ysim(:,i+nd),'linewidth',2,'color',colors(4,:));
%        
%         maxm = max(ysim(1,i),ysim(1,i+nd));
%         maxm = max(maxm,max(edata{selectd(i,2)}.y(1,2:2:36)));
%     elseif size(selectd,2) == 3
%         errorbar(-edata{selectd(i,1)}.x(1,2:2:36),edata{selectd(i,1)}.y(1,2:2:36),...
%             edata{selectd(i,1)}.v(1,2:2:36),'linewidth',2,'color',colors(1,:));        
%         hold on
%         errorbar(-edata{selectd(i,2)}.x(1,2:2:36),edata{selectd(i,2)}.y(1,2:2:36),...
%             edata{selectd(i,2)}.v(1,2:2:36),'linewidth',2,'color',colors(2,:));
%         errorbar(edata{selectd(i,3)}.x(1,2:2:36),edata{selectd(i,3)}.y(1,2:2:36),...
%             edata{selectd(i,3)}.v(1,2:2:36),'linewidth',2,'color',colors(1,:));
%         plot(xsim,ysim(:,i),'linewidth',2,'color',colors(3,:));
%         plot(xsim,ysim(:,i+nd),'linewidth',2,'color',colors(4,:));
%         plot(xsim,ysim(:,i+nd+nd),'linewidth',2,'color',colors(5,:));
%         maxm = max(ysim(1,i),ysim(1,i+nd));
%         maxm = max(maxm,max(edata{selectd(i,2)}.y));
%     else
%         errorbar(-edata{selectd(i)}.x(1,2:2:36),edata{selectd(i)}.y(1,2:2:36),...
%             edata{selectd(i)}.v(1,2:2:36),'linewidth',2,'color',colors(1,:));
%         hold on
%         plot(-xsim,ysim(:,i),'linewidth',2,'color',colors(2,:));
%         maxm = max(ysim(end,i),max(edata{selectd(i,1)}.y));
%     end    
    % labels1
    maxm = max(ysim(1,i),ysim(1,i+nd));
%     maxx = max(-xsim(,i));
    for j=1:size(selectd,2)
        x = -xsim(j,2+(i-1)*36:2:36+(i-1)*36)';
        lo = edata{selectd(i,j)}.y(1,2:2:36)' - edata{selectd(i,j)}.v(1,2:2:36)';
        hi = edata{selectd(i,j)}.y(1,2:2:36)' + edata{selectd(i,j)}.v(1,2:2:36)';
        hp = patch([x; x(end:-1:1); x(1)], [lo; hi(end:-1:1); lo(1)],colors(j,:),'facealpha',0.3,'DisplayName',strcat('Exp ',exp_name(selectd(i,j))));  
       
        hold on
        plot(-xsim(j,1+(i-1)*36:36+(i-1)*36),ysim(j,1+(i-1)*36:36+(i-1)*36),'linewidth',2,'color',colors(j,:),'DisplayName',strcat('Sim ',exp_name(selectd(i,j))));
%             edata{selectd(i,1)}.v(1,2:2:36)
        
        set(gca,'FontSize',10)
        legend() 
        maxm = max(ysim(end,i),max(edata{selectd(i,j)}.y));
        
        
    end
    
%     xlim([0.05 maxx]);
    ylim([0,1.2]);%ylim([-0.1 maxm+0.2]);

end
%tightfig
    
    
end
