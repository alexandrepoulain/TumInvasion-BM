%% Script figures for the article
close all
savefigs = 0; % yes we save the figs in files in the directory specified   
plotpath = '.' % please enter here a valid path

tunit = 60*60; % Time unit for plotting (s): 1 hour

hfig = figure;  % save the figure handle in a variable
plot(t_BM/tunit,ct_BM,'b-','LineWidth',3,'DisplayName','TIMP-2');
hold on 
plot(t_BM/tunit,ca_BM,'r-','LineWidth',3,'DisplayName','MMP-2');
plot(t_BM/tunit,cp_BM,'g-','LineWidth',3,'DisplayName','proMMP-2');
grid on;
xlabel('time $t$ (h)')
ylabel('concentration (nM)')
xlim([0, max(t_BM/tunit)])
ylim([0 max(max(ct_BM), max(cp_BM))+5])
fname = strcat(test , '_concentrations_BM.eps');
fname = strcat(plotpath,fname)
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
legend('Location', 'east');
if savefigs
    saveas(hfig,fname, 'epsc')
end
%% BM density 
hfig = figure(2);  

fname = strcat(test , '_BM_density.eps');
fname = strcat(plotpath,fname)

plot(t_BM/tunit,M_BM,'b-','LineWidth',3,'DisplayName','BM density');
grid on;
xlabel('time $t$ (h)')
ylabel('density (nM)')
xlim([0, max(t_BM/tunit)])
ylim([0 Mmax])

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.8; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
legend('Location', 'east');;

if savefigs
    saveas(hfig,fname, 'epsc')
end

%% complexes
if test == 2 | test == 3
    hfig = figure(3);  % save the figure handle in a variable
    plot(t_BM/tunit,c3_BM,'g-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/TIMP-2');    
    
    grid on;
    xlabel('time $t$ (h)')
    ylabel('concentration (nM)')
    xlim([0, max(t_BM/tunit)])
    ylim([0 max(c3_BM)+20])
    fname = strcat(test , '_complexe3_BM.eps');
    fname = strcat(plotpath,fname)
    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio = 0.8; % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
    
    set(findall(hfig,'-property','Box'),'Box','off') % optional
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-painters','-fillpage')
    %print(hfig,fname,'-dpng','-painters')
    legend('Location', 'east');
    
    if savefigs
        saveas(hfig,fname, 'epsc')
    end

    hfig = figure(5);  % save the figure handle in a variable
    plot(t_BM/tunit,c1_BM,'b-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2');
    hold on 
    plot(t_BM/tunit,c2_BM,'r-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/proMMP-2');
    hold off
    grid on;
    xlabel('time $t$ (h)')
    ylabel('concentration (nM)')
    xlim([0, max(t_BM/tunit)])
    ylim([0 max(max(c2_BM), max(c1_BM))+0.05])
    fname = strcat(test , '_complexes12_BM.eps');
    fname = strcat(plotpath,fname)
    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio = 0.8; % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
    
    set(findall(hfig,'-property','Box'),'Box','off') % optional
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-painters','-fillpage')
    %print(hfig,fname,'-dpng','-painters')
    legend('Location', 'northeast');
    
    if savefigs
        saveas(hfig,fname, 'epsc')
    end

end


%% TIMP and proMMP in conjunctive
for ii = 1:size(t_BM)
    hfig = figure(12);  
    
    %fname = strcat(test , '_BM_density.eps');
    %fname = strcat(plotpath,fname)
    
    plot(x,ct_conj(ii,:),'b-','LineWidth',3,'DisplayName','TIMP-2');
    hold on;
    plot(x,cp_conj(ii,:),'g-','LineWidth',3,'DisplayName','proMMP-2');
    hold off;
    grid on;
    xlabel('Conjunctive tissue $x$ (dm)')
    ylabel('density (nM)')
    xlim([0, max(x)])
    ylim([0 25])
    %title(t_BM(ii)/tunit)
    picturewidth = 20; % set this parameter and keep it forever
    hw_ratio = 0.8; % feel free to play with this ratio
    set(findall(hfig,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
    
    set(findall(hfig,'-property','Box'),'Box','off') % optional
    set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
    set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
    set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
    pos = get(hfig,'Position');
    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
    %print(hfig,fname,'-dpdf','-painters','-fillpage')
    %print(hfig,fname,'-dpng','-painters')
    legend;
    if (ii == 1 | ii == 50 | ii == 100)
        fname = strcat(test , '_concentrations_conjunctive.eps');
        fname = strcat(num2str(t_BM(ii)/tunit),fname);
        fname = strcat(plotpath,fname);
        if savefigs
            saveas(hfig,fname, 'epsc')
        end

    end

end

