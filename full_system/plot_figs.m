%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Script figures for the article.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
savefigs = 1; % yes we save the figs in files in the directory specified   
dynamic_plot = 0; % if yes plot figures at all time steps
plotpath = '../results/results_full/' % Please enter here a valid path


tunit = 60*60; % Time unit for plotting (s): 1 hour

if test == 1
    test_case_str = "healthy";
elseif test ==2
    test_case_str = "tumor";
elseif test == 3
    test_case_str = "tumor+SF";
end

if dim == 1  
    

    hfig = figure(1);  % save the figure handle in a variable
    plot(t_BM/tunit,ct_BM,'b-','LineWidth',3,'DisplayName','TIMP-2');
    hold on 
    plot(t_BM/tunit,ca_BM,'r-','LineWidth',3,'DisplayName','MMP-2');
    plot(t_BM/tunit,cp_BM,'g-','LineWidth',3,'DisplayName','proMMP-2');
    grid on;
    xlabel('time $t$ (h)')
    ylabel('concentration (nM)')
    xlim([0, max(t_BM/tunit)])
    ylim([0 max(max(ct_BM), max(cp_BM))+5])
    fname = strcat(test_case_str , '_concentrations_BM.eps');
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
    
    fname = strcat(test_case_str , '_BM_density.eps');
    fname = strcat(plotpath,fname)
    
    plot(t_BM/tunit,M_BM,'b-','LineWidth',3,'DisplayName','BM density');
    grid on;
    xlabel('time $t$ (h)')
    ylabel('density (nM)')
    xlim([0, max(t_BM/tunit)])
    ylim([min(M_BM) Mmax])
    
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

    %% MT1-MMP
    hfig = figure(55);  % save the figure handle in a variable
    plot(t_BM/tunit,cm_BM,'-','Color',[0.4940 0.1840 0.5560],'LineWidth',3,'DisplayName','monomeric MT1-MMP');
    hold on 
    plot(t_BM/tunit,cd_BM,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP');
    hold off
    grid on;
    xlabel('time $t$ (h)')
    ylabel('concentration (nM)')
    xlim([0, max(t_BM/tunit)])
    ylim([0 max(max(cm_BM), max(cd_BM))+0.05])
    fname = strcat(test_case_str , '_MT1_MMP_BM.eps');
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
    
    %% complexes
    
    hfig = figure(3);  % save the figure handle in a variable
    plot(t_BM/tunit,c3_BM,'g-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/TIMP-2');    
    
    grid on;
    xlabel('time $t$ (h)')
    ylabel('concentration (nM)')
    xlim([0, max(t_BM/tunit)])
    ylim([0 max(c3_BM)+20])
    fname = strcat(test_case_str , '_complexe3_BM.eps');
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
    fname = strcat(test_case_str , '_complexes12_BM.eps');
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
    
    %% TIMP and proMMP in conjunctive
    for ii = 1:size(t_BM)
        hfig = figure(12);  
        
        %fname = strcat(test_case_str , '_BM_density.eps');
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
        if (ii == 1 | ii == 51 | ii == 100)
            fname = strcat(test_case_str , '_concentrations_conjunctive.eps');
            fname = strcat(num2str(t_BM(ii)/tunit),fname);
            fname = strcat(plotpath,fname);
            if savefigs == "yes"
                saveas(hfig,fname, 'epsc')
            end
    
        end
    
    end
elseif dim == 2

    for ii = 1:size(t_BM)
        %% concentrations in BM
        if dynamic_plot
            hfig1 = figure(1);  % save the figure handle in a variable
            plot(x,ct_BM(ii,:),'b-','LineWidth',3,'DisplayName','TIMP-2');
            hold on 
            plot(x,ca_BM(ii,:),'r-','LineWidth',3,'DisplayName','MMP-2');
            plot(x,cp_BM(ii,:),'g-','LineWidth',3,'DisplayName','proMMP-2');
            grid on;
            hold off
            xlabel('basement membrane length (in dm)')
            ylabel('concentration (nM)')
    
            xlim([0, max(y)])
            ylim([0 max(max(ct_BM(ii,:)), max(cp_BM(ii,:)))+5])
            
            picturewidth = 20; % set this parameter and keep it forever
            hw_ratio = 0.8; % feel free to play with this ratio
            set(findall(hfig1,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
            
            set(findall(hfig1,'-property','Box'),'Box','off') % optional
            set(findall(hfig1,'-property','Interpreter'),'Interpreter','latex') 
            set(findall(hfig1,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
            set(hfig1,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
            pos = get(hfig1,'Position');
            set(hfig1,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
            %print(hfig,fname,'-dpdf','-painters','-fillpage')
            %print(hfig,fname,'-dpng','-painters')
            legend('Location', 'east');
        end

        if (ii == 1 || ii == 51 || ii == 100)
            hfig11 = figure(11);
            if ii == 1
                plot(x,ct_BM(ii,:),'b-','LineWidth',3,'DisplayName','TIMP-2');
                hold on 
                plot(x,ca_BM(ii,:),'r-','LineWidth',3,'DisplayName','MMP-2');
                plot(x,cp_BM(ii,:),'g-','LineWidth',3,'DisplayName','proMMP-2');
            elseif ii == 51
                plot(x,ct_BM(ii,:),'b--','LineWidth',3,'DisplayName','TIMP-2','HandleVisibility','off');
                plot(x,ca_BM(ii,:),'r--','LineWidth',3,'DisplayName','MMP-2','HandleVisibility','off');
                plot(x,cp_BM(ii,:),'g--','LineWidth',3,'DisplayName','proMMP-2','HandleVisibility','off');
            elseif ii == 100
                plot(x,ct_BM(ii,:),'b:','LineWidth',3,'DisplayName','TIMP-2','HandleVisibility','off');
                plot(x,ca_BM(ii,:),'r:','LineWidth',3,'DisplayName','MMP-2','HandleVisibility','off');
                plot(x,cp_BM(ii,:),'g:','LineWidth',3,'DisplayName','proMMP-2','HandleVisibility','off');

                grid on;
                hold off
                xlabel('basement membrane length (in dm)')
                ylabel('concentration (nM)')
        
                xlim([0, max(x)])
                ylim([0 max(max(ct_BM,[], 'all'), max(cp_BM,[], 'all'))+5])
                
                picturewidth = 20; % set this parameter and keep it forever
                hw_ratio = 0.8; % feel free to play with this ratio
                set(findall(hfig11,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                
                set(findall(hfig11,'-property','Box'),'Box','off') % optional
                set(findall(hfig11,'-property','Interpreter'),'Interpreter','latex') 
                set(findall(hfig11,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                set(hfig11,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                pos = get(hfig11,'Position');
                set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                %print(hfig,fname,'-dpdf','-painters','-fillpage')
                %print(hfig,fname,'-dpng','-painters')
                legend('Location', 'east');
                ii
                if savefigs
                    fname = strcat(test_case_str , '_concentrations_BM.eps');
                    fname = strcat(plotpath,fname);
                    saveas(hfig11,fname, 'epsc')
                end
            end
        end
        
        pause(0.01)


        % zooming on active MMP
        if test == 2 | test == 3
            if (ii == 1 || ii == 51 || ii == 100)
                hfig11 = figure(111);
                if ii == 1
                    hold on 
                elseif ii == 51
                    plot(x,ca_BM(ii,:),'r--','LineWidth',3,'DisplayName','MMP-2');
                elseif ii == 100
                    plot(x,ca_BM(ii,:),'r:','LineWidth',3,'DisplayName','MMP-2','HandleVisibility','off');
    
                    grid on;
                    hold off
                    xlabel('basement membrane length (in dm)')
                    ylabel('concentration (nM)')
            
                    xlim([0, max(x)])
                    ylim([0 max(ca_BM(2:end,:),[], 'all')+0.001])
                    
                    picturewidth = 20; % set this parameter and keep it forever
                    hw_ratio = 0.8; % feel free to play with this ratio
                    set(findall(hfig11,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                    
                    set(findall(hfig11,'-property','Box'),'Box','off') % optional
                    set(findall(hfig11,'-property','Interpreter'),'Interpreter','latex') 
                    set(findall(hfig11,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                    set(hfig11,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                    pos = get(hfig11,'Position');
                    set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                    %print(hfig,fname,'-dpdf','-painters','-fillpage')
                    %print(hfig,fname,'-dpng','-painters')
                    legend('Location', 'northeast');
                    ii
                    if savefigs 
                        fname = strcat(test_case_str , '_active_MMP_concentrations_BM.eps');
                        fname = strcat(plotpath,fname);
                        saveas(hfig11,fname, 'epsc')
                    end
                end
            end
        end
        %% BM density
        if dynamic_plot
            hfig2 = figure(2);  
            plot(x,M_BM(ii,:),'b-','LineWidth',3,'DisplayName','BM density');
            grid on;
            xlabel('basement membrane length (in dm)')
            ylabel('density (nM)')
            xlim([0, max(x)])
            ylim([0 Mmax])
            
            picturewidth = 20; % set this parameter and keep it forever
            hw_ratio = 0.8; % feel free to play with this ratio
            set(findall(hfig2,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
            
            set(findall(hfig2,'-property','Box'),'Box','off') % optional
            set(findall(hfig2,'-property','Interpreter'),'Interpreter','latex') 
            set(findall(hfig2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
            set(hfig2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
            pos = get(hfig2,'Position');
            set(hfig2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
            %print(hfig,fname,'-dpdf','-painters','-fillpage')
            %print(hfig,fname,'-dpng','-painters')
            legend('Location', 'east');
        end

        if (ii == 1 || ii == 51 || ii == 100)
            hfig22 = figure(22);
            if ii == 1
                plot(x,M_BM(ii,:),'b-','LineWidth',3,'DisplayName','BM density');
                hold on
            elseif ii == 51
                plot(x,M_BM(ii,:),'b--','LineWidth',3,'DisplayName','BM density','HandleVisibility','off');

            elseif ii == 100
                plot(x,M_BM(ii,:),'b:','LineWidth',3,'DisplayName','BM density','HandleVisibility','off');

                hold off
                grid on;
                xlabel('basement membrane length (in dm)')
                ylabel('density (nM)')
                xlim([0, max(x)])
                ylim([0 Mmax])
                
                picturewidth = 20; % set this parameter and keep it forever
                hw_ratio = 0.8; % feel free to play with this ratio
                set(findall(hfig22,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                
                set(findall(hfig22,'-property','Box'),'Box','off') % optional
                set(findall(hfig22,'-property','Interpreter'),'Interpreter','latex') 
                set(findall(hfig22,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                set(hfig22,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                pos = get(hfig22,'Position');
                set(hfig22,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                %print(hfig,fname,'-dpdf','-painters','-fillpage')
                %print(hfig,fname,'-dpng','-painters')
                legend('Location', 'east');

                if savefigs
                    fname = strcat(test_case_str , '_BM_density.eps');
                    fname = strcat(plotpath,fname)
                    saveas(hfig22,fname, 'epsc')
                end
            end
        end
        
        %% MT1-MMP
        if dynamic_plot

            hfig3 = figure(3);  % save the figure handle in a variable
            plot(x,cm_BM(ii,:),'-','Color',[0.4940 0.1840 0.5560],'LineWidth',3,'DisplayName','monomeric MT1-MMP');
            hold on 
            plot(x,cd_BM(ii,:),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP');
            hold off
            grid on;
            xlabel('basement membrane length (in dm)')
            ylabel('concentration (nM)')
            xlim([0, max(x)])
            ylim([0 max(max(cm_BM(ii,:)), max(cd_BM(ii,:)))+0.05])
            
            picturewidth = 20; % set this parameter and keep it forever
            hw_ratio = 0.8; % feel free to play with this ratio
            set(findall(hfig3,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
            
            set(findall(hfig3,'-property','Box'),'Box','off') % optional
            set(findall(hfig3,'-property','Interpreter'),'Interpreter','latex') 
            set(findall(hfig3,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
            set(hfig3,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
            pos = get(hfig3,'Position');
            set(hfig3,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
            %print(hfig,fname,'-dpdf','-painters','-fillpage')
            %print(hfig,fname,'-dpng','-painters')
            legend('Location', 'northeast');
        end
        if (ii == 1 || ii == 51 || ii == 100)
            hfig33 = figure(33);
            if ii == 1
                plot(x,cm_BM(ii,:),'-','Color',[0.4940 0.1840 0.5560],'LineWidth',3,'DisplayName','monomeric MT1-MMP');
                hold on
                plot(x,cd_BM(ii,:),'-','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP');

            elseif ii == 51
                plot(x,cm_BM(ii,:),'--','Color',[0.4940 0.1840 0.5560],'LineWidth',3,'DisplayName','monomeric MT1-MMP','HandleVisibility','off');
                plot(x,cd_BM(ii,:),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP','HandleVisibility','off');
            elseif ii == 100
                plot(x,cm_BM(ii,:),':','Color',[0.4940 0.1840 0.5560],'LineWidth',3,'DisplayName','monomeric MT1-MMP','HandleVisibility','off');
                plot(x,cd_BM(ii,:),':','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP','HandleVisibility','off');
                hold off
                grid on;
                xlabel('basement membrane length (in dm)')
                ylabel('concentration (nM)')
                xlim([0, max(x)])
                ylim([0 max(max(cm_BM,[], 'all'), max(cd_BM,[], 'all'))+0.05])
                
                picturewidth = 20; % set this parameter and keep it forever
                hw_ratio = 0.8; % feel free to play with this ratio
                set(findall(hfig33,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                
                set(findall(hfig33,'-property','Box'),'Box','off') % optional
                set(findall(hfig33,'-property','Interpreter'),'Interpreter','latex') 
                set(findall(hfig33,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                set(hfig33,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                pos = get(hfig33,'Position');
                set(hfig33,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                %print(hfig,fname,'-dpdf','-painters','-fillpage')
                %print(hfig,fname,'-dpng','-painters')
                legend('Location', 'northeast');

                if savefigs
                    fname = strcat(test_case_str , '_MT1_MMP_BM.eps');
                    fname = strcat(plotpath,fname);
                    saveas(hfig33,fname, 'epsc')
                end
            end
        end

        %% zoom on dimeric MT1-MMP
        if (ii == 1 || ii == 51 || ii == 100)
            hfig333 = figure(333);
            if ii == 51
                plot(x,cd_BM(ii,:),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP','HandleVisibility','on');
                hold on;
            elseif ii == 100
                plot(x,cd_BM(ii,:),':','Color',[0.9290 0.6940 0.1250],'LineWidth',3,'DisplayName','dimeric MT1-MMP','HandleVisibility','off');
                hold off
                grid on;
                xlabel('basement membrane length (in dm)')
                ylabel('concentration (nM)')
                xlim([0, max(x)])
                ylim([0  max(cd_BM(50:end),[], 'all')])
                
                picturewidth = 20; % set this parameter and keep it forever
                hw_ratio = 0.8; % feel free to play with this ratio
                set(findall(hfig333,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                
                set(findall(hfig333,'-property','Box'),'Box','off') % optional
                set(findall(hfig333,'-property','Interpreter'),'Interpreter','latex') 
                set(findall(hfig333,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                set(hfig333,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                pos = get(hfig33,'Position');
                set(hfig333,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                %print(hfig,fname,'-dpdf','-painters','-fillpage')
                %print(hfig,fname,'-dpng','-painters')
                legend('Location', 'northeast');

                if savefigs
                    fname = strcat(test_case_str , '_dimeric_MT1_MMP_BM.eps');
                    fname = strcat(plotpath,fname);
                    saveas(hfig333,fname, 'epsc')
                end
            end
        end

        %% complexes
        if dynamic_plot

            hfig4 = figure(4);  % save the figure handle in a variable
            plot(x,c3_BM(ii,:),'g-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/TIMP-2');    
            grid on;
            xlabel('basement membrane length (in dm)')
            ylabel('concentration (nM)')
            xlim([0, max(x)])
            ylim([0 max(c3_BM,[],'all')+0.1*max(c3_BM,[],'all')])
            
            picturewidth = 20; % set this parameter and keep it forever
            hw_ratio = 0.8; % feel free to play with this ratio
            set(findall(hfig4,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
            
            set(findall(hfig4,'-property','Box'),'Box','off') % optional
            set(findall(hfig4,'-property','Interpreter'),'Interpreter','latex') 
            set(findall(hfig4,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
            set(hfig4,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
            pos = get(hfig4,'Position');
            set(hfig4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
            %print(hfig,fname,'-dpdf','-painters','-fillpage')
            %print(hfig,fname,'-dpng','-painters')
            legend('Location', 'northeast');
        end

        if (ii == 1 || ii == 51 || ii == 100)
            hfig44 = figure(44);
            if ii == 1
                plot(x,c3_BM(ii,:),'g-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/TIMP-2');
                hold on

            elseif ii == 51
                plot(x,c3_BM(ii,:),'g--','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/TIMP-2','HandleVisibility','off');
            elseif ii == 100
                plot(x,c3_BM(ii,:),'g:','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/TIMP-2','HandleVisibility','off');
                hold off
                grid on;
                xlabel('basement membrane length (in dm)')
                ylabel('concentration (nM)')
                xlim([0, max(x)])
                ylim([0 max(c3_BM,[],'all')+0.1*max(c3_BM,[],'all')])
                
                picturewidth = 20; % set this parameter and keep it forever
                hw_ratio = 0.8; % feel free to play with this ratio
                set(findall(hfig44,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                
                set(findall(hfig44,'-property','Box'),'Box','off') % optional
                set(findall(hfig44,'-property','Interpreter'),'Interpreter','latex') 
                set(findall(hfig44,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                set(hfig44,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                pos = get(hfig44,'Position');
                set(hfig44,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                %print(hfig,fname,'-dpdf','-painters','-fillpage')
                %print(hfig,fname,'-dpng','-painters')
                legend('Location', 'northeast');
                fname = strcat(test_case_str , '_complexe3_BM.eps');
                fname = strcat(plotpath,fname)
                if savefigs 
                    saveas(hfig44,fname, 'epsc')
                end
            end
        end
        if dynamic_plot

            hfig5 = figure(5);  % save the figure handle in a variable
            plot(x,c1_BM(ii,:),'b-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2');
            hold on 
            plot(x,c2_BM(ii,:),'r-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/proMMP-2');
            hold off
            grid on;
            xlabel('basement membrane length (in dm)')
            ylabel('concentration (nM)')
            xlim([0, max(x)])
            ylim([0 max(max(c2_BM(ii,:)), max(c1_BM(ii,:)))+0.05])
    
            picturewidth = 20; % set this parameter and keep it forever
            hw_ratio = 0.8; % feel free to play with this ratio
            set(findall(hfig5,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
            
            set(findall(hfig5,'-property','Box'),'Box','off') % optional
            set(findall(hfig5,'-property','Interpreter'),'Interpreter','latex') 
            set(findall(hfig5,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
            set(hfig5,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
            pos = get(hfig5,'Position');
            set(hfig5,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
            %print(hfig,fname,'-dpdf','-painters','-fillpage')
            %print(hfig,fname,'-dpng','-painters')
            legend('Location', 'northeast');
        end

        if (ii == 1 || ii == 51 || ii == 100)
            hfig55 = figure(55);
            if ii == 1
                plot(x,c1_BM(ii,:),'b-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2');
                hold on
                plot(x,c2_BM(ii,:),'r-','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/proMMP-2');

            elseif ii == 51
                plot(x,c1_BM(ii,:),'b--','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2','HandleVisibility','off');
                plot(x,c2_BM(ii,:),'r--','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/proMMP-2','HandleVisibility','off');
            elseif ii == 100
                plot(x,c1_BM(ii,:),'b:','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2','HandleVisibility','off');
                plot(x,c2_BM(ii,:),'r:','LineWidth',3,'DisplayName','MT1-MMP/TIMP-2/proMMP-2','HandleVisibility','off');
                hold off
                grid on;
                xlabel('basement membrane (in dm)')
                ylabel('concentration (nM)')
                xlim([0, max(x)])
                ylim([0 max(max(c2_BM,[], 'all'), max(c1_BM,[], 'all'))+0.05])
        
                picturewidth = 20; % set this parameter and keep it forever
                hw_ratio = 0.8; % feel free to play with this ratio
                set(findall(hfig55,'-property','FontSize'),'FontSize',21) % adjust fontsize to your document
                
                set(findall(hfig55,'-property','Box'),'Box','off') % optional
                set(findall(hfig55,'-property','Interpreter'),'Interpreter','latex') 
                set(findall(hfig55,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
                set(hfig55,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
                pos = get(hfig55,'Position');
                set(hfig55,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
                %print(hfig,fname,'-dpdf','-painters','-fillpage')
                %print(hfig,fname,'-dpng','-painters')
                legend('Location', 'northeast');
                fname = strcat(test_case_str , '_complexes12_BM.eps');
                fname = strcat(plotpath,fname)
                if savefigs
                    saveas(hfig55,fname, 'epsc')
                end
            end
        end
    
        %% TIMP and proMMP in conjunctive
        if dynamic_plot
            hfig = figure(6);  
    
            % reshape for plotting
            ct_conj_plot = reshape(ct_conj(ii,:),[Nx,Nx]);
            
            surf(X,Yx,ct_conj_plot','DisplayName','TIMP-2'); 
            hold on 
            for iSF =1:size(posSF,1)
                plot(posSF(1),posSF(2),'+')

            end
            view(2)
            shading interp
            clbr = colorbar
    
            %xlabel('x axis')
            %ylabel('y axis')
    
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
        end
        if (ii == 1 | ii == 51 | ii == 100)
            hfig = figure(6);  
    
            % reshape for plotting
            ct_conj_plot = reshape(ct_conj(ii,:),[Nx,Nx]);
            
            surf(X,Yx,ct_conj_plot','DisplayName','TIMP-2');
            view(2)
            shading interp
            hold on 
            if test == 1
                plot3(x(:),ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r-', 'LineWidth', 10); % solid line
            elseif test == 2| test == 3
                pos_dotted = ones(length(x),1).*(rho0==0);
                x_solid = x(pos_dotted>0);   
                if ii == 1
                    plot3(x(:),ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r-', 'LineWidth', 10); % solid line
                elseif ii == 51
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r:', 'LineWidth', 10); % dashed line
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 10); % solid line
                elseif ii == 100
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r.', 'LineWidth', 10); % dotted line
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 10); % solid line
                end
            end
            if test == 3
                for iSF =1:size(posSF,1)
                    plot3(posSF(1),posSF(2),100, 'k+', 'LineWidth', 3, 'MarkerSize', 10)
   
                end
            end
            hold off
            
            clbr = colorbar
    
            %xlabel('x axis')
            %ylabel('y axis')
    
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
            fname = strcat(test_case_str , '_TIMP-2_conjunctive.eps');
            fname = strcat(num2str(t_BM(ii)/tunit), fname );
            fname = strcat(plotpath,fname);
            if savefigs
                saveas(hfig,fname, 'epsc')
            end
    
        end
        if dynamic_plot
            hfig = figure(7);  
            % reshape for plotting
            cp_conj_plot = reshape(cp_conj(ii,:),[Nx,Nx]);
            
            surf(X,Yx,cp_conj_plot','DisplayName','proMMP-2');
            
            view(2)
            shading interp
            hold on 
            plot3(x,y(end), 100, 'r-', 'LineWidth', 3);
            if test == 3
                for iSF =1:size(posSF,1)
                    plot3(posSF(1),posSF(2),100, 'k+', 'LineWidth', 3, 'MarkerSize', 10)
   
                end
            end
            hold off
            clbr = colorbar
            
            %xlabel('x axis')
            %ylabel('y axis')
    
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
            %legend;
        end

        if (ii == 1 | ii == 51 | ii == 100)
            hfig = figure(7);  
            % reshape for plotting
            cp_conj_plot = reshape(cp_conj(ii,:),[Nx,Nx]);
            
            surf(X,Yx,cp_conj_plot','DisplayName','proMMP-2');
            view(2)
            shading interp
            hold on 
            if test == 1
                plot3(x(:),ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r-', 'LineWidth', 10); % solid line
            elseif test == 2| test == 3
                pos_dotted = ones(length(x),1).*(rho0==0);
                x_solid = x(pos_dotted>0);   
                if ii == 1
                    plot3(x(:),ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r-', 'LineWidth', 10); % solid line
                elseif ii == 51
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r:', 'LineWidth', 10); % dashed line
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 10); % solid line
                elseif ii == 100
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100, 'r.', 'LineWidth', 10); % dotted line
                    plot3(x,ones(size(x))*y(end)-0.0013, ones(size(x))*100.*pos_dotted, 'r-', 'LineWidth', 10); % solid line
                end
            end

            if test == 3
                for iSF =1:size(posSF,1)
                    plot3(posSF(1),posSF(2),100, 'k+', 'LineWidth', 3, 'MarkerSize', 10)
   
                end
            end
            hold off
            clbr = colorbar
    
            %xlabel('x axis')
            %ylabel('y axis')
    
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
            %legend;
            fname = strcat(test_case_str , '_proMMP-2_conjunctive.eps');
            fname = strcat(num2str(t_BM(ii)/tunit), fname );

            fname = strcat(plotpath,fname);
            if savefigs
                saveas(hfig,fname, 'epsc')
            end
    
        end

        
    end


end