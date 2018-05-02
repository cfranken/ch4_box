function [ ] = plotMCF( St, model, data, baseName )
%%% Add errorbars?
add_error = true;
obsCol    = [  0,   0,   0]./256;
nhCol     = [204, 179, 102]./256;
shCol     = [ 58, 106, 176]./256;
eOpts     = {'-','Color','k','LineWidth',1};
yrs       = datevec(St);
yrs       = yrs(:,1);
xLims     = [datenum(yrs(1),1,1),datenum(yrs(end),1,1)];
%hold on

disp('MCF')
semilogy(St,data.nh_mcf,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
hold on
semilogy(St,data.sh_mcf,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
semilogy(St,model.nh_mcf,'-','Color',nhCol,'LineWidth',2)
semilogy(St,model.sh_mcf,'-','Color',shCol,'LineWidth',2)
%box on
set(gca,'YGrid','on','LineWidth',2,'TickDir','out','FontSize',12,'FontName','Helvetica')
%ylabel('CH_3CCl_3 (ppt)','FontSize',16)
xlim(xLims)
datetick('x','yyyy','keeplimits')
end