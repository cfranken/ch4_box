%%% =======================================================================
%%% = plotObs.m
%%% = Alex Turner
%%% = 04/12/2016
%%% =----------------------------------------------------------------------
%%% = NOTES
%%% =  ( 1): Plot the observations and the box-model output.
%%% =----------------------------------------------------------------------
%%% = INPUTS
%%% =  ( 1): St       -- Our time vector.
%%% =  ( 2): model    -- Structure containing the model results.
%%% =  ( 3): data     -- Structure containing the observations.
%%% =  ( 4): baseName -- Prefix for the plots.
%%% =----------------------------------------------------------------------
%%% = OUTPUTS
%%% =   N/A
%%% =======================================================================

function [ ] = plotNewObs( St, model, data, baseName )

global interactive_OH

%%% Get file extension
fExten = strsplit(baseName,'.');
fExten = fExten(end);
if strcmp(fExten,'tif')
    printOpts = {'-dtiff','-r300'};
elseif strcmp(fExten,'eps')
    printOpts = {'-depsc2'};
elseif strcmp(fExten,'pdf')
    printOpts = {'-dpdf'};
elseif strcmp(fExten,'png')
    printOpts = {'-dpng'};
end

%%% Get parameters
params = getParameters(St);

%%% Add errorbars?
add_error = true;
obsCol    = [  0,   0,   0]./256;
nhCol     = [204, 179, 102]./256;
shCol     = [ 58, 106, 176]./256;
eOpts     = {'-','Color','k','LineWidth',1};
yrs       = datevec(St);
yrs       = yrs(:,1);
xLims     = [datenum(yrs(1),1,1),datenum(yrs(end),1,1)];

%%% All
h = figure();
set(gcf, 'Position', [100, 100, 800, 1000])
% CH4
subplot(5,1,1);
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4(i); yE = data.nh_ch4_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_ch4(i); yE = data.sh_ch4_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St,data.nh_ch4,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
plot(St,data.sh_ch4,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
plot(St,model.nh_ch4,'-','Color',nhCol,'LineWidth',2)
plot(St,model.sh_ch4,'-','Color',shCol,'LineWidth',2)
box on
set(gca,'YGrid','on','LineWidth',2,'TickDir','out','FontSize',12,'FontName','Helvetica')
ylabel('CH_4 (ppb)','FontSize',16)
xlim(xLims)
datetick('x','yyyy','keeplimits')
% CO
subplot(5,1,2);
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_co(i); yE = data.nh_co_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_co(i); yE = data.sh_co_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St,data.nh_co,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
plot(St,data.sh_co,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
plot(St,model.nh_co,'-','Color',nhCol,'LineWidth',2)
plot(St,model.sh_co,'-','Color',shCol,'LineWidth',2)
box on
set(gca,'YGrid','on','LineWidth',2,'TickDir','out','FontSize',12,'FontName','Helvetica')
ylabel('CO (ppb)','FontSize',16)
xlim(xLims)
datetick('x','yyyy','keeplimits')
% OH
subplot(5,1,3);
hold on
if interactive_OH
    plot(St,(model.nh_oh-params.gmOH)/params.gmOH*100,'-','Color',nhCol,'LineWidth',2)
    plot(St,(model.sh_oh-params.gmOH)/params.gmOH*100,'-','Color',shCol,'LineWidth',2)
else
    plot(St,(model.nh_oh-1)*100,'-','Color',nhCol,'LineWidth',2)
    plot(St,(model.sh_oh-1)*100,'-','Color',shCol,'LineWidth',2)
end
box on
set(gca,'YGrid','on','LineWidth',2,'TickDir','out','FontSize',12,'FontName','Helvetica')
ylabel('OH anomaly (%)','FontSize',16)
xlim(xLims)
datetick('x','yyyy','keeplimits')
% d13C
subplot(5,1,4);
hold on
if add_error
    for i = 1:length(St)
        yV = data.nh_ch4c13(i); yE = data.nh_ch4c13_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_ch4c13(i); yE = data.sh_ch4c13_err(i);
        if ~isnan(yV) && ~isnan(yE)
            plot([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
plot(St,data.nh_ch4c13,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
plot(St,data.sh_ch4c13,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
plot(St,model.nh_ch4c13,'-','Color',nhCol,'LineWidth',2)
plot(St,model.sh_ch4c13,'-','Color',shCol,'LineWidth',2)
box on
set(gca,'YGrid','on','LineWidth',2,'TickDir','out','FontSize',12,'FontName','Helvetica')
ylabel(sprintf('\\delta^{13}CH_{4} (%s)',char(8240)),'FontSize',16)
xlim(xLims)
datetick('x','yyyy','keeplimits')
% MCF
subplot(5,1,5);

if add_error
    for i = 1:length(St)
        yV = data.nh_mcf(i); yE = data.nh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            semilogy([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
        yV = data.sh_mcf(i); yE = data.sh_mcf_err(i);
        if ~isnan(yV) && ~isnan(yE)
            semilogy([St(i),St(i)],[yV-yE,yV+yE],eOpts{:})
        end
    end
end
hold on
%disp('MCF')
semilogy(St,data.nh_mcf,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
semilogy(St,data.sh_mcf,'^','Color',obsCol,'MarkerSize',6,'MarkerFaceColor',obsCol,'MarkerEdgeColor','none')
semilogy(St,model.nh_mcf,'-','Color',nhCol,'LineWidth',2)
semilogy(St,model.sh_mcf,'-','Color',shCol,'LineWidth',2)
%box on
set(gca,'YGrid','on','LineWidth',2,'TickDir','out','FontSize',12,'FontName','Helvetica')
ylabel('CH_3CCl_3 (ppt)','FontSize',16)
xlim(xLims)
datetick('x','yyyy','keeplimits')
print(h,printOpts{:},sprintf(baseName,'Plot'));

end


%%% =======================================================================
%%% = END
%%% =======================================================================
