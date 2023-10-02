clc
clear all
close all

simdata = readtable('results/design_var_raninit_20231001_194514.xlsx');

simdata.Obj = round(simdata.Obj,5);

figure1 = figure;


coordvars = {'capture_width', 'pen_diameter', 'pen_height', 'stock_density', 'Obj', 'fish_yield_cons_ineq', 'sustainable_power_operation_cons'};

% simdata_new = simdata(~((simdata.Obj == 0) | (simdata.Obj > 10) | (simdata.fish_yield_cons_ineq > 1e-4)),:); % filter and remove 
% simdata_removal = simdata(((simdata.Obj == 0) | (simdata.Obj > 10)),:); % filter and remove 
% simdata_unsuccessful = simdata((simdata.Obj == 0),:); % filter and remove 
simdata_successful = simdata((simdata.success == 1),:); 

p = parallelplot(simdata_successful,'CoordinateVariables',coordvars, 'LineWidth',2, 'Color',{'#5F5F5F'});
p.CoordinateTickLabels = {'Capture Width [m]', 'Pen Diameter [m]', 'Pen Height [m]', 'Stocking Density [kg/m^3]', 'Cost / Fish Yield [$/kg]', 'Power Supply [W]', 'Fish Yield [kg]'};
p.CoordinateTickLabels = {''};
p.DataNormalization = 'range';
p.Jitter = 0;
p.LineAlpha = [0.2 0.8 0.4];
p.FontName = 'Arial';
p.FontSize = 14;


% Create textbox
annotation(figure1,'textbox',...
    [0.58 0.94 0.098 0.052],...
    'String','Objective Function',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor',[0.901960784313726 0.901960784313726 0.901960784313726]);

% Create textbox
annotation(figure1,'textbox',...
    [0.72 0.94 0.14 0.052],...
    'String','Constraints',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor',[0.901960784313726 0.901960784313726 0.901960784313726]);

% Create textbox
annotation(figure1,'textbox',...
    [0.16 0.94 0.37 0.05],...
    'String','Design Variables',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor',[0.901960784313726 0.901960784313726 0.901960784313726]);


% X labels:

% Create textbox
annotation(figure1,'textbox',[0.131333333333333 0.86 0.11 0.052],...
    'String','Capture Width [m]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.241333333333333 0.86 0.11 0.052],...
    'String','Pen Diameter [m]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.351333333333333 0.86 0.11 0.052],...
    'String','Pen Height [m]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.456666666666665 0.87 0.12 0.052],...
    'String','Stocking Density [kg/m^3]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.57333333333333 0.86 0.11 0.052],...
    'String','Cost / Fish Yield [$/kg]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.683333333333331 0.86 0.11 0.052],...
    'String','Norm. Power Supply [-]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.793999999999998 0.86 0.11 0.052],...
    'String','Norm. Fish Yield [-]',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Create textbox
annotation(figure1,'textbox',[0.16 0.94 0.37 0.05],...
    'String','Design Variables',...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'BackgroundColor',[0.901960784313726 0.901960784313726 0.901960784313726]);


x0=100;
y0=100;
width=1500;
height=500;
set(gcf,'position',[x0,y0,width,height])

set(findobj(gcf,'type','axes'),'FontName','Arial','FontWeight','Bold', 'FontSize', 14);

pos=get(gca,'position');  % retrieve the current values
pos(4)=0.9*pos(4);        % try reducing width 10%
set(gca,'position',pos);  % write the new values