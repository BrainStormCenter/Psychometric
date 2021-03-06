function createfigure(cdata1, cdata2, cdata3)
%CREATEFIGURE(CDATA1, CDATA2, CDATA3)
%  CDATA1:  image cdata
%  CDATA2:  image cdata
%  CDATA3:  image cdata

%  Auto-generated by MATLAB on 01-Feb-2018 13:38:10

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,...
    'Color',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'YTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
    'Layer','top',...
    'YDir','reverse',...
    'Position',[0.13 0.583837209302326 0.334659090909091 0.341162790697675]);
box(axes1,'on');
hold(axes1,'on');

% Create image
image(cdata1,'Parent',axes1);

% Create title
title({'HC'},'FontSize',11);

% Create subplot
subplot1 = subplot(2,2,2,'Parent',figure1,...
    'Color',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'YTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
    'Layer','top',...
    'YDir','reverse');
box(subplot1,'on');
hold(subplot1,'on');

% Create image
image(cdata2,'Parent',subplot1);

% Create title
title({'CLBP'},'FontSize',11);

% Create axes
axes2 = axes('Parent',figure1,...
    'Color',[0.501960813999176 0.501960813999176 0.501960813999176],...
    'YTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
    'XTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16],...
    'Layer','top',...
    'YDir','reverse',...
    'Position',[0.127142857142857 0.11 0.334659090909091 0.341162790697674]);
box(axes2,'on');
hold(axes2,'on');

% Create image
image(cdata3,'Parent',axes2,'CDataMapping','scaled');

% Create title
title({'Group differences in pain node connectivity','brighter = stronger HC, darker = stronger CLBP'},...
    'FontSize',11);

% Create textbox
annotation(figure1,'textbox',...
    [0.571616113744078 0.117381993041789 0.110055423594616 0.334519572953737],...
    'String',{'Pain_RightAmy','Pain_LeftAmy','Pain_RightThal','Pain_LeftThal','Pain_RightACC','Pain_LeftACC','Pain_RightAIns','Pain_LeftAIns','Pain_RightMCC','Pain_LeftMCC','Pain_RightPCC','Pain_LeftPCC','Pain_RightPIns','Pain_LeftPIns','Pain_RightSMA','Pain_LeftSMA'},...
    'Interpreter','none',...
    'FontSize',16);

