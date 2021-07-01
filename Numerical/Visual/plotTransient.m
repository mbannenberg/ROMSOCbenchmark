% ------------------------------------------------------------------------------
% Function for plotting the transient
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

function h = plotTransient(t,y)
    h = figure(); hold on;
    plot(t,y,'LineWidth',2);
    grid on;
    xlabel('Time in ms');
    ylabel('Nodal voltages and currents');
    set(gca, 'FontName', 'Times New Roman','FontSize',14);
end