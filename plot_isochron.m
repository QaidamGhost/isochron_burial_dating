function plot_isochron(a,b,sigma_b,linearized_data,data,removed_data,...
    init_Rinh)

% Plot isochron line, production ratio line, and errorbars of samples
%
% Arguments:
% a: intercept of isochron line (atom/g/yr; scalar)
% b: slope of isochron line (unitless; scalar)
% sigma_b: 1s error of slope (unitless; scalar)
% linearized_data must have fields as:
%   linearized_data.x: linearized 10Be concentration (atom/g; nx1 vector)
%   linearized_data.dx: 1s error of linearized 10Be (atom/g; nx1 vector)
%   linearized_data.y: measured 26Al concentration (atom/g; nx1 vector)
%   linearized_data.dy: 1s error of 26Al (atom/g; nx1 vector)
% data: measured data precluding reworked clasts (ditto)
% removed_data: measured data of reworked clasts (ditto)
% init_Rinh: initial guess of Rinh (unitless; scalar)

    % find the range of axis x and y
    total_x =[linearized_data.x data.x removed_data.x 0];
    total_dx=[linearized_data.dx data.dx removed_data.dx 0];
    total_y =[data.y removed_data.y a 0];
    total_dy=[data.dy removed_data.dy 0 0];
    [max_x,~]=max(total_x+total_dx);
    [min_x,~]=min(total_x-total_dx);
    [max_y,~]=max(total_y+total_dy);
    [min_y,~]=min(total_y-total_dy);

    if min_x==0
        X = linspace(0,max_x*1.1);
    else
        X = linspace(min_x*1.1-max_x*.1,max_x*1.1-min_x*.1);
    end
    % isochron line
    Y1 = b*X+a;
    % upper error line
    Y2 = (b+sigma_b)*X+a;
    % lower error line
    Y3 = (b-sigma_b)*X+a;
    % production ratio line
    Y4 = init_Rinh*X;
    plot(X,Y1,'k',X,Y2,'r--',X,Y3,'b--',X,Y4,'m'),
    xlabel('10Be Concentration (atom/g)'),
    ylabel('26Al Concentration (atom/g)'),
    grid on;

    % reset the range of axis x and y
    if min_x==0
        if min_y==0
            axis([0,max_x*1.1,0,max_y*1.1]);
        else
            axis([0,max_x*1.1,min_y*1.1-max_y*.1,max_y*1.1-min_y*.1]);
        end
    else
        if min_y==0
            axis([min_x*1.1-max_x*.1,max_x*1.1-min_x*.1,0,max_y*1.1]);
        else
            axis([min_x*1.1-max_x*.1,max_x*1.1-min_x*.1,min_y*1.1-max_y...
                *.1,max_y*1.1-min_y*.1]);
        end
    end
    % plot errorbar
    hold on;
    errorbar(linearized_data.x,linearized_data.y,linearized_data.dx,...
        'horizontal','ko');
    errorbar(linearized_data.x,linearized_data.y,linearized_data.dy,'ko');
    errorbar(data.x,data.y,data.dx,'horizontal','o','Color',[.5 .5 .5]);
    errorbar(data.x,data.y,data.dy,'o','Color',[.5 .5 .5]);
    errorbar(removed_data.x,removed_data.y,removed_data.dx,'horizontal',...
        ':o','Color',[.5 .5 .5]);
    errorbar(removed_data.x,removed_data.y,removed_data.dy,':o','Color',...
        [.5 .5 .5]);
    hold off;
end