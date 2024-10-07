function plot_isochron(a,sigma_a,b,sigma_b,linearized_data,data,removed_data,Rs,option)

%% Plot the isochron line, production ratio line, and errorbars of samples and upper limit of the post-burial concentration

%% Arguments:
% a: the intercept of the isochron line (atom/g/yr; scalar)
% sigma_a: 1 sigma absolute error of the intercept (atom/g/yr; scalar)
% b: the slope of the isochron line (unitless; scalar)
% sigma_b: 1 sigma absolute error of the slope (unitless; scalar)
% linearized_data must have fields as:
%   linearized_data.x: linearized 10Be concentration (atom/g; nx1 vector)
%   linearized_data.dx: 1 sigma absolute error of 10Be (atom/g; nx1 vector)
%   linearized_data.y: linearized 26Al concentration (atom/g; nx1 vector)
%   linearized_data.dy: 1 sigma absolute error of 26Al (atom/g; nx1 vector)
% data: measured data precluding reworked clasts (ditto)
% removed_data: measured data of reworked clasts (ditto)
% Rs: production rate ratio of 26Al versus 10Be at the surface (unitless;
% scalar)
% option have fields as:
%   option.flag: "0" for default usage, "1" for lower limit , and "2" for 
%   upper limit (unitless; scalar)
%   option.Npb:
%       option.Npb.x: post-burial 10Be concentration (atom/g; scalar)
%       option.Npb.dx: 1 sigma absolute error of 10Be (atom/g; scalar)
%       option.Npb.y: post-burial 26Al concentration (atom/g; scalar)
%       option.Npb.dy: 1 sigma absolute error of 26Al (atom/g; scalar)
%   option.Rd: production rate ratio of 26Al versus 10Be at the sampling
%   depth (unitless; scalar)

%% Output:
% void

    load consts.mat simulation_times;

    % find the range of axis x and y
    if removed_data.x(:,1)==-1
        total_x =[linearized_data.x data.x 0];
        total_dx=[linearized_data.dx data.dx 0];
        total_y =[data.y a 0];
        total_dy=[data.dy 0 0];
    else
        total_x =[linearized_data.x data.x removed_data.x 0];
        total_dx=[linearized_data.dx data.dx removed_data.dx 0];
        total_y =[data.y removed_data.y a 0];
        total_dy=[data.dy removed_data.dy 0 0];
    end
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
    % production ratio line
    Y2 = Rs*X;

    % step 1: set the title and axis and plot isochron line and production
    % rate ratio line
    if option.flag==0
        fig=figure('Name','Original Isochron Burial Dating');
        title('Original Isochron Burial Dating'),
        plot(X,Y1,'k',X,Y2,'m'),
    elseif option.flag==1
        fig=figure('Name','Lower Limit of the Burial Age');
        title('Lower Limit of the Burial Age'),
        plot(X,Y1,'b',X,Y2,'m'),
    elseif option.flag==2
        Y0 = option.Rd*X;
        fig=figure('Name','Upper Limit of the Burial Age');
        title('Upper Limit of the Burial Age'),
        plot(X,Y1,'r',X,Y2,'m',X,Y0,'g'),
    end
    xlabel('10Be Concentration (atom/g)'),
    ylabel('26Al Concentration (atom/g)'),
    grid on;

    hold on;
    % step 2: plot upper and lower errors of the isochron line
    % initializating
    rand_X=zeros(size(linearized_data.x));
    rand_Y=rand_X;
    Y3=zeros(size(X));    % upper error line
    Y3(1)=nan;
    for i=1:simulation_times
        while true
            for j=1:max(size(linearized_data.x))
                rand_X(j)=normrnd(linearized_data.x(j),linearized_data.dx(j));
                rand_Y(j)=normrnd(linearized_data.y(j),linearized_data.dy(j));
            end
            % determine random intercept and slope
            if option.flag==0 || option.flag==2
                coef=ones(1,length(rand_X));
                rand_X2=[rand_X;coef]';
                rand_Y2=rand_Y';
                tmp=(rand_X2'*rand_X2)\rand_X2'*rand_Y2;
                rand_b=tmp(1);
                rand_a=tmp(2);
            elseif option.flag==1
                % isochron line fixed at the origin
                rand_X3=rand_X';
                rand_Y3=rand_Y';
                rand_b=(rand_X3'*rand_X3)\rand_X3'*rand_Y3;
                rand_a=0;
            end
            % admit the random intercept and slope if they both fall within
            % their 1 sigma errors
            if rand_a<=a+sigma_a && rand_a>=a-sigma_a && rand_b<=b+sigma_b && rand_b>=b-sigma_b
                break;
            end
        end
        % define the upper and lower error lines
        Y=X*rand_b+rand_a;
        if isnan(Y3(1))
            Y3=Y;
            Y4=Y;
        else
            for k=1:length(X)
                if Y(k)>Y3(k)
                    Y3(k)=Y(k);
                end
                if Y(k)<Y4(k)
                    Y4(k)=Y(k);
                end
            end
        end
    end
    if option.flag==0   % rgb 204, 204, 204
        plot(X,Y3,'Color',[.8 .8 .8]);
        plot(X,Y4,'Color',[.8 .8 .8]);
    elseif option.flag==1   % rgb 204, 204, 255
        plot(X,Y3,'Color',[.8 .8 1]);
        plot(X,Y4,'Color',[.8 .8 1]);
    elseif option.flag==2   % rgb 255, 204, 204
        plot(X,Y3,'Color',[1 .8 .8]);
        plot(X,Y4,'Color',[1 .8 .8]);
    end

    % step 3: reset the range of axis x and y
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

    % step 4: plot errorbar
    if option.flag~=2
        error_ellipse(linearized_data,'k');
    else
        m=length(linearized_data.x);
        linearized_data.x(:,m)=[];
        linearized_data.dx(:,m)=[];
        linearized_data.y(:,m)=[];
        linearized_data.dy(:,m)=[];
        error_ellipse(linearized_data,'k');
        pb_data.x=data.x(:,m);
        pb_data.dx=data.dx(:,m);
        pb_data.y=data.y(:,m);
        pb_data.dy=data.dy(:,m);
        error_ellipse(pb_data,'g');
        data.x(:,m)=[];
        data.dx(:,m)=[];
        data.y(:,m)=[];
        data.dy(:,m)=[];
    end
    error_ellipse(data,[.5 .5 .5]);
    if removed_data.x(:,1)~=-1
        error_ellipse(removed_data,[1 0 0]);
    end
    % finish plotting
    hold off;

    % export vector paintings of the results
    if option.flag==0
        print(fig,'-painters','-dpdf','isochron');
    elseif option.flag==1
        print(fig,'-painters','-dpdf','lower limit');
    elseif option.flag==2
        print(fig,'-painters','-dpdf','upper limit');
    end
    
    % sub-function
    function error_ellipse(data,color)
        t=linspace(0,2*pi,1000);
        for l=1:length(data.x)
            x=data.x(l)+data.dx(l)*cos(t);
            y=data.y(l)+data.dy(l)*sin(t);
            if ischar(color)
                plot(x,y,color);
            else
                plot(x,y,'Color',color);
            end
        end
    end
end
