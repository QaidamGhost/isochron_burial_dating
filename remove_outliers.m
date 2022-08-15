function [new_data,removed_data,MSWD] = remove_outliers(data,alpha,option)

%% Statistically remove the reworked clast(s) following Odom, 2020
% Note that this script might not findout the outliers whose values are
% extremly far away from the others. These outliers can easily exert their
% influence on the isochron line and occasionally escape from the
% criterion, while the MSWD remains high during iteration. Please try to
% remove these outliers manually at first if it helps dramatically reduce
% the MSWD (see Odom, 2020).
%
%% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; 1xn vector)
%   data.dx: 1 sigma absolute error of 10Be (atom/g; 1xn vector)
%   data.y: measured 26Al concentration (atom/g; 1xn vector)
%   data.dy: 1 sigma absolute error of 26Al (atom/g; 1xn vector)
% alpha: cutoff value for confidence intervals (unitless; scalar)
% option have fields as:
%   option.flag: "0" for default usage, "1" for york regression whose
%   intercept isfixed at zero, and "2" for protecting the added post-burial
%   data point from being removed (unitless; scalar)
%   option.Npb: maximum post-burial concentration and have fields as:
%       option.Npb.x: post-burial 10Be concentration (atom/g; scalar)
%       option.Npb.dx: 1 sigma absolute error of 10Be (atom/g; scalar)
%       option.Npb.y: post-burial 26Al concentration (atom/g; scalar)
%       option.Npb.dy: 1 sigma absolute error of 26Al (atom/g; scalar)

%% Output:
% new_data: ditto
% removed_data: ditto
% MSWD (unitless; scalar)

    % if no sample is removed
    removed_data.x(:,1)=-1;
    removed_data.dx(:,1)=-1;
    removed_data.y(:,1)=-1;
    removed_data.dy(:,1)=-1;
    MSWD=-1;
    
    if option.flag==0
        out = york(data,alpha);
    elseif option.flag==1
        out = york_fixed_intercept(data,alpha,0);
    % combine the post-burial concentration with the measured data
    % the post-burial concentration is the last elements in each field in
    % data2
    elseif option.flag==2
        data2.x=[data.x,option.Npb.x];
        data2.dx=[data.dx,option.Npb.dx];
        data2.y=[data.y,option.Npb.y];
        data2.dy=[data.dy,option.Npb.dy];
        out = york(data2,alpha);
    end
    if (out.mswd<=1)
        new_data=data;
        disp('No reworked clast is removed.');
        return;
    end
    
    count=0;
    while(1)
        b=out.b;
        a=out.a;
        if option.flag==0||option.flag==1
            X=data.x;
            sX=data.dx;
            Y=data.y;
            sY=data.dy;
        elseif option.flag==2
            X=data2.x;
            sX=data2.dx;
            Y=data2.y;
            sY=data2.dy;
        end
        w=1./(sY.*sY+b.*b.*sX.*sX);
        % weighted deviation for each sample (and the post-burial point)
        wd=w.*(Y-b.*X-a).^2;
        n=size(X,2);
        % degrees of freedom
        df=n-2;
        max_dev=tinv(1-alpha,df);

        [max_wd,i]=max(wd);
        %wd2=wd;
        %wd2(:,i)=[];
        %[max_wd2,i2]=max(wd2);
        
        if (max_wd<=max_dev)
            if option.flag==0||option.flag==1
                new_data=data;
            elseif option.flag==2
                m=size(data2.x,2);
                % the post-burial concentration is protected from being
                % removed
                % export data without the added post-burial concentration
                data2.x(:,m)=[];
                data2.dx(:,m)=[];
                data2.y(:,m)=[];
                data2.dy(:,m)=[];
                new_data=data2;
            end
            if (count==0)
                disp('No reworked clast is removed.');
            else
                if (count==1)
                    disp('1 reworked clast is removed.');
                else
                    fprintf('%d reworked clasts are removed.\n',count);
                end
            end
            if out.mswd>1
                fprintf('Warning! MSWD is beyond "1".\n');
            end
            return;
        end
        if (b*X(i)+a-Y(i)<=0) % the outlier is above the isochron line
            if option.flag==0||option.flag==1
                new_data=data;
            elseif option.flag==2
                m=size(data2.x,2);
                % the post-burial concentration is protected from being
                % removed
                % export data without the added post-burial concentration
                data2.x(:,m)=[];
                data2.dx(:,m)=[];
                data2.y(:,m)=[];
                data2.dy(:,m)=[];
                new_data=data2;
            end
            if (count==0)
                disp('No reworked clast is removed.');
            else
                if (count==1)
                    disp('1 reworked clast is removed.');
                else
                    fprintf('%d reworked clasts are removed.\n',count);
                end
            end
            if out.mswd>1
                fprintf('Warning! MSWD is beyond "1".\n');
            end
            return;
        end
        % protect the post-burial concentration data point 
        if option.flag==2
            m=size(data2.x,2);
            if m==i
                disp('Warning! The post-burial concentration data point has the highest weighted deviation.');
                data2.x(:,m)=[];
                data2.dx(:,m)=[];
                data2.y(:,m)=[];
                data2.dy(:,m)=[];
                new_data=data2;
                if out.mswd>1
                    fprintf('Warning! MSWD is beyond "1".\n');
                end
                return;
            end
            
        end
        % remove the outlier from data or data2
        count=count+1;
        if option.flag==0||option.flag==1
            removed_data.x(:,count)=data.x(:,i);
            removed_data.dx(:,count)=data.dx(:,i);
            removed_data.y(:,count)=data.y(:,i);
            removed_data.dy(:,count)=data.dy(:,i);
            data.x(:,i)=[];
            data.dx(:,i)=[];
            data.y(:,i)=[];
            data.dy(:,i)=[];
        elseif option.flag==2
            removed_data.x(:,count)=data2.x(:,i);
            removed_data.dx(:,count)=data2.dx(:,i);
            removed_data.y(:,count)=data2.y(:,i);
            removed_data.dy(:,count)=data2.dy(:,i);
            data2.x(:,i)=[];
            data2.dx(:,i)=[];
            data2.y(:,i)=[];
            data2.dy(:,i)=[];
        end

        if option.flag==0
            out = york(data,alpha);
        elseif option.flag==1
            out = york_fixed_intercept(data,alpha,0);
        elseif option.flag==2
            out = york(data2,alpha);
        end

        MSWD=out.mswd;
        if(out.mswd<=1)
            if (count==1)
                disp('1 reworked clast is removed.');
            else
                fprintf('%d reworked clasts are removed.\n',count);
            end
            if option.flag==0||option.flag==1
                new_data=data;
            elseif option.flag==2
                m=size(data2.x,2);
                % the post-burial concentration is protected from being
                % removed
                % export data without the added post-burial concentration
                data2.x(:,m)=[];
                data2.dx(:,m)=[];
                data2.y(:,m)=[];
                data2.dy(:,m)=[];
                new_data=data2;
            end
            return;
        end
    end
end
