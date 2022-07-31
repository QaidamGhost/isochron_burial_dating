function [new_data,removed_data,MSWD] = kill_outliers_Odom(data,alpha)

% Statistically remove reworked clasts following Odom, 2020
%
% Note that this script might not properly findout the outliers whose 
% values are extremly far away from the others. These outliers can easily 
% exert their influence on the isochron line and occasionally escape from 
% the criterion, while the MSWD remains high during iteration. Please 
% try to remove these outliers manually at first if it helps dramatically 
% reduce the MSWD (see Odom, 2020).
%
% References: 
% William III, E. Odom. Dating the Cenozoic incision history of the 
%   Tennessee and Shenandoah Rivers with cosmogenic nuclides and 40Ar/39Ar 
%   in manganese oxides. Diss. Purdue University Graduate School, 2020. 
%   (https://hammer.purdue.edu/articles/thesis/Dating_the_Cenozoic_incision
%   _history_of_the_Tennessee_and_Shenandoah_Rivers_with_cosmogenic
%   _nuclides_and_40Ar_39Ar_in_manganese_oxides/13275017
%
% Arguments:
% data must have fields as:
%   data.x: measured 10Be concentration (atom/g; nx1 vector)
%   data.dx: absolute error of 10Be (atom/g; nx1 vector)
%   data.y: measured 26Al concentration (atom/g; nx1 vector)
%   data.dy: absolute error of 26Al (atom/g; nx1 vector)
%
% Output:
% new_data: ditto
% removed_data: ditto
% MSWD (unitless; scalar)

    % if no sample is removed
    removed_data.x(:,1)=-1;
    removed_data.dx(:,1)=-1;
    removed_data.y(:,1)=-1;
    removed_data.dy(:,1)=-1;
    MSWD=-1;

    if (nargin<2)
        alpha=0.05;
    end

    out = york_Vermeesch(data,alpha);
    if (out.mswd<=1)
        new_data=data;
        disp('No reworked clast is precluded.\n');
        return;
    end
    
    count=0;
    while(1)
        b=out.b;
        a=out.a;
        X=data.x;
        sX=data.dx;
        Y=data.y;
        sY=data.dy;
        w=1./(sY.*sY+b.*b.*sX.*sX);
        % weighted deviation for each sample (York, 2004
        wd=w.*(Y-b.*X-a).^2;
        n=size(X,2);
        % degrees of freedom
        df=n-2;
        max_dev=tinv(1-alpha,df);

        [max_wd,i]=max(wd);
        if (max_wd<=max_dev)
            new_data=data;
            if (count==0)
                disp('No reworked clast is precluded.\n');
            else
                if (count==1)
                    disp('1 reworked clast is precluded.\n');
                else
                    fprintf('%d reworked clasts are precluded.\n',count);
                end
            end
            return;
        end
        if (b*X(i)+a-Y(i)<=0) % the outlier is above the isochron line
            new_data=data;
            if (count==0)
                disp('No reworked clast is precluded.\n');
            else
                if (count==1)
                    disp('1 reworked clast is precluded.\n');
                else
                    fprintf('%d reworked clasts are precluded.\n',count);
                end
            end
            fprintf('MSWD is beyond "1" and might still be high.\n');
            return;
        end
        % remove the outlier
        count=count+1;
        removed_data.x(:,count)=data.x(:,i);
        removed_data.dx(:,count)=data.dx(:,i);
        removed_data.y(:,count)=data.y(:,i);
        removed_data.dy(:,count)=data.dy(:,i);
        data.x(:,i)=[];
        data.dx(:,i)=[];
        data.y(:,i)=[];
        data.dy(:,i)=[];

        out = york_Vermeesch(data,alpha);
        MSWD=out.mswd;
        if(out.mswd<=1)
            if (count==1)
                disp('1 reworked clast is precluded.\n');
            else
                fprintf('%d reworked clasts are precluded.\n',count);
            end
            new_data=data;
            return;
        end
    end
end
