function [outputArg1,outputArg2] = followContour(array,bounds)
%FOLLOWCONTOUR follows the contour of an array region of value 1 defining 
% by bounds.

% Constants
[ROWS,COLS,P] = getConstantsLW;

% Definition
start_pt = bounds(1,:);
end_pt = bounds(2,:);

% Sweep down and up in parallel
pt_contour = []; % Contains all the countour
array_current = array;

dir = [1,-1];
i = 1;
for j = 1:2
    dx = dir(j);
    dy = 1;
    pt = start_pt; % Current point
    while ~isequal(pt,end_pt)
        % If not hitting the edges
        if pt(1)>1 && pt(1)<ROWS && pt(2)>0 && pt(2)<COLS
            if array(pt(1)+dx,pt(2))==1 && ~ismember([pt(1)+dx,pt(2)],pt_contour,'rows') % up
                pt = [pt(1)+dx,pt(2)]; 
            elseif array(pt(1)+dx,pt(2)+dy)==1 && ~ismember([pt(1)+dx,pt(2)+dy],pt_contour,'rows') % diagonal up
                pt = [pt(1)+dx,pt(2)+dy];
            elseif array(pt(1),pt(2)+dy)==1 % right
                pt = [pt(1),pt(2)+dy];
            elseif array(pt(1)-dx,pt(2)+dy)==1 % diagonal down
                pt = [pt(1)-dx,pt(2)+dy];
            elseif array(pt(1)-dx,pt(2))==1 && i>3 % down
                pt = [pt(1)-dx,pt(2)];
            elseif i>2 
                pt = pt_contour(i-2,:);
            end
        % If on horizontal edges
        elseif (pt(1)==1 || pt(1)==ROWS) && (pt(2)>1) && (pt(2)<COLS)
            if array(pt(1),pt(2)+dy)==1 % right
                pt = [pt(1),pt(2)+dy];
            elseif array(pt(1)-dx,pt(2)+dy)==1 % diagonal down
                pt = [pt(1)-dx,pt(2)+dy];
            elseif array(pt(1)-dx,pt(2))==1 % down
                pt = [pt(1)-dx,pt(2)];
            end
        % If on vertical edges
        elseif pt(2)==COLS
            if pt(1) == 1% Lower right corner
                pt = [pt(1)+dx,pt(2)]; 
            elseif pt(1) == ROWS && ~ismember([pt(1)-dx,pt(2)],pt_contour,'rows') % Upper right corner
                pt = [pt(1)-dx,pt(2)];
            elseif array(pt(1)+dx,pt(2))==1 && ~ismember([pt(1)+dx,pt(2)],pt_contour,'rows') % up
                pt = [pt(1)+dx,pt(2)]; 
            elseif array(pt(1)-dx,pt(2))==1 % down
                pt = [pt(1)-dx,pt(2)];
            end
        end
        pt_contour(i,:) = pt;
        % Debug
        array_current(pt(1),pt(2)) = max(max(array))+1;
        figclf(20);
        surf(array_current); view(2); axis equal;

        % Increment
        i = i+1;
    end
    disp('End of contour');
end
