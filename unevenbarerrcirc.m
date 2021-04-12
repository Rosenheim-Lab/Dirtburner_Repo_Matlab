function [RR]=unevenbarerrcirc_beta(xmark,y,base,err,figaspect,EdgeColor,FaceColor)
%
%This function plots rectangles of unequal width on the active axis of the
%active graph.  The output is an "uneven" bar graph.
%
%INPUT VARIABLES
%
%xmark [n] is a vector of x values that mark the vertical edges of the
%rectangles.  Naturally, the vector will contain n elements for n-1
%rectangles, because 2 sides characterize each rectangle and all but 2 of
%the edges are shared.
%
%y [n-1] is the height of the rectangles.  Values must be positive
%integers.
%
%base [scalar] is the y value of the base of all rectangles.
%
%err [n-1] is a vector containing an uncertainty for each y (rectangle
%height).  If entered, the rectangles will have a line through the top of
%the rectangle depicting the uncertainty error bars.  If not entered, no
%line will be present.
%
%figaspect (scalar) is the aspect ratio of the figure in which these
%rectangles will be plotted.  It is calculated in what ever routine calls
%this one, by using get(handle,'Position') and dividing the 3rd element by
%the 4th element of the answer.  This ratio is w/h - y values must be
%multiplied by it to equate to x values and vice versa.
%
%EdgeColor [1X3] is the color of the edges of the rectangles where all
%elements of the vector are between 0 and 1
%
%FaceColor [1X3] is the color of the fill of the rectangles where all
%elements of the vector are between 0 and 1
%
%Defaults:  Black Edge Color, no FaceColor, base of 0 and no error bars.


if nargin==7;fc=1; end
if nargin<7; fc=0; end
if nargin<6; fc=0; EdgeColor=[0 0 0]; figaspect=1; end
if nargin<5; fc=0; EdgeColor=[0 0 0]; figaspect=1; err=zeros(length(y)); end
if nargin<4; fc=0; EdgeColor=[0 0 0]; figaspect=1; err=zeros(length(y));base=0; end
if numel(base)==1; base=base*ones(size(y)); end

for h=2:(length(xmark))
    w(h)=xmark(h)-xmark(h-1);       %find rectangle widths
    mp(h)=(xmark(h-1)+xmark(h))/2;  %find rectangle midpoints
end




for i=2:(length(xmark))
    if fc==1                         %If face color is expressed
        if (y(i-1)-base(i-1))>=0     %If base is positive or zero
            R(i)=rectangle('Position',[xmark(i-1) base(i-1) (xmark(i))-xmark(i-1) (roundn(y(i-1)-base(i-1),-3))],'EdgeColor',EdgeColor,'FaceColor',FaceColor);
            C(i)=rectangle('Position',[mp(i)-0.5*0.02*max(xlim) y(i-1)-0.5*0.02*max(ylim)*figaspect 0.02*max(xlim) 0.02*max(ylim)*figaspect], 'Curvature',[1 1]);  
        else                         %if base is negative
            R(i)=rectangle('Position',[xmark(i-1) base(i-1) (xmark(i))-xmark(i-1) (roundn(y(i-1)-base(i-1),-3))],'EdgeColor',EdgeColor,'FaceColor',FaceColor);
            C(i)=rectangle('Position',[mp(i)-0.5*0.02*max(xlim) base(i-1)-0.5*0.02*max(ylim)*figaspect 0.02*max(xlim) 0.02*max(ylim)*figaspect], 'Curvature',[1 1]);  
        end
    else                             %If face color is not expressed
        if (y(i-1)-base(i-1))>=0     %if base is positive or zero
            R(i)=rectangle('Position',[xmark(i-1) base(i-1) (xmark(i))-xmark(i-1) (roundn(y(i-1)-base(i-1),-3))],'EdgeColor',EdgeColor);
            C(i)=rectangle('Position',[mp(i)-0.5*0.02*max(xlim) y(i-1)-0.5*0.02*max(ylim)*figaspect 0.02*max(xlim) 0.02*max(ylim)*figaspect], 'Curvature',[1 1]);  
        else                         %if base is negative
            R(i)=rectangle('Position',[xmark(i-1) base(i-1) (xmark(i))-xmark(i-1) (roundn(y(i-1)-base(i-1),-3))],'EdgeColor',EdgeColor);
            C(i)=rectangle('Position',[mp(i)-0.5*0.02*max(xlim) base(i-1)-0.5*0.02*max(ylim)*figaspect 0.02*max(xlim) 0.02*max(ylim)*figaspect], 'Curvature',[1 1]);  
        end
    end
    eb(i)=line([mp(i) mp(i)],[(base(i-1)+y(i-1)-err(i-1)) (base(i-1)+y(i-1)+err(i-1))],'Color',EdgeColor);
    topwing(i)=line([mp(i)-0.5*0.02*max(xlim) mp(i)+0.5*0.02*max(xlim)],[(base(i-1)+y(i-1)+err(i-1)) (base(i-1)+y(i-1)+err(i-1))],'Color',EdgeColor);
    bottomwing(i)=line([mp(i)-0.5*0.02*max(xlim) mp(i)+0.5*0.02*max(xlim)],[(base(i-1)+y(i-1)-err(i-1)) (base(i-1)+y(i-1)-err(i-1))],'Color',EdgeColor);
end    

