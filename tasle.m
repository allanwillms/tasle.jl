function [Tout,Yout,height]=tasle(t,y,user_tmin)
% tasle  Top-down Algorithm Splitting at Largest Error.  Piecewise linear bounds on data.
%
%   [Tout,Yout,height] = tasle(t,y,user_tmin) will take the data (t,y) and produce 
%      a near minimax continuous piecewise linear band enclosing the data.  The band's 
%      lower bound is defined by the points (Tout,Yout) and the upper bound by 
%      (Tout,Yout+height).  No segment whose slope does not lie between its neighbour's 
%      slopes will have length less than user_tmin.

%  The algorithm is described in
%  Emily K. Szusz, Allan R. Willms, "A linear time algorithm for near minimax
%  continuous piecewise linear representations of discrete data", to appear in
%  SIAM J. Sci. Comput.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2010 Allan R. Willms
%
% tasle is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.
% 
% Allan Willms
% Dept. of Mathematics and Statistics
% University of Guelph
% Guelph, ON N1G 2W1
% Canada
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of data points.
n = length(t);
if length(y) ~= n
    error('t and y must be equal length');
end

% Determine the minimum distance between data points.
mindt = min(t(2:n) - t(1:n-1));
% Compare with user_tmin.
if user_tmin <= mindt
    % optimal solution is simply to connect the data
    Tout = t(:);
    Yout = y(:);
    height = 0;
    return;
end

%%%%%%%%%%% Initialization of the first segment and pivot points %%%%%%%%%%
%
% The segments are structures containing the absolute indices of the first 
% and last data points within the segment, the slope of the band around the 
% segment, the height of the segment, and a structure of 3 pivot points. 
% The pivot point types are:
%   -1 = pivot is on lower boundary
%    0 = unknown type
%    1 = pivot is on upper boundary
% The third pivot point is initially unknown until the segment is
% sandwiched.  The first and second pivots are initially set to LOWER and UPPER
% respectively, although the other way around would work too.
LOWER = -1;
UNKNOWN = 0;
UPPER = 1;

% The first and last points are the initial pivots.
segment(1) = struct('first',1,'last',n,'slope',[],'height',[],...
    'pivot',struct('loc',{1,n,0},'type',{UPPER,LOWER,UNKNOWN}));
% Process the first segment to minimize its height and determine its slope
segment(1) = remez(segment(1),t,y);

%%%%%%%%%%% Initialization of set of solutions %%%%%%%%%%%%%%%%%%%%%%%%%
% Structure to store tsort, tjoin arrays and the tmin value for intermediate 
% solutions.
% tsort keeps track of the segments sorted by time, tsort(i) is the i'th 
% segment from the left
% tjoin is an array (with nseg+1 entries) that stores the join points.  
% tjoin(i) is the join between the (i-1)'st and i'th segments from the left.  
% By definition the join point for the first segment is the first data t value, and
% the tjoin(nseg+1) is t(n).
soln(1) = struct('tsort',1,'tjoin',[t(1);t(n)],'height',segment(1).height,'tmin',inf);
%%%%%%%%%%%%%%%%%%%%%%%% End of initializations %%%%%%%%%%%%%%%%%%%%%%%%

% nseg keeps track of the number of segments in a solution
nseg = 1;
% maxseg is the maximum number of segments allowed
maxseg = segment(1).last - 1;
% allseg keeps track of the total number of elements in the segment array.  This is
% in general larger than nseg since each time a segment is split, the old one is
% kept in the array and two new ones are added to the end.
allseg = 1;

% The tmin array holds constraint information.
% There are two constraints:
%  1. tjoin values are strictly increasing (no segment is swallowed by its
%     neighbours.
%  2. If a segment's slope is above or below both its neighbours' slopes then the
%     segment must be no shorter than user_tmin.
% If segment(soln(nseg).tsort(i+1)).tjoin - segment(soln(nseg).tsort(i)).tjoin is 
% nonpositive then tmin(i) is set to mindt. (A holding value which will ensure this
% intermediate solution is not used but allows the main loop to continue in the
% event that the segment gets "unswallowed" due to subsequent segment splitting.)
% If the i'th segment from the left violates constraint 2 then tmin(i) holds the
% value: segment(soln(nseg).tsort(i+1)).tjoin - segment(soln(nseg).tsort(i)).tjoin
% If both constraints are satisfied, then tmin(i) holds the value inf.
tmin(1) = inf;
eff_tmin = inf;

% Since we only have one segment, it contains the biggest height
height = segment(1).height;
current = 1;

% Main while loop 
% In mode=1 (reduction mode), the data is processed, always using the largest height
% segment, until eff_tmin is less than user_tmin and the sum of the lengths of the
% worst violator and its two neighbors is also less than user_tmin.
% At this point, the set of intermediate solutions is checked to see which one
% obeying the two constraints has the smallest height.
% Mode then becomes 2 (smoothing mode) and segments of the candidate soln are checked 
% starting from the left to determine which can be split to reduce changes in the slope 
% without swallowing a segment.  Only segments which, due to their geometry, will
% never cause a violation of constraint 2 are split in this mode.
mode = 1;
% tmin_stop = max(0.9*mindt,0.1*user_tmin);
while (mode == 1 || (mode == 2 && current < nseg))

    % If in smoothing mode simply increase current by one if the slopes of the
    % neighbours of the current segment are not correctly related depending on
    % the type of pivot(2), or if the segment's height is already zero.
    if mode == 2 && (segment(soln(nseg).tsort(current)).height < eps || ...
	((segment(soln(nseg).tsort(current)).pivot(2).type == LOWER) && ...
	(segment(soln(nseg).tsort(current-1)).slope > segment(soln(nseg).tsort(current)).slope || ...
         segment(soln(nseg).tsort(current+1)).slope < segment(soln(nseg).tsort(current)).slope)) || ...
	((segment(soln(nseg).tsort(current)).pivot(2).type == UPPER) && ...
        (segment(soln(nseg).tsort(current-1)).slope < segment(soln(nseg).tsort(current)).slope || ...
         segment(soln(nseg).tsort(current+1)).slope > segment(soln(nseg).tsort(current)).slope)))
	current = current + 1;
	continue;  % return to top of while loop.
    end
    % Increase counters for splitting and initialize the next soln.
    nseg = nseg + 1;
    allseg = allseg + 2;
    soln(nseg).tsort = [soln(nseg-1).tsort(1:current-1),allseg-1,allseg,...
	soln(nseg-1).tsort(current+1:nseg-1)];
    soln(nseg).tjoin = [soln(nseg-1).tjoin(1:current);soln(nseg-1).tjoin(current:nseg)];
    
    % Split the current segment and sandwich both parts to minimize their heights.
    [segment(soln(nseg).tsort(current)), segment(soln(nseg).tsort(current+1))] = ...
	split_segment(segment(soln(nseg-1).tsort(current)));
    segment(soln(nseg).tsort(current)) = remez(segment(soln(nseg).tsort(current)),t,y);
    segment(soln(nseg).tsort(current+1)) = remez(segment(soln(nseg).tsort(current+1)),t,y);
    
    % newtjoin keeps track of which segments need to have their tjoin value
    % calculated.  Segment 1 always has tjoin=t(1) and tjoin(nseg+1) is always t(n),
    % so they never need recalculating.
    newtjoin = max(2,current):min(nseg,current+2);
    
    % Calculate new tjoin values.
    for i = newtjoin
        % Expand the segment to largest height and join with segment to the left to
        % produce a piecewise smooth band.
        % The equations for the bottom line of segment 1 and 2 respectively are:
        %  Y = Y1 - (H + type1*H1)/2 + M1*(T-T1)
        %  Y = Y2 - (H + type2*H2)/2 + M2*(T-T2)
        % where the known values for each line [(T1,Y1) and (T2,Y2)] are the data at
        % pivot point 1 for each segment, H1 and H2 are the heights of the bands,
        % H is the largest height of all the bands, and type is the pivot type (-1 if
        % lower, 1 if upper).
        % Equating the Y values and solving for T we get 
        % T = (Y2 - Y1 + M1*T1 - M2*T2 + type1*H1/2 - type2*H2/2)/(M1 - M2), 
        % (The pivotting process guarantees that M1 is not equal to M2.)
        soln(nseg).tjoin(i) = (y(segment(soln(nseg).tsort(i)).pivot(1).loc) - ...
            y(segment(soln(nseg).tsort(i-1)).pivot(1).loc) + ...
            segment(soln(nseg).tsort(i-1)).slope*t(segment(soln(nseg).tsort(i-1)).pivot(1).loc) - ...
            segment(soln(nseg).tsort(i)).slope*t(segment(soln(nseg).tsort(i)).pivot(1).loc) + ...
            segment(soln(nseg).tsort(i-1)).pivot(1).type*segment(soln(nseg).tsort(i-1)).height*0.5 - ...
            segment(soln(nseg).tsort(i)).pivot(1).type*segment(soln(nseg).tsort(i)).height*0.5)/...
            (segment(soln(nseg).tsort(i-1)).slope - segment(soln(nseg).tsort(i)).slope);
    end
    
    
    
    % newtmin keeps track of which segments will have had their widths altered due
    % to changes in tjoin values.
    newtmin = max(1,current-1):min(nseg,current+2);

    % Shift old values of tmin over one spot
    if mode == 2
	tminprevious = tmin;
    end
    undo = 0;
    tmin(current+2:nseg) = tmin(current+1:nseg-1);

    for i=newtmin
	delta = soln(nseg).tjoin(i+1) - soln(nseg).tjoin(i);
	if (delta <= 0)
	    if mode == 1
		tmin(i) = mindt;
	    else % mode == 2
		undo = 1;
	    end
	elseif (i>1 && i<nseg) && ...
	    ((segment(soln(nseg).tsort(i)).slope > segment(soln(nseg).tsort(i-1)).slope && ...
	      segment(soln(nseg).tsort(i)).slope > segment(soln(nseg).tsort(i+1)).slope) || ...
	     (segment(soln(nseg).tsort(i)).slope < segment(soln(nseg).tsort(i-1)).slope && ...
	      segment(soln(nseg).tsort(i)).slope < segment(soln(nseg).tsort(i+1)).slope))
	    if mode == 1
		tmin(i) = delta;
	    else % mode == 2
		undo = 1;
	    end
	else
	    tmin(i) = inf;
	end
	if undo % mode == 2
	    % this split caused a violation so undo it and move one segment along
	    nseg = nseg - 1;
	    tmin = tminprevious;
	    current = current + 1;
	    break
	end
    end
    
    if mode == 1
	% Effective tmin value this pass through the while loop is the minimum of 
	% all tmin values
	[eff_tmin,loc] = min(tmin);
	% Find the segment with the biggest height for the next round.
	[height,current] = find_biggest(segment,nseg,soln(nseg).tsort);
	soln(nseg).height = height;
	soln(nseg).tmin = eff_tmin;
    
	if (height < eps || (eff_tmin < user_tmin && soln(nseg).tjoin(loc+2) - ...
	    soln(nseg).tjoin(loc-1) < user_tmin))
	    % Determine the best soln found that satisfies the user_tmin value.
	    while soln(nseg).tmin < user_tmin
		nseg = nseg-1;
	    end
	    mode = 2;
	    current = 2;
	end
    end
end

[height,current] = find_biggest(segment,nseg,soln(nseg).tsort);
% Calculate the output T and Y values
Tout = soln(nseg).tjoin;
Yout = zeros(nseg+1,1);
for i=1:nseg
    Yout(i) = y(segment(soln(nseg).tsort(i)).pivot(1).loc) - 0.5*(height + ...
	segment(soln(nseg).tsort(i)).pivot(1).type*segment(soln(nseg).tsort(i)).height) + ...
	segment(soln(nseg).tsort(i)).slope*...
	(Tout(i) - t(segment(soln(nseg).tsort(i)).pivot(1).loc));
end
Yout(nseg+1) = Yout(nseg) + segment(soln(nseg).tsort(nseg)).slope*(Tout(nseg+1) - Tout(nseg));
