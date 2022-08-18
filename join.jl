function joinseg!(tjoin,segsort,loc,segment,t,y)
	# update tjoin values for the solution.  tjoin[i] is the intersection of the midlines
	# of segment i-1 and segment i.  (thus the left edge of segment i) 
	# tjoin values for segments loc to loc+2 must be recomputed, since loc and loc+1 are the new
	# segments resulting from the last split.
	# If n is the number of segments, then tjoin[1] and tjoin[n+1] are always t[1] and t[n], 
	# respectively, hence they never need updating.
	pt = Array{Tuple{Float64,Float64}}(undef,2)
	for i=max(2,loc):min(length(tjoin)-1,loc + 2)
		# find the intersection of the midline of segment i-1 and segment i
		(tjoin[i],yint) = calc_seg_intersection!(pt,segment[segsort[i-1:i]],t,y)
	end
	return nothing
end

function calc_seg_intersection!(pt,segments,t,y)
	# calculate the intersection of the midlines of two 
	# given segments, pt is an array of two tuples of 2 floats, used for storage
	for j=1:2
		ptloc = segments[j].first - 1 + segments[j].pivot[1]
		pt[j] = (t[ptloc], y[ptloc] - segments[j].pivottype[1]*segments[j].height*0.5)
	end
	return lineintersect(pt[1],segments[1].slope,pt[2],segments[2].slope)
end

function lineintersect(pt1,slope1,pt2,slope2)
	# find the intersection of two lines, each specified by a point (x,y) and a slope y'.
	# returns (x_int, y_int)
	return ((pt2[2] - pt1[2] - slope2*pt2[1] + slope1*pt1[1])/(slope1 - slope2),
			(pt2[2]*slope1 - pt1[2]*slope2 + slope1*slope2*(pt1[1] - pt2[1]))/(slope1 - slope2))
end

function computetmin(seglen,segsort,segment)
	# compute the tmin value for this set of tjoin points.  The tmin value is the smallest
	# segment length (as defined by seglen) that has a slope above or below both 
	# its neighbours.  Segment i is swallowed if seglen[i] < 0.  If a swallowed segment
	# has a slope above or below both its neighbours, then tmin will be <= 0
	# Returns (tmin, tminloc), where tminloc is the location of the constraining segment.
	tmin = Inf 
	tminloc = 0
	for i=2:length(segsort)-1
		if seglen[i] < tmin && !xor(segment[segsort[i]].slope > segment[segsort[i-1]].slope,
				 segment[segsort[i]].slope > segment[segsort[i+1]].slope)
			tmin = seglen[i]
			tminloc = i
		end
	end
	return (tmin, tminloc)
end
