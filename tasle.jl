module Tasle

export tasle

# Structure to store segment information.
# A Segment has a first and last index into the t and y arrays indicating the first and last
# points in this segment.  It also stores the segment's slope and height, plus an array
# of three indices for the pivot points.  These indices are relative to the segment, so if the 
# first pivot is 1, this refers to the point t[segment.first].  Finally there is an array of 
# three integers that specify the type for each of the pivots.  0 = undefined, 1 = upper, and
# -1 = lower (meaning the pivot is on the lower boundary of the band).
struct Segment
	first::Int64  # index into t,y arrays of first point in this segment
	last::Int64   # index into t,y arrays of last point in this segment
	splittable::Bool
	slope::Float64 
	height::Float64
	pivot::Array{Int64,1}  # array of size 3, indices relative to this segment
	pivottype::Array{Int64,1}  # array of size 3, 0=unknown, 1=upper, -1=lower
end

# Structure to store solution information.
# A Soln has segsort, and tjoin arrays, and the height and tmin values.
# segsort keeps track of the segments sorted by time, segment[segsort[i]] is the i'th 
# segment from the left for this solution.
# tjoin is an array (with nseg+1 entries, where nseg is the number of segments in this
# solution) that stores the join points.  
# tjoin[i] is the intersection of the midlines of the (i-1)'st and i'th segments.
# By definition the join point for the first segment is t[1], and tjoin[nseg+1] is t[n].
struct Soln
	segsort::Array{Int64,1}
	tjoin::Array{Float64}
	height::Float64
	tmin::Float64
end

include("sandwich.jl")
include("split.jl")
include("join.jl")

@enum Mode reduction recovery smoothing

"""
	(LB,UB,height,tmin) = tasle(t,y,user_tmin,[minheightfrac=0.8,variableheight=false])

Compute a piecewise linear bound on the data (`t`,`y`).

Tasle stands for Top-down Algorithm Splitting at the Largest Error.  It produces a near
minimax continuous piecewise linear band enclosing the data.  No segment of this band
whose slope does not lie between its neighbour's slopes will have midline length less than
`user_tmin`.  The returned value `tmin` is the minimum such length in the actual solution
and the returned value `height` is the maximal height of the band.  The band's lower bound
is defined by connecting the points in the mx2 matrix `LB` (first column is the t values,
second column is the y values) and the upper bound by the nx2 matrix `UB`.   The band may
have constant height, or each linear segment may have its own height (if
`varableheight=true`).  In the case that a constant height band is requested, the t values
in both `LB` and 'UB' will be identical.  The variable `minheightfrac` controls how much
splitting of segments is done after the largest height segment that cannot be split
without causing a violation of the constraints is found.  Segments whose height is below
`minheightfrac*height` will be split.  In the case that `variableheight=false`, all
segments are expanded vertically to have the same height as the constraining segment, but
the result is a smoother looking band.

# Examples #
```julia-repl
julia> t = range(0.0,10.0,50); y = @. sin(t) + 0.2*cos(7*t);
julia> (LB,UB,height,tmin) = tasle(t,y,0.4);
julia> height
0.3776184236760187
julia> p=plot(t,y,seriescolor=:blue,markersize=2,seriestype=:scatter)
julia> plot!(p,LB[:,1],LB[:,2],seriescolor=:red)
julia> plot!(p,UB[:,1],UB[:,2],seriescolor=:red)
```

---

The algorithm is described in
Emily K. Szusz, Allan R. Willms, 2010 "A linear time algorithm for near minimax
continuous piecewise linear representations of discrete data", SIAM J. Sci.
Comput. 32 (5), pp. 2584-2602, doi=10.1137/090769077.

This implementation has made a few adjustments to the algorithm compared to the 
description in that paper.  In particular, segments whose slopes are between their
neighbours' slopes are allowed to be swallowed.  Also, if a violation of the 
constraints is encountered, the algorithm focuses on splitting nearby segments to
see if the violation can be reversed.  If not the minimun height is considered to
be the height of the band whose splitting caused the violation.  The algorithm may
then continue to split less tall bands, up to `minheightfrac` of the minimum height.

# Copyright 2022 (Julia version) Allan R. Willms #

tasle is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

Dr. Allan R. Willms,
Dept. of Mathematics and Statistics,
University of Guelph,
Guelph, ON N1G 2W1,
Canada
"""
function tasle(t,y,user_tmin,minheightfrac=0.8,variableheight=false)

	# Number of data points.
	n = length(t)
	if length(y) != n
		error("t and y must be equal length")
	end
	if user_tmin < 0.0
		error("user_tmin must be nonnegative")
	end

	# Determine the minimum distance between data points and compare with user_tmin.
	if n < 3 || user_tmin <= minimum(t[2:n] .- t[1:n-1])
		# optimal solution is simply to connect the data
		return ([t y],[t y],0.0,tmincalc0(t,y))
	end

	# maxseg is the maximum number of segments allowed in a solution
	maxseg = n - 1
	# allseg keeps track of the total number of elements in the segment array.  This is
	# in general larger than the total number of segments since each time a segment is 
	# split, the old one is kept in the array and two new ones are added to the end.
	allseg = 1
	# recoverlevelmax is the maximum number of times we try splitting neighbouring
	# segments that are not helpful when working in recovery mode
	recoverlevelmax = 4
	########### Initialization of the first segment and pivot points ##########
	#
	# The segments are structures containing the absolute indices of the first 
	# and last data points within the segment, the slope of the band around the 
	# segment, the height of the segment, an array of 3 pivot point locations 
	# (relative to the segment) and an array of 3 pivottypes.

	# The first, middle, and last points are the initial pivots.  Their types are 
	# unknown=0.
	segmentarraysize = 2*n
	segment = Array{Segment}(undef,segmentarraysize)
	segment[1] = Segment(1,n,true,0.0,0.0,[1,div(n+1,2),n],zeros(Int64,3))
	# Process the first segment to minimize its height and determine its slope
	segment[1] = sandwich!(segment[1],t,y)

	# Initialize the array of intermediate solution structures.
	solnarraysize = 20
	soln = Array{Soln}(undef,solnarraysize)
	soln[1] = Soln([1],[t[1],t[n]],segment[1].height,Inf)

	# There are three Modes:
	# reduction : the current solution has tmin >= user_tmin, split the tallest
	#	  segment that is not marked as unsplittable, and has positive height, until no
	#	  such segments remain.  If recovery mode has been entered once and failed, then
	#	  the min height of the solution has been found.  In this case split any segment
	#	  that is below minheightfrac*soln.height, or would have split lengths (in terms 
	#	  of number of data points) long enough.
	# recovery : the current solution has tmin < user_tmin, split segments near
	#     the violating segment to try to alter this situation, if unsuccessful, mark the
	#     violating segment as unsplittable, return to reduction mode
	# smoothing : go through all segments splitting those that are around corners.
	mode = reduction
	current = 1  # index into soln array of the current solution being worked on
	seg2split = 1
	recoverlevel = 0
	recover = []
	violatingsegment = 0
	minheightfound = false
	while seg2split != 0
		# allocate more segment space if needed
		if allseg > segmentarraysize - 2
			segmentarraysize += n
			resize!(segment,segmentarraysize)
		end
		# split the selected segment of the current solution
		segment[allseg .+ [1,2]] = splitseg(segment[soln[current].segsort[seg2split]])
		# sandwich both new segments
		for i=1:2
			segment[allseg + i] = sandwich!(segment[allseg + i],t,y)
		end
		allseg += 2
		# construct the segsort and tjoin arrays for the new solution
		newsegsort = [soln[current].segsort[1:seg2split-1]; allseg-1; allseg; 
					  soln[current].segsort[seg2split+1:end]]
		newtjoin = [soln[current].tjoin[1:seg2split]; NaN; soln[current].tjoin[seg2split+1:end]]
		# compute the join points for seg2split, seg2split+1, and seg2split+2 
		# (other join points are not altered)
		joinseg!(newtjoin,newsegsort,seg2split,segment,t,y)
		# compute segment lengths
		seglen = @view(newtjoin[2:end]) .- @view(newtjoin[1:end-1])
		# compute tmin for this new solution
		(tmin, tminloc) = computetmin(seglen,newsegsort,segment)
		# find the height of the new solution 
		height = maximum(s -> s.height, segment[newsegsort])

		if tmin >= user_tmin
			# This is the best solution so far, replace the solution array with just this
			# solution.
			soln[1] = Soln(newsegsort,newtjoin,height,tmin)
			current = 1
		else
			# This solution is invalid but may yet be recovered, record to solution array.
			if current == solnarraysize  # make bigger if necessary
				solnarraysize += 20
				resize!(soln,solnarraysize)
			end
			current = current + 1
			soln[current] = Soln(newsegsort,newtjoin,height,tmin)
		end
		# Determine if we need to change modes and next segment to split.
		if mode != smoothing
			if tmin >= user_tmin
				mode = reduction
				recoverlevel = 0
				recover = []
				seg2split = findnextreductionsplit(soln[1],segment,minheightfrac,minheightfound)
				if seg2split == 0 || length(soln[current].segsort) == maxseg || height == 0.0
					mode = smoothing   # move to smoothing mode
					seg2split = findnextsmoothingsplit(soln[1],segment,2)
				end
			else
				if mode == reduction
					violatingsegment = seg2split  # record the violating segment
				end
				mode = recovery
				recoverlevel += 1
				if recoverlevel <= recoverlevelmax  # find next split
					segs = findnextrecoverysplit(soln[current],segment,tminloc)
					# If both entries of segs are not 0, then the first will be split
					# now and the second added to the recover array as a possible train.
					seg2split = segs[1]
					if segs[2] != 0  
						push!(recover,(recoverlevel,current,segs[2]))
					end
				end
				# If recoverlevel is above max or both entries of segs are 0, then 
				# nothing more can be done on this train.
				if recoverlevel > recoverlevelmax || segs[1] == 0 # give up on this train, pop another 
					if length(recover) > 0
						(recoverlevel,current,seg2split) = popfirst!(recover)
					else  # no recovery trains left, record violating segment as unsplittable
						mode = reduction
						current = 1
						recoverlevel = 0
						segment[soln[1].segsort[violatingsegment]] = 
								makeunsplittable(segment[soln[1].segsort[violatingsegment]])
						minheightfound = true
						seg2split = findnextreductionsplit(soln[1],segment,minheightfrac,minheightfound)
						if seg2split == 0 || length(soln[1].segsort) == maxseg || soln[1].height == 0.0
							mode = smoothing   # move to smoothing mode
							seg2split = findnextsmoothingsplit(soln[1],segment,2)
						end
					end
				end
			end
		else   # mode == smoothing 
			# if the split was successful just find the next to split
			if tmin >= user_tmin
				seg2split = findnextsmoothingsplit(soln[1],segment,seg2split)
			else  # unsuccessful split, so ignore new solution and find the next segment to split
				current = 1
				seg2split = findnextsmoothingsplit(soln[1],segment,seg2split+1)
			end
		end
	end
	# Calculate (t,y) values for the lower and upper boundaries of the band
	(LB,UB) = calculateband(soln[1],segment,t,y,variableheight,minheightfrac)

	#plotsoln(user_tmin,LB,UB,soln[1],segment,t,y,true)
	return (LB,UB,soln[1].height,soln[1].tmin)
end

function makeunsplittable(seg)
	# replace the segment seg with an equivlant one that is marked unsplittable
	return Segment(seg.first,seg.last,false,seg.slope,seg.height,seg.pivot,seg.pivottype)
end

function calculateband(soln,segment,t,y,variableheight,minheightfrac)
	# calculate the lower bound, LB, and upper bound, UB, of a band defined by this solution.
	# We cannot simply use the tjoin points already computed because some segments may be swallowed.
	# The first column of both LB and UB are the t values, and the second columns are the y values.
	# In the variable height case the t values will be different.
	LB = Array{Float64,2}(undef,length(soln.segsort) + 1,2)
	m = calculatebound!(LB,-1,soln,segment,t,y,variableheight,minheightfrac)
	LB = LB[1:m,:]
	UB = Array{Float64,2}(undef,length(soln.segsort) + 1,2)
	m = calculatebound!(UB,1,soln,segment,t,y,variableheight,minheightfrac)
	UB = UB[1:m,:]
	return (LB,UB)
end

function calculatebound!(bnd,boundtype,soln,segment,t,y,variableheight,minheightfrac)
	# calculate the boundary (lower or upper depending on whether boundtype is -1 or +1) of a band 
	# defined by this solution.
	# We cannot simply use the tjoin points already computed because some segments may be swallowed.
	# The first column of bnd is the t values, and the second column is the y values.
	if variableheight
		height = i -> max(segment[i].height,minheightfrac*soln.height)
	else
		height = i -> soln.height
	end
	allind = collect(1:length(soln.segsort))  # This will keep the unswalled segment indices.
	pt = Array{Tuple{Float64,Float64}}(undef,2)  # pre-allocated space for calc_seg_intersection
	i = 1
	while i <= length(allind)+1
		valid = false
		while !valid && i <= length(allind)+1
			# If i=1 or i=last+1 then we compute the end point otherwise we join segments
			if i==1 
				bnd[i,1:2] = calculate_endpoint(soln.segsort[allind[1]],segment,boundtype,height,t[1],t,y)
			elseif i==length(allind)+1
				bnd[i,1:2] = calculate_endpoint(soln.segsort[allind[end]],segment,boundtype,height,t[end],t,y)
			else
				# Join segment allind[i-1] to segment allind[i]
				inds = soln.segsort[allind[i-1:i]]
				(midt,midy) = calc_seg_intersection!(pt,segment[inds],t,y)
				# calc_seg_intersection finds the intersection of the midlines of the segments, adjust
				# to get the upper or lower boundary intersections
				bnd[i,1] = midt - boundtype*(height(inds[1]) - height(inds[2]))/
										 (2.0*(segment[inds[1]].slope - segment[inds[2]].slope))
				bnd[i,2] = midy + boundtype*height(inds[1])*0.5 + segment[inds[1]].slope*(bnd[i,1] - midt)
			end
			if i>1 && bnd[i,1] <= bnd[i-1,1]  # segment allind[i-1] is swallowed, ignore it
				i -= 1
				popat!(allind,i)
			else
				valid = true
			end
		end
		i += 1
	end
	return i-1
end

function calculate_endpoint(ind,segment,boundtype,height,tedge,t,y)
	loc = segment[ind].first - 1 + segment[ind].pivot[1]
	return [tedge, y[loc] - segment[ind].pivottype[1]*segment[ind].height*0.5 +
			boundtype*height(ind)*0.5 + segment[ind].slope*(tedge - t[loc])]
end

function tmincalc0(t,y)
	# compute tmin in the case that the optimal solution is simply to connect all data
	# points in a zero height band.
	slopecalc(i) = (y[i] - y[i-1])/(t[i] - t[i-1])
	slope = [0.0, slopecalc(3), slopecalc(2)]
	tmin = Inf
	for i=4:length(t)
		slope[1:3] = [slope[2:3];slopecalc(i)]
		if (slope[1] <= slope[2] && slope[2] <= slope[3]) || (slope[1] >= slope[2] && slope[2] >= slope[3])
			len = t[i-1] - t[i-2]
			if len < tmin
				tmin = len
			end
		end
	end
	return tmin
end

end
