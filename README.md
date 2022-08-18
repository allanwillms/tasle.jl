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

