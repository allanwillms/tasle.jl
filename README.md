# TASLE
TASLE (Top-down Algorithm Splitting at Largest Error) is software for bounding time series data with a near minimax, continuous, piecewise linear band.

Copyright 2009 Allan Willms and Emily Szusz.

This Octave/MATLAB software source code is distributed under the GNU General Public License, Version 3.

Bug reports and comments should be sent to Allan Willms. 

<h3>Description</h3>

TASLE constructs a continuous piecewise linear band (two piecewise linear curves differing by a constant vertical offset) which bounds a given data set {(t<sub>i</sub>, y<sub>i</sub>)}. The quality of the fit is governed by a parameter which determines the minimum length between opposite sign slope discontinuities. The algorithm scales linearly with the number of data points.

TASLE is an improvement on an earlier related algorithm called linenvelope.

The TASLE algorithm is described in
<ul>
    <li>E.K. Szusz, A.R. Willms, <cite>A linear time algorithm for near minimax continuous piecewise linear representations of discrete data</cite>, SIAM J. Sci. Comput. 32 (5) (2010) 2584-2602. 
</ul>
All of the data sets used in the above paper may be downloaded via one of the two formats below.
<ul>
    <li>datasets.tar.gz (45 KB)
    <li>datasets.zip (45 KB) 
</ul>
These data sets and their original sources are:
<ol>
    <li>Echocardiogram (ECG) data with measurements taken roughly every 0.008 seconds over a period of five seconds. Units: seconds and millivolts. Source: PhysioNet, the first five seconds from the record "learning-set/n01" of their "AF Termination Challenge Database."
    <li>The average monthly exchange rate between the American and Canadian dollar from 1979 through 2009. Units: days and US dollars. Source: The Pacific Exchange Rate Service, Sauder School of Business, University of British Columbia. Retrieval settings: base currency: Canadian dollars, target currency: U.S. dollars, start date: Jan. 1, 1979, end date: Dec. 31, 2009, volume notation.
    <li>The relative brightness of the variable star X Cygni recorded roughly once or twice per day over a five month period. Units: days and relative brightness. Source: The American Association of Variable Star Observers' (AAVSO). Retrieval settings: star name: "X Cyg" (AUID: 000-BCM-291), start date: 29 April 2009, stop date: 24 September 2009.
    <li>Voltage recordings from a pyloric dilator neuron of a spiny lobster, where the voltage in millivolts is measured regularly every 0.2 milliseconds for 2 seconds. Source: Jack Peck, Ithaca College. 
      </ol>
