# Argo program -- ocean temperature and salinity

The Argo program uses autonomous remote floats to measure temperature
and salinity throughout the world's oceans.  Each float is dropped
into the ocean by a ship, then drifts for several years collecting
data which are periodically sent via satellite to land.  The floats
periodically dive to several thousand meters, obtaining a high resolution
vertical profile of temperature, pressure and salinity for each dive.

The overall ARGO data is a very large collection.  We restrict our analyses
to the Pacific ocean in 2020, but this is still a big dataset and some of
the scripts will take a while to run.  You may need to randomly subsample
the data if your computer is slow.

Pressure here is used as a proxy for depth, so the main interest is in the
relationships between temperature and depth, salinity and depth,
and the covariation of temperature and salinity as depth varies.  For each
recorded profile, we know the latitude, longitude, and date where
the float was located.

The raw data are not "gridded", meaning that each dive has data for a
different set of pressures.  Also, the least and greatest pressure varies
from dive to dive.  To address this, we interpolate the data onto a uniform grid of
pressures and exclude all dives that do not have data covering this grid.
We exclude the data from the shallowest depths (near the surface) since
it is complicated to relate depth and pressure in this region due to waves.
You should review the prep.py or prep.jl scripts to see how this processing
is done in detail.

You do not need to dive deeply into the underlying science, but if you
are curious the following two links may be informative:

https://www.frontiersin.org/articles/10.3389/fmars.2020.00700/full

https://argo.ucsd.edu/science/argo-and-climate-change

To obtain the data run one of the "get_data" scripts, then run
one of the "prep" scripts (there is no prep.R so you will have to
use Python or Julia for this step).  As always for these case studies,
you will need to adjust the paths in the scripts.

*About the measurement units:* Pressure is measured in decibars
(dbar), and 1 dbar very closely corresponds to 1 meter of depth.
Temperature is measured in centigrade units and salinity is
measured in "practical salinity units" (psu), which is the
grams of disolved ions per kilogram of water.  Typical seawater
has a psu of around 35.
