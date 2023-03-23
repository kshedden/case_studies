# Argo program -- ocean temperature and salinity

The Argo program uses autonomous remote floats to measure temperature
and salinity throughout the world's oceans.  Each float is dropped
into the ocean by a ship, then drifts for several years collecting
data which are periodically sent via satellite to land.  The floats
submerge to several thousand meters, obtaining a high resolution
vertical profile of temperature, pressure and salinity.  Pressure here
is used as a proxy for depth, so the main interest is in the
relationships between temperature and depth, salinity and depth,
and the covariation of temperature and salinity.  For each
recorded profile, we know the latitude, longitude, and date where
the float was located.

You do not need to dive deeply into the underlying science, but if you
are curious the following two links may be informative:

https://www.frontiersin.org/articles/10.3389/fmars.2020.00700/full

https://argo.ucsd.edu/science/argo-and-climate-change

To obtain the data run one of the "get_data" scripts, then run
one of the "prep" scripts.  You will need to adjust the paths
in these scripts.

*About the measurement units:* Pressure is measured in decibars
(dbar), and 1 dbar very closely corresponds to 1 meter of depth.
Temperature is measured in centigrade units and salinity is
measured in "practical salinity units" (psu), which is the
grams of disolved ions per kilogram of water.  Typical seawater
has a psu of around 35.

