# Space weather and X-ray flux

The sun emits streams of charged particles in all directions called
the [solar wind](https://en.wikipedia.org/wiki/Solar_wind) -- see
[here](https://www.youtube.com/watch?v=GX5FbXX-hks)
for a visualization.  The intensity of the solar wind varies according
to an 11-year [solar
cycle](https://en.wikipedia.org/wiki/Solar_cycle), and is also
influenced by occasional [coronal mass
ejections](https://en.wikipedia.org/wiki/Coronal_mass_ejection) and by
[sun spots](https://en.wikipedia.org/wiki/Sunspot).

The Earth is surrounded by a magnetic field called the
[magnetosphere](https://en.wikipedia.org/wiki/Magnetosphere), which
extends to about 6-10 times the Earth's radius in the direction of the
sun, and to about 60 times the Earth's radius in the direction
opposite the sun.  The magnetosphere largely deflects the solar wind,
and exists in a dynamic state reflecting what is known as [space
weather](https://en.wikipedia.org/wiki/Space_weather).  There is great
concern that major space weather events could cause extensive damage
and extended outages of power systems and other critical
infrastructure.

A series of satellites called
[GOES](https://en.wikipedia.org/wiki/Geostationary_Operational_Environmental_Satellite)
(for *Geostationary Operational Environmental Satellite*) have been
operated for several decades to monitor many characteristics of the
Earth, including the magnetosphere.  Here we analyze X-ray flux data
captured in two variables named flux-1 and flux-2.
These are time series measured at approximately a two-second sampling
rate.  However the sampling is not consistently at a two-second cadence,
and there are occasional large gaps.

X-ray flux has units of watts per square meter. A strong X-ray flux event is
called a "flare", and flares can be classified in increasing severity as
B, C, M, and X class flares according to these [thresholds](http://solar-center.stanford.edu/SID/activities/flare.html#:~:text=Scientists%20classify%20solar%20flares%20according,M9%2C%20and%20X1%20to%20X9.).

There are a number of interesting questions that can be addressed with
the GOES X-ray flux data.  Here are a few examples

* What is the temporal correlation structure of the flux levels?

* Are there extreme events in the flux levels, and if so how are
they distributed in time?

* Are there periodicities in the flux levels?

* What is the shape of the marginal flux distribution, and what can be
said about the right tail of the flux distribution (corresponding to the largest flux levels)?

* Can we predict flux levels in the near future from current flux
  levels?
