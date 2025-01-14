Sebastian's idea was to use the settling time or the max eigen value in continous time to find the sampling time t. This is possible,the problem is that in our case the system is simply stable and so the settling time is infinite (not defined if you want) while the max eig is in 0. This leads to an infinity sampling frequency and so MATLAB not being able to properly run the code.

My "solution" was to use the fft and find the armoninc content of the signal, with this I obtain the dominant frequency that I then used to find the minimun sampling time necessary to avoid aliasing (according to Shannon-Nyquist theorem).

h is around 2,8:2 seconds depending on the initial conditions, in particular the less symmetric they're the more needy (faster) is the sampling time. To get a good signal recontruction I used 1/10 of the minimum h and to avoid any problem related to possible increasing of the dominant frequency I set a condtion in which if h_min > 0.1 then use 0.1, otherwhise use the h_min.

Using h = 0.1 should allow for a very good signal reconstruction for dominant frequency up to 0,5 Hz (h = 1/(20*f_dominant)) while the frequency content even in very asymmetric cases should be confined in 0,4 (I,ve seen a maximum of 0,39 for x0= [-1000,0,1,0,-1,0]).

I've checked the open loop signal for a symmetric case and the frequency of the signal seems to be around 0,23/0,24 Hz that is consitent with the dominant one and the period of around 4,1 s that I found in the fft analysis.
Also you should see that as the initial conditions become asymmetric a new faster signal is superimposed to the old one (that is consitent with the faster frequencies' picks found in the fft analysis), this can be expalied physically with the fact that the small masses oscillates by themself other than the global system oscillation.

Try to check if what I've done is reasonable and if something is wrong do not esitate to correct it.