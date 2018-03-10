El Nino Recharge Oscillation
-------------------------------------------------------
Copyright(c) 2018 Lingqiang He

A code for a simple recharge oscillator with white noise. The time series is ensembled 10 times. 

Execution
-------------------------------------------------------
Compile with g++ and input the total time of simulation. Add "> filename.dat" to save the time series. 

It outputs columns of number of 
    
    1) ensemble, 2) shift of recharge index, 3) time for current ensemble and shift (40 min per iteration), 4) Eastward SST for current       ensemble and shift, 5) Eastward thermocline depth for current ensemble and shift, 6) combined noise for current ensemble and shift,       7) white noise for current ensemble and shift
