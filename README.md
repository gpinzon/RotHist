REFUGEE 	[ Rotational historiEs oF yoUnG stEllar objEcts ]
Tool for the analysis of global trends between rotation and stellar parameters such as accretion, magnetic field strenghts, disc timelifes and stellar winds.

Input : File "input.dat" containing :

Col[1] Disk timelife in Myr,
Col[2]  End of simulations in Myr,
Col[3] Characteristic timescale for accretion in yr,
Col[4] Magnetic field strenght in G

For running just type ./refugee.sh
Output : File "OUT.dat" containing 26 columns describing the rotational track (mass-vsini) over the mass interval 0.1 - 7.0 Mo

The file "refugee.dat" contain only mass (mo) vs equatorial velocity (km/s). This velocity is computed at TFINAL (second column input file)

If you want to run just for one mass, you should use the source c file corresponding to desired mass which operates under initial conditions given by input.dat as well.



