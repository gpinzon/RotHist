REFUGEE 	[ Rotational historiEs oF yoUnG stEllar objEcts ]
Tool for the analysis of global trends between rotation and stellar
parameters such as accretion, magnetic field strenghts, disc timelifes
and stellar winds.

README FILE

Files in this directory:

*Input : File "input.dat" containing :

Col[1] Disk timelife in Myr,

Col[2]  End of simulations in Myr,

Col[3] Characteristic timescale for accretion in yr,

Col[4] Magnetic field strenght in G

*Stellar Models (Siess+2000) in files :
mXX.zo02.var1, mXX.zo02.var2, mXX.zo02.hrd, mXX.zo02.xsurf

*other input files inside directory: (LEAVE THERE!!)
 labolo.txt  mcore.txt mradi.txt nat*.dat  KH95_BC_TEFF_VI.dat edad.dat


For running just type ./refugee.sh

Output : File "refugee.dat" 
containing 26 columns describing the rotational track  over the mass interval 0.1 - 7.0 Mo

Col[1] t in Myr 
Col[2] Angular frequency in solar units
Col[3] The accretion rate in Mo / yr
Col[4] Stellar radius (Siess+2000) in Ro
Col[5] Equatorial rotational velocity in km/s
Col[6] Stellar period in days
Col[7] Convective turnover time ( Interpolated from Landin +2010 )
Col[8] Trunctation radius (Ro)
Col[9] Corotation radius (Ro)
Col[10] Stellar wind torque (in cgs)
Col[11] Accretion torque (in cgs)
Col[12] Magnetic torque (in cgs)
Col[13] Escape velocity (in cgs)
Col[14] Inertia moment core (solar units)
Col[15] Total torque ( in cgs)
Col[16] Initial Mass (Mo)
Col[17] Stellar mass (Mo)
Col[18] Mass loss rate (Mo/yr)
Col[19] dR*/dt (computed from Siess+2010)
Col[20] dIc*/dt (computed from Siess+2010)
Col[21] Mass loss rate (in cgs)
Col[22] Stellar dipolar magnetic moment (cgs)
Col[23] ksquare
Col[24] Total torque ( incgs )
Col[25] Magnetic field strength (G)
Col[26] Alvén radius in solar units


The file "refugee.dat" 
contains only mass (Mo) vs equatorial velocity (km/s). This velocity is computed at TFINAL (second column input file)

If you wish to run just for one mass, you should use the source c file corresponding to desired mass which operates under initial conditions given by input.dat as well.

Any question does not hessitate to contact us.
gapinzone@unal.edu.co

2021
_________________________________________________
You can  run it on https://mybinder.org/
Do not forget change permissions before using : chmod 777 refugee.sh 




