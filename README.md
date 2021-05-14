# Active space model [CAS(113e,76o)] for the FeMoco in nitrogenase

For details, see the following reference:

"The electronic complexity of the ground-state of the FeMo cofactor of nitrogenase as relevant to quantum simulations"
Zhendong Li, Junhao Li, Nikesh S. Dattani, C. J. Umrigar, and Garnet Kin-Lic Chan, J. Chem. Phys. 150, 024302 (2019).
doi: https://aip.scitation.org/doi/10.1063/1.5063376
arXiv: https://arxiv.org/abs/1809.10307

Some explanations about how the FCIDUMP file is generated:

The DFT calculations and integral transformation were performed with the PYSCF package.

Geometry:

Fe 10.018000 -7.568000 56.109001 \
Fe 12.528000 -8.326000 55.692001 \
Fe 11.277000 -6.902000 53.853001 \
Fe 11.990000 -5.771000 56.166000 \
Fe 14.820000 -7.601000 54.671001 \
Fe 13.550000 -6.198000 52.848999 \
Fe 14.226000 -5.072000 55.094002 \
Mo 16.132000 -5.638000 53.306000 \
C 13.044000 -6.665000 54.723999 \
S 11.719000 -7.459000 57.646000 \
S 16.169001 -6.043000 55.618000 \
S 10.001000 -5.522000 55.134998 \
S 13.162000 -3.960000 56.627998 \
S 14.407000 -9.544000 55.750999 \
S 14.337000 -4.109000 53.098000 \
S 10.718000 -8.998000 54.466999 \
S 15.203000 -7.626000 52.453999 \
S 11.674000 -6.465000 51.712002 \
S 8.197000 -8.169000 57.320999 \
C 13.298000 -2.024000 50.055000 \
N 12.366000 -2.944000 50.620998 \
C 11.087000 -2.665000 50.879002 \
N 10.465000 -1.648000 50.319000 \
N 10.437000 -3.395000 51.761002 \
C 15.000000 -1.023000 59.625000 \
C 16.169001 -0.765000 58.679001 \
C 17.408001 -1.557000 59.002998 \
O 17.870001 -1.613000 60.158001 \
N 18.021999 -2.150000 57.976002 \
C 8.265000 -1.347000 60.679001 \
C 9.399000 -1.996000 59.963001 \
N 9.933000 -3.152000 60.439999 \
C 10.065000 -1.670000 58.792999 \
C 10.903000 -3.530000 59.608002 \
N 10.987000 -2.647000 58.599998 \
C 8.211000 -4.616000 49.907001 \
C 7.638000 -4.371000 51.291000 \
C 6.881000 -3.235000 51.567001 \
C 7.762000 -5.325000 52.312000 \
C 6.315000 -3.036000 52.811001 \
C 7.180000 -5.134000 53.570000 \
C 6.480000 -3.971000 53.827999 \
O 5.954000 -3.723000 55.074001 \
N 5.234000 -10.420000 55.743999 \
C 6.175000 -9.825000 56.689999 \
C 5.492000 -9.417000 57.998001 \
O 4.993000 -8.311000 58.167999 \
C 6.936000 -8.649000 56.117001 \
N 5.537000 -10.342000 58.952000 \
C 18.375999 -9.779000 57.895000 \
O 18.938999 -9.535000 58.952999 \
N 17.193001 -9.276000 57.611000 \
C 16.604000 -8.174000 58.390999 \
C 15.700000 -8.561000 59.548000 \
O 16.066999 -8.490000 60.719002 \
N 14.482000 -8.960000 59.175999 \
C 13.413000 -9.170000 60.129002 \
C 12.448000 -10.321000 59.819000 \
O 11.585000 -10.619000 60.638000 \
N 12.608000 -10.986000 58.669998 \
C 11.669000 -12.033000 58.242001 \
C 12.397000 -13.165000 57.548000 \
O 12.311000 -14.342000 57.938000 \
N 13.074000 -12.829000 56.459999 \
C 13.702000 -13.811000 55.596001 \
C 14.324000 -13.116000 54.386002 \
C 13.275000 -12.471000 53.491001 \
C 13.911000 -11.621000 52.422001 \
N 12.878000 -11.123000 51.516998 \
C 13.118000 -10.372000 50.460999 \
N 14.374000 -10.095000 50.109001 \
N 12.094000 -9.945000 49.705002 \
O 22.684000 -7.070000 47.084999 \
C 18.180000 -7.369000 50.759998 \
C 18.711000 -7.463000 52.162998 \
N 18.159000 -6.788000 53.228001 \
C 19.879999 -8.011000 52.627998 \
C 18.976999 -6.904000 54.275002 \
N 20.003000 -7.668000 53.919998 \
C 17.480000 -1.176000 54.216000 \
C 17.136999 -1.687000 52.838001 \
C 17.757999 -3.064000 52.515999 \
C 19.298000 -2.979000 52.311001 \
C 19.787001 -4.275000 51.637001 \
C 21.243001 -4.658999 51.755000 \
C 17.111000 -3.651000 51.271999 \
O 17.582001 -2.027000 55.126999 \
O 17.659000 0.041000 54.360001 \
O 21.753999 -5.267001 50.771999 \
O 21.844000 -4.408001 52.838001 \
O 16.628000 -4.810000 51.337002 \
O 17.166000 -2.989000 50.195999 \
O 17.490000 -3.979000 53.599998 \
O 15.888000 -4.907000 46.811001 \
O 19.811001 -2.775000 44.562000 \
O 21.295000 -1.675000 56.474998 \
O 14.797000 -5.542000 49.339001 \
O 20.999001 -7.558000 49.209999 \
O 21.018000 -4.234000 48.202000 \
O 18.273001 -3.799000 47.834000 \
O 24.264000 -5.111000 53.854000 \
O 18.688999 -0.721000 49.472000 \
O 20.888000 -4.190000 55.507000 \
O 19.867001 0.773000 55.806000 \
O 19.983999 2.139000 49.183998 \
O 19.601999 0.653000 51.790001 \
H 14.253909 -2.551439 49.902245 \
H 13.485333 -1.148637 50.707842 \
H 12.976245 -1.666277 49.058362 \
H 12.756517 -3.749057 51.133080 \
H 10.948382 -1.006864 49.696509 \
H 9.418611 -3.350350 51.840380 \
H 10.851817 -4.301103 52.064489 \
H 14.176828 -0.309351 59.447134 \
H 15.333467 -0.925316 60.673517 \
H 14.594556 -2.042045 59.487892 \
H 15.865440 -0.939269 57.629728 \
H 16.465566 0.302217 58.741955 \
H 17.747578 -2.051739 56.983900 \
H 18.886142 -2.654352 58.166876 \
H 7.943033 -0.415330 60.180924 \
H 8.542154 -1.100693 61.721569 \
H 7.393391 -2.026762 60.730390 \
H 11.631869 -2.750219 57.800192 \
H 9.950464 -0.842408 58.089483 \
H 11.523221 -4.427460 59.671513 \
H 8.474673 -3.670967 49.398976 \
H 7.479231 -5.140569 49.262567 \
H 9.119194 -5.242955 49.957698 \
H 6.731721 -2.478799 50.782532 \
H 8.369443 -6.224645 52.146371 \
H 5.741910 -2.119104 53.011911 \
H 7.324840 -5.870775 54.367803 \
H 5.546271 -2.840425 55.063551 \
H 5.750165 -10.676231 54.895864 \
H 4.564723 -9.698382 55.448024 \
H 6.922962 -10.606213 56.936088 \
H 6.237039 -7.808570 55.941648 \
H 7.401622 -8.936745 55.157491 \
H 6.037121 -11.218482 58.832310 \
H 5.164447 -10.123761 59.874931 \
H 18.781655 -10.466135 57.105626 \
H 16.704789 -9.500843 56.734897 \
H 16.043744 -7.544498 57.677036 \
H 17.429698 -7.586986 58.824307 \
H 14.219156 -8.808225 58.189655 \
H 12.783401 -8.263772 60.234252 \
H 13.872680 -9.367031 61.112602 \
H 13.256381 -10.608473 57.970429 \
H 11.174293 -12.420930 59.144023 \
H 10.895271 -11.603903 57.570696 \
H 13.174849 -11.839973 56.192474 \
H 14.473224 -14.370958 56.160710 \
H 12.945500 -14.559086 55.282907 \
H 15.032972 -12.341634 54.741837 \
H 14.913117 -13.854269 53.810118 \
H 12.641772 -13.254547 53.028674 \
H 12.605520 -11.821374 54.090680 \
H 14.436696 -10.768534 52.906517 \
H 14.648390 -12.199351 51.832233 \
H 11.938028 -11.063509 51.914966 \
H 14.531875 -9.489172 49.305241 \
H 15.069805 -9.986501 50.852478 \
H 12.200422 -9.064588 49.201927 \
H 23.601459 -6.857982 47.330021 \
H 22.249074 -6.194096 47.039440 \
H 18.955188 -7.698977 50.048433 \
H 17.905390 -6.332504 50.518737 \
H 17.264032 -7.973340 50.626379 \
H 20.628332 -8.602807 52.103837 \
H 18.853442 -6.425326 55.245736 \
H 16.033435 -1.787606 52.792646 \
H 17.426546 -0.972948 52.052632 \
H 19.528063 -2.102232 51.681236 \
H 19.759202 -2.868484 53.309073 \
H 19.247123 -5.115528 52.106350 \
H 19.502098 -4.271284 50.571493 \
H 16.780810 -4.624223 47.121885 \
H 15.433213 -5.151276 47.644086 \
H 19.980766 -3.723341 44.435067 \
H 19.321145 -2.766364 45.406260 \
H 21.109337 -1.750623 57.425786 \
H 20.759935 -0.882264 56.201074 \
H 15.501236 -5.269415 49.977894 \
H 14.144043 -5.981686 49.917945 \
H 21.703906 -7.489930 48.526249 \
H 21.149790 -6.782293 49.795777 \
H 21.232006 -4.616723 49.099190 \
H 21.241122 -3.291069 48.291377 \  
H 19.195384 -4.089918 48.024181 \
H 17.876980 -3.527801 48.695535 \
H 24.697333 -5.483950 53.066151 \
H 23.382051 -4.837434 53.479899 \
H 18.953931 -0.321612 50.327834 \
H 18.115647 -1.481067 49.723219 \
H 20.992458 -3.241281 55.753084 \
H 21.160164 -4.239411 54.564410 \
H 20.330760 1.130055 55.030296 \
H 18.979135 0.500481 55.434638 \
H 19.940635 1.936802 50.139503 \
H 19.518351 1.363575 48.815374 \
H 20.354371 0.130718 52.119534 \
H 18.942524 0.622904 52.521728 \
H 20.831763 -7.852887 54.530690 \
H 9.545331 -1.367341 50.646139 \
H 11.144726 -10.137694 50.018567 \
C 23.292999 -8.105000 51.974998 \
C 24.523001 -7.559000 51.272999 \
O 25.587000 -8.169000 51.252998 \
N 24.365000 -6.351000 50.728001 \
C 25.365000 -5.756000 49.866001 \
C 24.761000 -4.521000 49.146000 \
O 22.260000 -7.893000 55.609001 \
H 17.491148 -3.347688 54.445403 \
H 22.481544 -8.230248 51.233996 \
H 22.943402 -7.372156 52.723626 \
H 23.534473 -9.061675 52.460204 \
H 23.421240 -5.901816 50.778903 \
H 26.261472 -5.456441 50.446366 \
H 25.720615 -6.508632 49.132117 \
H 24.470248 -3.748936 49.879990 \
H 23.848896 -4.798400 48.587383 \
H 25.488867 -4.081178 48.443117 \
H 21.937137 -7.439342 56.409331 \
H 22.855842 -7.224751 55.210211

Basis:
{'Fe':'TZP-DKH','Mo':'TZP-DKH','S':'TZP-DKH',\
'C':'def2-SVP','H':'def2-SVP','O':'def2-SVP','N':'def2-SVP'}

Charge = -3, Spin = 35 (PYSCF convention: 2*S) for UKS/B3LYP with sf-X2C Hamiltonian and COSMO eps=4.0

Molecular orbitals: Split-Localized MOs, see MOs_plots (ordered by Genetic Algorithm implemented in the BLOCK code)

Active space information: CAS(113e,76o) 
