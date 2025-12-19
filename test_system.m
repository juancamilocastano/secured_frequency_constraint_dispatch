%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                                                                  			 	  %%%%%
%%%%    IEEE PES Power Grid Library - Optimal Power Flow with HVDC Lines - v23.09     %%%%%
%%%%          (https://github.com/power-grid-lib/pglib-opf-hvdc)      				  %%%%%
%%%%               Benchmark Group - Typical Operations               				  %%%%%
%%%%                         23 - August - 2023                       				  %%%%%
%%%%                                                                  				  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   AC/DC grid OPF test case based on:
%   Sass, F., Sennewald, T., Marten, A.-K. and Westermann, D. (2017), 
%   Mixed AC high-voltage direct current benchmark test system for security constrained optimal power flow calculation. 
%   IET Gener. Transm. Distrib., 11: 447-455. https://doi.org/10.1049/iet-gtd.2016.0993
%
%   Copyright (c) 2023 by Hakan Ergun, Vaishally Bhardwaj ({hakan.ergun, vaishally.bhardwaj}@kuleuven.be)
%   Licensed under the Creative Commons Attribution 4.0
%   International license, http://creativecommons.org/licenses/by/4.0/
%

function mpc = case67
mpc.version = '2';
mpc.baseMVA = 100.0;
mpc.baseKG=31034.48;

%After Vm, they are generic values

%% bus data
%	bus_i	type		Pd		Qd		  Gs	 Bs	   area		Vm	   Va		   baseKV	zone	  Vmax			 Vmin
mpc.bus = [
	1	 	1	 		0.0	 	0.0	 	  0.0	 0.0	 1	    1.0	    0	 		380.0	  1	      1.10000	    0.90000;
	2	 	2	 		0.0	    0.0	      0.0	 0.0	 1	    1.0	    0	 		380.0	  1	      1.10000	    0.90000;
	3	 	3	 		220.0	0.0	 	  0.0	 0.0	 1	    1.0	    0	 		380.0	  1	      1.10000	    0.90000;
	4	 	1	 		0.0	 	0.0	 	  0.0	 0.0	 2	    1.0	    0	 		380.0	  2	      1.10000	    0.90000;
	5	 	2	 	  	0.0	    0.0	      0.0	 0.0	 2	    1.0526	0	 		380.0	  2	      1.10000	    0.90000;
	6	 	3	 		190.0	0.0	      0.0	 0.0	 2	    1.0	    0	 		380.0	  2	      1.10000	    0.90000;
    ];



%% generator data
%	bus	   Pg	  Qg		Qmax		 Qmin	  Vg	   mBase	status		Pmax	  Pmin	 Pc1	 Pc2   Qc1min	Qc1max    Qc2min	Qc2max	ramp_agc	    ramp_10	      ramp_30	  ramp_q	   apf         alpha       IC
mpc.gen = [
	1	  0.0	  0.0	 	 0.0	    	 0.0	 0.0	     100.0	 0	 		120.0	  12.0	 0.0	 0.0	 0.0	  0.0	    0.0	     0.0	  0.0	 			0.0	 		0.0	 		0.0	 		0.0			0.0		   4;
	2	  0.0	  0.0	 	 0.0	    	 0.0	 0.0	     100.0	 0	 		 50.0	   5.0	 0.0	 0.0	 0.0	  0.0	    0.0	     0.0	  0.0	 			0.0	 		0.0	 		0.0	 		0.0			0.0		   4;
	3	  0.0	  0.0	 	 0.0	    	 0.0	 0.0	     100.0	 0	 		 50.0	   5.0	 0.0	 0.0	 0.0	  0.0	    0.0	     0.0	  0.0	 			0.0	 		0.0	 		0.0	 		0.0			0.0		   4;
	4	  0.0	  0.0	 	 0.0	    	 0.0	 0.0	     100.0	 0	 		120.0	  12.0	 0.0	 0.0	 0.0	  0.0	    0.0	     0.0	  0.0	 			0.0	 		0.0	 		0.0	 		0.0			0.0		   4;
	5	  0.0	  0.0	 	 0.0	    	 0.0	 0.0	     100.0	 0	 		 50.0	   5.0	 0.0	 0.0	 0.0	  0.0	    0.0	     0.0	  0.0	 			0.0	 		0.0	 		0.0	 		0.0			0.0		   4;
	6	  0.0	  0.0	 	 0.0	    	 0.0	 0.0	     100.0	 0	 		 50.0	   5.0	 0.0	 0.0	 0.0	  0.0	    0.0	     0.0	  0.0	 			0.0	 		0.0	 		0.0	 		0.0			0.0		   4;
];


%% bus mapping: n1=1, n2=2, n3=3, n4=4, n5=5, n6=6

%% branch data r, x, b assigned randomly using random data for matpower case 14
%  1-fbus  2-tbus   3-r       4-x       5-b       6-rateA  7-rateB  8-rateC  9-ratio  10-angle  11-status  12-angmin  13-angmax
mpc.branch = [
     1       2     0.0235    0.0712    0.0356     100      100      100       0         0         1         -60         60;  % n1–n2
     1       3     0.0479    0.1925    0.0423     100      100      100       0         0         1         -60         60;  % n1–n3
     2       3     0.0526    0.1583    0.0317     100      100      100       0         0         1         -60         60;  % n2–n3
     4       5     0.0142    0.0453    0.0110     100      100      100       0         0         1         -60         60;  % n4–n5
     4       6     0.0668    0.1339    0.0125     100      100      100       0         0         1         -60         60;  % n4–n6
     5       6     0.1213    0.2661    0.0498     100      100      100       0         0         1         -60         60;  % n5–n6
];

%% dc grid topology
%colunm_names% dcpoles
mpc.dcpol=2;
% numbers of poles (1=monopolar grid, 2=bipolar grid)
%% bus data
%column_names%  busdc_i    grid    Pdc     Vdc   basekVdc    Vdcmax   Vdcmin   Cdc  area
mpc.busdc = [
    			1       1       0       1       500       1.1     0.9     0	  1;
    			2       1       0       1       500       1.1     0.9     0	  1;
				3       1       0       1       500       1.1     0.9     0	  2;
				4       1       0       1       500       1.1     0.9     0	  2;
];

%% converters
%column_names%  busdc_i busac_i   type_dc   type_ac   P_g   Q_g   islcc      Vtar   rtf 	xtf     transformer  tm      bf   filter   rc    xc     reactor   basekVac      Vmmax   Vmmin    Imax   status   LossA  LossB  LossCrec LossCinv    droop    Pdcset      Vdcset  dVdcset    Pacmax  Pacmin   Qacmax         Qacmin
mpc.convdc = [
                2       2           2       1        -577.5  0    0            1     0.01  0.01        1          1     0.01   0       0.01   0.01       0        500         1.1     0.9     1.1     1        1.103 0.887  2.885    2.885      0.0050    -465.9871   0.9999        0           20     20       20             -20;
                3       3           3       1         1000   0    0            1     0.01  0.01        1          1     0.01   0       0.01   0.01       0        500         1.1     0.9     1.1     1        1.103 0.887  2.885    2.885      0.0050     500.0000   0.9947        0           20     20       20             -20;
                4       4           3       1         -550   0    0            1     0.01  0.01        1          1     0.01   0       0.01   0.01       0        500         1.1     0.9     1.1     1        1.103 0.887  2.885    2.885      0.0050    -517.1051   1.0022        0           20     20       20             -20;
                1       6           3       1         -600   0    0            1     0.01  0.01        1          1     0.01   0       0.01   0.01       0        500         1.1     0.9     1.1     1        1.103 0.887  2.885    2.885      0.0050    -560.8855   0.9972        0           20     20       20             -20;
];


%% branches
%column_names%  fbusdc  tbusdc  r        l   c   rateA   rateB   rateC   status   DT
mpc.branchdc = [
      2        4       0.0012   0    0      20       20       20        1         0.2;
      3        1       0.0012   0    0      20       20       20        1         0.2;
];

%% generator cost data
% 1-model  2-startup  3-shutdown  4-n   5-c2   6-c1($/MWh)  7-c0
mpc.gencost = [
    2        0            0         3     0        0          0;   % G1
    2        0            0         3     0        0          0;   % G2
    2        0            0         3     0        0          0;   % G3
    2        0            0         3     0        25          0;   % G4
    2        0            0         3     0        35          0;   % G5
    2        0            0         3     0        45          0;   % G6
];

%% generator extra
%   bus     DT      IC      upramp [Mw/h]   downramp [Mw/h]   MUT [h]   MDT [h]
mpc.genextra = [
    1       15      4       720             720               4         4;
    2       15      4       720             720               4         4;
    3       15      4       720             720               4         4;
    4       15      4       720             720               4         4;
    5       15      4       720             720               4         4;
    6       15      4       720             720               4         4;
];


%% BESS data (custom)
%  1-AC_bus   2-area   3-Pmax_MW   4-Emax_MWh   5-DODmax   6-eta_c   7-eta_d   8-E_init_MWh   9-E_final_MWh   10-DeployTime_s   11-ReserveCost_$_per_MW_h
mpc.bess = [
    3           1        40         93.3        0.793      0.9       0.9        18.6            18.6             0.2                4.0;
    6           2        30         70.0        0.793      0.9       0.9        14.             1.0              0.2                4.0;
];




%% Electrolyzer data
%  1-AC_bus   2-area   3-Pmax_MW   4-Pmin_MW   5-eta_ely_MWh_per_kg   6-LoadFactor   7-max_hydrog_flow_kg   8-Smax_kg    9-Smin_kg   10-S0_kg   11-ST_kg   12-DeployTime_s   13-ReserveCost_£_per_MW_h   14-eta_c   15-eta_d   16-StartupCost_£   17-CompressEnergy_MWh_per_kg   18-Hydrogen_market_price_Eur_per_kg
mpc.elec = [
    3   1   75   7.5   0.058   0.8   1293.10   31034.48   310.34   7758.62   7758.62   0.2   4   0.9   0.9   855.5   0.00167   0;
    6   2   65   6.5   0.058   0.8   1120.68   26896.50   268.96   6724.13   6724.13   0.2   4   0.9   0.9   855.5   0.00167   0;
];

%% Wind data (reduced to 3 columns)
% 1-AC_bus   2-area   3-Pmax_MW
mpc.wind = [
    3   1   100;
    6   2   33.3;
];
