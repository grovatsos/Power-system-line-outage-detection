function [bus line] = buslinedata(num)
% 3 bus system
% Type 1= slack, type 2=PV, type 3= PQ
%         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi |  Qmin | Qmax |
bus3   = [ 1      1      1      0     1.9    0     0     0      0      0;
           2      3    1.01    0.01    0     0     1     0.4     -2      2;
           3      3    0.99    0.01    0     0    0.9    0.6     -2      2];
type = bus3(:,2);
bus3 = [bus3(:,1), bus3(:,3:8), zeros(size(bus3,1),2), type, bus3(:,10) bus3(:,9)];

        
%         |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%         |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |
line3     =[ 1      2         0       0.0504       0         1
             1      3         0       0.0636       0         1
             2      3         0       0.0372       0         1];
line3(:,5) = line3(:,5)*2;
line3 = [line3 zeros(size(line3,1),1)];



% IEEE 14 bus system

%         |Bus | Type | Vsp | theta | PGi | QGi | PLi | QLi |  Qmin | Qmax |
bus14   = [ 1     1    1.060   0     232.4 -16.9     0     0       0       0;
            2     2    1.045   -4.98   40   42.4 21.7   12.7    -40     50;
            3     2    1.010   -12.72  0    23.4 94.2   19.0     0      40;
            4     3    1.019   -10.33  0     0   47.8   -3.9     0       0;
            5     3    1.02    -8.78   0     0    7.6    1.6     0       0;
            6     2    1.070   -14.22  0    12.2  11.2    7.5    -6      24;
            7     3    1.062   -13.37  0     0    0.0    0.0     0       0;
            8     2    1.090   -13.36  0    17.4  0.0    0.0    -6      24;
            9     3    1.056   -14.94  0     0   29.5   16.6     0       0;
            10    3    1.051   -15.10  0     0    9.0    5.8     0       0;
            11    3    1.057   -14.79  0     0    3.5    1.8     0       0;
            12    3    1.055   -15.07  0     0    6.1    1.6     0       0;
            13    3    1.050   -15.16  0     0   13.5    5.8     0       0;
            14    3    1.036   -16.04  0     0   14.9    5.0     0       0 ];
bus14(:,5:10) = bus14(:,5:10)/100;
type = bus14(:,2);
bus14 = [bus14(:,1), bus14(:,3:8), zeros(size(bus14,1),2), type, bus14(:,10) bus14(:,9)];

        
%         |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
%         |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |
line14     =[1      2       0.01938   0.05917    0.0264         1
             1      5       0.05403   0.22304    0.0246         1
             2      3       0.04699   0.19797    0.0219         1
             2      4       0.05811   0.17632    0.0170         1
             2      5       0.05695   0.17388    0.0173         1
             3      4       0.06701   0.17103    0.0064         1
             4      5       0.01335   0.04211    0.0            1
             4      7       0.0       0.20912    0.0        0.978
             4      9       0.0       0.55618    0.0        0.969
             5      6       0.0       0.25202    0.0        0.932
             6     11       0.09498   0.19890    0.0            1
             6     12       0.12291   0.25581    0.0            1
             6     13       0.06615   0.13027    0.0            1
             7      8       0.0       0.17615    0.0            1
             7      9       0.0       0.11001    0.0            1
             9     10       0.03181   0.08450    0.0            1
             9     14       0.12711   0.27038    0.0            1
            10     11       0.08205   0.19207    0.0            1
            12     13       0.22092   0.19988    0.0            1
            13     14       0.17093   0.34802    0.0            1 ];
line14(:,5) = line14(:,5)*2;
line14 = [line14 zeros(size(line14,1),1)];

% IEEE 30 bus system

bus30 = [   1     1    1.06     0       0     0     0     0       0       0;
            2     2    1.043   -5.48   40   50.0  21.7   12.7    -40     50;
            3     3    1.021   -7.96    0     0    2.4    1.2     0       0;
            4     3    1.012   -9.62    0     0    7.6    1.6     0       0;
            5     2    1.01    -14.37   0   37.0  94.2   19.0    -40     40;
            6     3    1.01    -11.34   0     0    0.0    0.0     0       0;
            7     3    1.002   -13.12   0     0   22.8   10.9     0       0;
            8     2    1.01    -12.10   0   37.3  30.0   30.0    -10     40;
            9     3    1.051   -14.38   0     0    0.0    0.0     0       0;
            10    3    1.045   -15.97   0     0    5.8    2.0     0       0;
            11    2    1.082   -14.39   0   16.2   0.0    0.0    -6      24;
            12    3    1.057   -15.24   0     0   11.2    7.5     0       0;
            13    2    1.071   -15.24   0   10.6   0.0    0.0    -6      24;
            14    3    1.042   -16.13   0     0    6.2    1.6     0       0;
            15    3    1.038   -16.22   0     0    8.2    2.5     0       0;
            16    3    1.045   -15.83   0     0    3.5    1.8     0       0;
            17    3    1.040   -16.14   0     0    9.0    5.8     0       0;
            18    3    1.028   -16.82   0     0    3.2    0.9     0       0;
            19    3    1.026   -17.00   0     0    9.5    3.4     0       0;
            20    3    1.030   -16.80   0     0    2.2    0.7     0       0;
            21    3    1.033   -16.42   0     0   17.5   11.2     0       0;
            22    3    1.033   -16.41   0     0    0.0    0.0     0       0;
            23    3    1.027   -16.61   0     0    3.2    1.6     0       0;
            24    3    1.021   -16.78   0     0    8.7    6.7     0       0;
            25    3    1.017   -16.35   0     0    0.0    0.0     0       0;
            26    3    1.000   -16.77   0     0    3.5    2.3     0       0;
            27    3    1.023   -15.82   0     0    0.0    0.0     0       0;
            28    3    1.007   -11.97   0     0    0.0    0.0     0       0;
            29    3    1.003   -17.06   0     0    2.4    0.9     0       0;
            30    3    0.992   -17.94   0     0   10.6    1.9     0       0 ];
bus30(:,5:10) = bus30(:,5:10)/100;
type = bus30(:,2);
bus30 = [bus30(:,1), bus30(:,3:8), zeros(size(bus30,1),2), type, bus30(:,9:10)];
% bus30 = [bus30(:,1), bus30(:,3:8), zeros(size(bus30,1),2), type];


line30 =    [1      2       0.0192    0.0575     0.0264         1   0
             1      3       0.0452    0.1652     0.0204         1   0
             2      4       0.0570    0.1737     0.0184         1   0
             3      4       0.0132    0.0379     0.0042         1   0
             2      5       0.0472    0.1983     0.0209         1   0
             2      6       0.0581    0.1763     0.0187         1   0
             4      6       0.0119    0.0414     0.0045         1   0
             5      7       0.0460    0.1160     0.0102         1   0
             6      7       0.0267    0.0820     0.0085         1   0
             6      8       0.0120    0.0420     0.0045         1   0
             6      9       0.0       0.2080     0.0        0.978   0
             6     10       0.0       0.5560     0.0        0.969   0
             9     11       0.0       0.2080     0.0            1   0
             9     10       0.0       0.1100     0.0            1   0
             4     12       0.0       0.2560     0.0        0.932   0
            12     13       0.0       0.1400     0.0            1   0
            12     14       0.1231    0.2559     0.0            1   0
            12     15       0.0662    0.1304     0.0            1   0
            12     16       0.0945    0.1987     0.0            1   0
            14     15       0.2210    0.1997     0.0            1   0
            16     17       0.0824    0.1923     0.0            1   0
            15     18       0.1073    0.2185     0.0            1   0
            18     19       0.0639    0.1292     0.0            1   0
            19     20       0.0340    0.0680     0.0            1   0
            10     20       0.0936    0.2090     0.0            1   0
            10     17       0.0324    0.0845     0.0            1   0
            10     21       0.0348    0.0749     0.0            1   0
            10     22       0.0727    0.1499     0.0            1   0
            21     23       0.0116    0.0236     0.0            1   0
% According to diagram, there's no line from 21 to 23, rather 21 to 22
%             21     22       0.0116    0.0236     0.0            1   0
            15     23       0.1000    0.2020     0.0            1   0
            22     24       0.1150    0.1790     0.0            1   0
            23     24       0.1320    0.2700     0.0            1   0
            24     25       0.1885    0.3292     0.0            1   0
            25     26       0.2544    0.3800     0.0            1   0
            25     27       0.1093    0.2087     0.0            1   0
            28     27       0.0       0.3960     0.0        0.968   0
            27     29       0.2198    0.4153     0.0            1   0
            27     30       0.3202    0.6027     0.0            1   0
            29     30       0.2399    0.4533     0.0            1   0
             8     28       0.0636    0.2000     0.0214         1   0
             6     28       0.0169    0.0599     0.065          1   0 ];
line30(:,5) = line30(:,5)*2;
line30 = [line30 zeros(size(line30,1),1)];

% A 3-machine 9-bus system from Chow's book pp.70
% data3m9b.m
% modified for non-diagonalizable resonance
% bus data format
% bus: number, voltage(pu), angle(degree), p_gen(pu), q_gen(pu),
%      p_load(pu), q_load(pu),G shunt,B shunt, bus_type
%      bus_type - 1, swing bus
%               - 2, generator bus (PV bus)
%               - 3, load bus (PQ bus)

bus9 = [...
    1 1.04    0.00   0.716  0.27    0.00  0.00  0.00  0.00 1;
	2 1.025   9.3    1.63   0.067   0.00  0.00  0.00  0.00 2;
	3 1.025   4.7    0.85   -0.109  0.00  0.00  0.00  0.00 2;
	4 1.026   -2.2   0.00   0.00    0.00  0.00  0.00  0.00 3;
	5 0.996   -4.0   0.00   0.00    1.25  0.5   0.00  0.00 3;
	6 1.013   -3.7   0.00   0.00    0.90  0.3   0.00  0.00 3;
	7 1.026   3.7    0.00   0.00    0.00  0.00  0.00  0.00 3;
	8 1.016   0.7    0.00   0.00    1.00  0.35  0.00  0.00 3;
	9 1.032   2.0    0.00   0.00    0.00  0.00  0.00  0.00 3];

% line data format
% line: from bus, to bus, resistance(pu), reactance(pu),
%       line charging(pu), tap ratio

line9 = [...
    1 4 0.0    0.0576 0.     1. 0. ;
	4 5 0.01   0.085  0.176  1. 0. ;
	5 7 0.032  0.161  0.306  1. 0. ;
	4 6 0.017  0.092  0.158  1. 0. ;
	6 9 0.039  0.17   0.358  1. 0. ;
	7 8 0.0085 0.072  0.149  1. 0. ;
	3 9 0.0    0.0586 0.     1. 0. ;
	8 9 0.0119 0.1008 0.209  1. 0. ;
	2 7 0.00   0.0625 0.     1. 0. ];


%% bus data 118, note the TYPE from MATPOWER and the version we use are different 1 is PQ bus in Matpower
% Also, in our code, we always assumed bus 1 to be slack. Hence, we moved original bus 10 as slack to bus 1
% bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus118 = [
    1	1	0	0	0	0	1	1.05	35.61	345	1	1.06	0.94;
	2	3	20	9	0	0	1	0.971	11.22	138	1	1.06	0.94;
	3	3	39	10	0	0	1	0.968	11.56	138	1	1.06	0.94;
	4	2	39	12	0	0	1	0.998	15.28	138	1	1.06	0.94;
	5	3	0	0	0	-40	1	1.002	15.73	138	1	1.06	0.94;
	6	2	52	22	0	0	1	0.99	13	    138	1	1.06	0.94;
	7	3	19	2	0	0	1	0.989	12.56	138	1	1.06	0.94;
	8	2	28	0	0	0	1	1.015	20.77	345	1	1.06	0.94;
	9	3	0	0	0	0	1	1.043	28.02	345	1	1.06	0.94;
% 	10	2	0	0	0	0	1	1.05	35.61	345	1	1.06	0.94;
    10	2	51	27	0	0	1	0.955	10.67	138	1	1.06	0.94;
	11	3	70	23	0	0	1	0.985	12.72	138	1	1.06	0.94;
	12	2	47	10	0	0	1	0.99	12.2	138	1	1.06	0.94;
	13	3	34	16	0	0	1	0.968	11.35	138	1	1.06	0.94;
	14	3	14	1	0	0	1	0.984	11.5	138	1	1.06	0.94;
	15	2	90	30	0	0	1	0.97	11.23	138	1	1.06	0.94;
	16	3	25	10	0	0	1	0.984	11.91	138	1	1.06	0.94;
	17	3	11	3	0	0	1	0.995	13.74	138	1	1.06	0.94;
	18	2	60	34	0	0	1	0.973	11.53	138	1	1.06	0.94;
	19	2	45	25	0	0	1	0.963	11.05	138	1	1.06	0.94;
	20	3	18	3	0	0	1	0.958	11.93	138	1	1.06	0.94;
	21	3	14	8	0	0	1	0.959	13.52	138	1	1.06	0.94;
	22	3	10	5	0	0	1	0.97	16.08	138	1	1.06	0.94;
	23	3	7	3	0	0	1	1	21	138	1	1.06	0.94;
	24	2	13	0	0	0	1	0.992	20.89	138	1	1.06	0.94;
	25	2	0	0	0	0	1	1.05	27.93	138	1	1.06	0.94;
	26	2	0	0	0	0	1	1.015	29.71	345	1	1.06	0.94;
	27	2	71	13	0	0	1	0.968	15.35	138	1	1.06	0.94;
	28	3	17	7	0	0	1	0.962	13.62	138	1	1.06	0.94;
	29	3	24	4	0	0	1	0.963	12.63	138	1	1.06	0.94;
	30	3	0	0	0	0	1	0.968	18.79	345	1	1.06	0.94;
	31	2	43	27	0	0	1	0.967	12.75	138	1	1.06	0.94;
	32	2	59	23	0	0	1	0.964	14.8	138	1	1.06	0.94;
	33	3	23	9	0	0	1	0.972	10.63	138	1	1.06	0.94;
	34	2	59	26	0	14	1	0.986	11.3	138	1	1.06	0.94;
	35	3	33	9	0	0	1	0.981	10.87	138	1	1.06	0.94;
	36	2	31	17	0	0	1	0.98	10.87	138	1	1.06	0.94;
	37	3	0	0	0	-25	1	0.992	11.77	138	1	1.06	0.94;
	38	3	0	0	0	0	1	0.962	16.91	345	1	1.06	0.94;
	39	3	27	11	0	0	1	0.97	8.41	138	1	1.06	0.94;
	40	2	66	23	0	0	1	0.97	7.35	138	1	1.06	0.94;
	41	3	37	10	0	0	1	0.967	6.92	138	1	1.06	0.94;
	42	2	96	23	0	0	1	0.985	8.53	138	1	1.06	0.94;
	43	3	18	7	0	0	1	0.978	11.28	138	1	1.06	0.94;
	44	3	16	8	0	10	1	0.985	13.82	138	1	1.06	0.94;
	45	3	53	22	0	10	1	0.987	15.67	138	1	1.06	0.94;
	46	2	28	10	0	10	1	1.005	18.49	138	1	1.06	0.94;
	47	3	34	0	0	0	1	1.017	20.73	138	1	1.06	0.94;
	48	3	20	11	0	15	1	1.021	19.93	138	1	1.06	0.94;
	49	2	87	30	0	0	1	1.025	20.94	138	1	1.06	0.94;
	50	3	17	4	0	0	1	1.001	18.9	138	1	1.06	0.94;
	51	3	17	8	0	0	1	0.967	16.28	138	1	1.06	0.94;
	52	3	18	5	0	0	1	0.957	15.32	138	1	1.06	0.94;
	53	3	23	11	0	0	1	0.946	14.35	138	1	1.06	0.94;
	54	2	113	32	0	0	1	0.955	15.26	138	1	1.06	0.94;
	55	2	63	22	0	0	1	0.952	14.97	138	1	1.06	0.94;
	56	2	84	18	0	0	1	0.954	15.16	138	1	1.06	0.94;
	57	3	12	3	0	0	1	0.971	16.36	138	1	1.06	0.94;
	58	3	12	3	0	0	1	0.959	15.51	138	1	1.06	0.94;
	59	2	277	113	0	0	1	0.985	19.37	138	1	1.06	0.94;
	60	3	78	3	0	0	1	0.993	23.15	138	1	1.06	0.94;
	61	2	0	0	0	0	1	0.995	24.04	138	1	1.06	0.94;
	62	2	77	14	0	0	1	0.998	23.43	138	1	1.06	0.94;
	63	3	0	0	0	0	1	0.969	22.75	345	1	1.06	0.94;
	64	3	0	0	0	0	1	0.984	24.52	345	1	1.06	0.94;
	65	2	0	0	0	0	1	1.005	27.65	345	1	1.06	0.94;
	66	2	39	18	0	0	1	1.05	27.48	138	1	1.06	0.94;
	67	3	28	7	0	0	1	1.02	24.84	138	1	1.06	0.94;
	68	3	0	0	0	0	1	1.003	27.55	345	1	1.06	0.94;
% 	69	3	0	0	0	0	1	1.035	30	138	1	1.06	0.94;
    69	2	0	0	0	0	1	1.035	30	138	1	1.06	0.94;
	70	2	66	20	0	0	1	0.984	22.58	138	1	1.06	0.94;
	71	3	0	0	0	0	1	0.987	22.15	138	1	1.06	0.94;
	72	2	12	0	0	0	1	0.98	20.98	138	1	1.06	0.94;
	73	2	6	0	0	0	1	0.991	21.94	138	1	1.06	0.94;
	74	2	68	27	0	12	1	0.958	21.64	138	1	1.06	0.94;
	75	3	47	11	0	0	1	0.967	22.91	138	1	1.06	0.94;
	76	2	68	36	0	0	1	0.943	21.77	138	1	1.06	0.94;
	77	2	61	28	0	0	1	1.006	26.72	138	1	1.06	0.94;
	78	3	71	26	0	0	1	1.003	26.42	138	1	1.06	0.94;
	79	3	39	32	0	20	1	1.009	26.72	138	1	1.06	0.94;
	80	2	130	26	0	0	1	1.04	28.96	138	1	1.06	0.94;
	81	3	0	0	0	0	1	0.997	28.1	345	1	1.06	0.94;
	82	3	54	27	0	20	1	0.989	27.24	138	1	1.06	0.94;
	83	3	20	10	0	10	1	0.985	28.42	138	1	1.06	0.94;
	84	3	11	7	0	0	1	0.98	30.95	138	1	1.06	0.94;
	85	2	24	15	0	0	1	0.985	32.51	138	1	1.06	0.94;
	86	3	21	10	0	0	1	0.987	31.14	138	1	1.06	0.94;
	87	2	0	0	0	0	1	1.015	31.4	161	1	1.06	0.94;
	88	3	48	10	0	0	1	0.987	35.64	138	1	1.06	0.94;
	89	2	0	0	0	0	1	1.005	39.69	138	1	1.06	0.94;
	90	2	163	42	0	0	1	0.985	33.29	138	1	1.06	0.94;
	91	2	10	0	0	0	1	0.98	33.31	138	1	1.06	0.94;
	92	2	65	10	0	0	1	0.993	33.8	138	1	1.06	0.94;
	93	3	12	7	0	0	1	0.987	30.79	138	1	1.06	0.94;
	94	3	30	16	0	0	1	0.991	28.64	138	1	1.06	0.94;
	95	3	42	31	0	0	1	0.981	27.67	138	1	1.06	0.94;
	96	3	38	15	0	0	1	0.993	27.51	138	1	1.06	0.94;
	97	3	15	9	0	0	1	1.011	27.88	138	1	1.06	0.94;
	98	3	34	8	0	0	1	1.024	27.4	138	1	1.06	0.94;
	99	2	42	0	0	0	1	1.01	27.04	138	1	1.06	0.94;
	100	2	37	18	0	0	1	1.017	28.03	138	1	1.06	0.94;
	101	3	22	15	0	0	1	0.993	29.61	138	1	1.06	0.94;
	102	3	5	3	0	0	1	0.991	32.3	138	1	1.06	0.94;
	103	2	23	16	0	0	1	1.001	24.44	138	1	1.06	0.94;
	104	2	38	25	0	0	1	0.971	21.69	138	1	1.06	0.94;
	105	2	31	26	0	20	1	0.965	20.57	138	1	1.06	0.94;
	106	3	43	16	0	0	1	0.962	20.32	138	1	1.06	0.94;
	107	2	50	12	0	6	1	0.952	17.53	138	1	1.06	0.94;
	108	3	2	1	0	0	1	0.967	19.38	138	1	1.06	0.94;
	109	3	8	3	0	0	1	0.967	18.93	138	1	1.06	0.94;
	110	2	39	30	0	6	1	0.973	18.09	138	1	1.06	0.94;
	111	2	0	0	0	0	1	0.98	19.74	138	1	1.06	0.94;
	112	2	68	13	0	0	1	0.975	14.99	138	1	1.06	0.94;
	113	2	6	0	0	0	1	0.993	13.74	138	1	1.06	0.94;
	114	3	8	3	0	0	1	0.96	14.46	138	1	1.06	0.94;
	115	3	22	7	0	0	1	0.96	14.46	138	1	1.06	0.94;
	116	2	184	0	0	0	1	1.005	27.12	138	1	1.06	0.94;
	117	3	20	8	0	0	1	0.974	10.67	138	1	1.06	0.94;
	118	3	33	15	0	0	1	0.949	21.92	138	1	1.06	0.94;
];


%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
gen118 = [
	1	450	0	200	-147	1.05	100	1	550	0	0	0	0	0	0	0	0	0	0	0	0;
	4	0	0	300	-300	0.998	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	6	0	0	50	-13	0.99	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	8	0	0	300	-300	1.015	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
    10  0	0	15	-5	0.955	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	12	85	0	120	-35	0.99	100	1	185	0	0	0	0	0	0	0	0	0	0	0	0;
	15	0	0	30	-10	0.97	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	18	0	0	50	-16	0.973	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	19	0	0	24	-8	0.962	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	24	0	0	300	-300	0.992	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	25	220	0	140	-47	1.05	100	1	320	0	0	0	0	0	0	0	0	0	0	0	0;
	26	314	0	1000	-1000	1.015	100	1	414	0	0	0	0	0	0	0	0	0	0	0	0;
	27	0	0	300	-300	0.968	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	31	7	0	300	-300	0.967	100	1	107	0	0	0	0	0	0	0	0	0	0	0	0;
	32	0	0	42	-14	0.963	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	34	0	0	24	-8	0.984	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	36	0	0	24	-8	0.98	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	40	0	0	300	-300	0.97	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	42	0	0	300	-300	0.985	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	46	19	0	100	-100	1.005	100	1	119	0	0	0	0	0	0	0	0	0	0	0	0;
	49	204	0	210	-85	1.025	100	1	304	0	0	0	0	0	0	0	0	0	0	0	0;
	54	48	0	300	-300	0.955	100	1	148	0	0	0	0	0	0	0	0	0	0	0	0;
	55	0	0	23	-8	0.952	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	56	0	0	15	-8	0.954	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	59	155	0	180	-60	0.985	100	1	255	0	0	0	0	0	0	0	0	0	0	0	0;
	61	160	0	300	-100	0.995	100	1	260	0	0	0	0	0	0	0	0	0	0	0	0;
	62	0	0	20	-20	0.998	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	65	391	0	200	-67	1.005	100	1	491	0	0	0	0	0	0	0	0	0	0	0	0;
	66	392	0	200	-67	1.05	100	1	492	0	0	0	0	0	0	0	0	0	0	0	0;
	69	516.4	0	300	-300	1.035	100	1	805.2	0	0	0	0	0	0	0	0	0	0	0	0;
	70	0	0	32	-10	0.984	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	72	0	0	100	-100	0.98	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	73	0	0	100	-100	0.991	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	74	0	0	9	-6	0.958	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	76	0	0	23	-8	0.943	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	77	0	0	70	-20	1.006	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	80	477	0	280	-165	1.04	100	1	577	0	0	0	0	0	0	0	0	0	0	0	0;
	85	0	0	23	-8	0.985	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	87	4	0	1000	-100	1.015	100	1	104	0	0	0	0	0	0	0	0	0	0	0	0;
	89	607	0	300	-210	1.005	100	1	707	0	0	0	0	0	0	0	0	0	0	0	0;
	90	0	0	300	-300	0.985	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	91	0	0	100	-100	0.98	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	92	0	0	9	-3	0.99	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	99	0	0	100	-100	1.01	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	100	252	0	155	-50	1.017	100	1	352	0	0	0	0	0	0	0	0	0	0	0	0;
	103	40	0	40	-15	1.01	100	1	140	0	0	0	0	0	0	0	0	0	0	0	0;
	104	0	0	23	-8	0.971	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	105	0	0	23	-8	0.965	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	107	0	0	200	-200	0.952	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	110	0	0	23	-8	0.973	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	111	36	0	1000	-100	0.98	100	1	136	0	0	0	0	0	0	0	0	0	0	0	0;
	112	0	0	1000	-100	0.975	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	113	0	0	200	-100	0.993	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	116	0	0	1000	-1000	1.005	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% line data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
line118 = [
	10	2	0.0303	0.0999	0.0254  9900	0	0	0	0	1	-360	360;
	10	3	0.0129	0.0424	0.01082	9900	0	0	0	0	1	-360	360;
	4	5	0.00176	0.00798	0.0021	9900	0	0	0	0	1	-360	360;
	3	5	0.0241	0.108	0.0284	9900	0	0	0	0	1	-360	360;
	5	6	0.0119	0.054	0.01426	9900	0	0	0	0	1	-360	360;
	6	7	0.00459	0.0208	0.0055	9900	0	0	0	0	1	-360	360;
	8	9	0.00244	0.0305	1.162	9900	0	0	0	0	1	-360	360;
	8	5	0	0.0267	0	9900	0	0	0.985	0	1	-360	360;
	9	1	0.00258	0.0322	1.23	9900	0	0	0	0	1	-360	360;
	4	11	0.0209	0.0688	0.01748	9900	0	0	0	0	1	-360	360;
	5	11	0.0203	0.0682	0.01738	9900	0	0	0	0	1	-360	360;
	11	12	0.00595	0.0196	0.00502	9900	0	0	0	0	1	-360	360;
	2	12	0.0187	0.0616	0.01572	9900	0	0	0	0	1	-360	360;
	3	12	0.0484	0.16	0.0406	9900	0	0	0	0	1	-360	360;
	7	12	0.00862	0.034	0.00874	9900	0	0	0	0	1	-360	360;
	11	13	0.02225	0.0731	0.01876	9900	0	0	0	0	1	-360	360;
	12	14	0.0215	0.0707	0.01816	9900	0	0	0	0	1	-360	360;
	13	15	0.0744	0.2444	0.06268	9900	0	0	0	0	1	-360	360;
	14	15	0.0595	0.195	0.0502	9900	0	0	0	0	1	-360	360;
	12	16	0.0212	0.0834	0.0214	9900	0	0	0	0	1	-360	360;
	15	17	0.0132	0.0437	0.0444	9900	0	0	0	0	1	-360	360;
	16	17	0.0454	0.1801	0.0466	9900	0	0	0	0	1	-360	360;
	17	18	0.0123	0.0505	0.01298	9900	0	0	0	0	1	-360	360;
	18	19	0.01119	0.0493	0.01142	9900	0	0	0	0	1	-360	360;
	19	20	0.0252	0.117	0.0298	9900	0	0	0	0	1	-360	360;
	15	19	0.012	0.0394	0.0101	9900	0	0	0	0	1	-360	360;
	20	21	0.0183	0.0849	0.0216	9900	0	0	0	0	1	-360	360;
	21	22	0.0209	0.097	0.0246	9900	0	0	0	0	1	-360	360;
	22	23	0.0342	0.159	0.0404	9900	0	0	0	0	1	-360	360;
	23	24	0.0135	0.0492	0.0498	9900	0	0	0	0	1	-360	360;
	23	25	0.0156	0.08	0.0864	9900	0	0	0	0	1	-360	360;
	26	25	0	0.0382	0	9900	0	0	0.96	0	1	-360	360;
	25	27	0.0318	0.163	0.1764	9900	0	0	0	0	1	-360	360;
	27	28	0.01913	0.0855	0.0216	9900	0	0	0	0	1	-360	360;
	28	29	0.0237	0.0943	0.0238	9900	0	0	0	0	1	-360	360;
	30	17	0	0.0388	0	9900	0	0	0.96	0	1	-360	360;
	8	30	0.00431	0.0504	0.514	9900	0	0	0	0	1	-360	360;
	26	30	0.00799	0.086	0.908	9900	0	0	0	0	1	-360	360;
	17	31	0.0474	0.1563	0.0399	9900	0	0	0	0	1	-360	360;
	29	31	0.0108	0.0331	0.0083	9900	0	0	0	0	1	-360	360;
	23	32	0.0317	0.1153	0.1173	9900	0	0	0	0	1	-360	360;
	31	32	0.0298	0.0985	0.0251	9900	0	0	0	0	1	-360	360;
	27	32	0.0229	0.0755	0.01926	9900	0	0	0	0	1	-360	360;
	15	33	0.038	0.1244	0.03194	9900	0	0	0	0	1	-360	360;
	19	34	0.0752	0.247	0.0632	9900	0	0	0	0	1	-360	360;
	35	36	0.00224	0.0102	0.00268	9900	0	0	0	0	1	-360	360;
	35	37	0.011	0.0497	0.01318	9900	0	0	0	0	1	-360	360;
	33	37	0.0415	0.142	0.0366	9900	0	0	0	0	1	-360	360;
	34	36	0.00871	0.0268	0.00568	9900	0	0	0	0	1	-360	360;
	34	37	0.00256	0.0094	0.00984	9900	0	0	0	0	1	-360	360;
	38	37	0	0.0375	0	9900	0	0	0.935	0	1	-360	360;
	37	39	0.0321	0.106	0.027	9900	0	0	0	0	1	-360	360;
	37	40	0.0593	0.168	0.042	9900	0	0	0	0	1	-360	360;
	30	38	0.00464	0.054	0.422	9900	0	0	0	0	1	-360	360;
	39	40	0.0184	0.0605	0.01552	9900	0	0	0	0	1	-360	360;
	40	41	0.0145	0.0487	0.01222	9900	0	0	0	0	1	-360	360;
	40	42	0.0555	0.183	0.0466	9900	0	0	0	0	1	-360	360;
	41	42	0.041	0.135	0.0344	9900	0	0	0	0	1	-360	360;
	43	44	0.0608	0.2454	0.06068	9900	0	0	0	0	1	-360	360;
	34	43	0.0413	0.1681	0.04226	9900	0	0	0	0	1	-360	360;
	44	45	0.0224	0.0901	0.0224	9900	0	0	0	0	1	-360	360;
	45	46	0.04	0.1356	0.0332	9900	0	0	0	0	1	-360	360;
	46	47	0.038	0.127	0.0316	9900	0	0	0	0	1	-360	360;
	46	48	0.0601	0.189	0.0472	9900	0	0	0	0	1	-360	360;
	47	49	0.0191	0.0625	0.01604	9900	0	0	0	0	1	-360	360;
	42	49	0.0715	0.323	0.086	9900	0	0	0	0	1	-360	360;
	42	49	0.0715	0.323	0.086	9900	0	0	0	0	1	-360	360;
	45	49	0.0684	0.186	0.0444	9900	0	0	0	0	1	-360	360;
	48	49	0.0179	0.0505	0.01258	9900	0	0	0	0	1	-360	360;
	49	50	0.0267	0.0752	0.01874	9900	0	0	0	0	1	-360	360;
	49	51	0.0486	0.137	0.0342	9900	0	0	0	0	1	-360	360;
	51	52	0.0203	0.0588	0.01396	9900	0	0	0	0	1	-360	360;
	52	53	0.0405	0.1635	0.04058	9900	0	0	0	0	1	-360	360;
	53	54	0.0263	0.122	0.031	9900	0	0	0	0	1	-360	360;
	49	54	0.073	0.289	0.0738	9900	0	0	0	0	1	-360	360;
	49	54	0.0869	0.291	0.073	9900	0	0	0	0	1	-360	360;
	54	55	0.0169	0.0707	0.0202	9900	0	0	0	0	1	-360	360;
	54	56	0.00275	0.00955	0.00732	9900	0	0	0	0	1	-360	360;
	55	56	0.00488	0.0151	0.00374	9900	0	0	0	0	1	-360	360;
	56	57	0.0343	0.0966	0.0242	9900	0	0	0	0	1	-360	360;
	50	57	0.0474	0.134	0.0332	9900	0	0	0	0	1	-360	360;
	56	58	0.0343	0.0966	0.0242	9900	0	0	0	0	1	-360	360;
	51	58	0.0255	0.0719	0.01788	9900	0	0	0	0	1	-360	360;
	54	59	0.0503	0.2293	0.0598	9900	0	0	0	0	1	-360	360;
	56	59	0.0825	0.251	0.0569	9900	0	0	0	0	1	-360	360;
	56	59	0.0803	0.239	0.0536	9900	0	0	0	0	1	-360	360;
	55	59	0.04739	0.2158	0.05646	9900	0	0	0	0	1	-360	360;
	59	60	0.0317	0.145	0.0376	9900	0	0	0	0	1	-360	360;
	59	61	0.0328	0.15	0.0388	9900	0	0	0	0	1	-360	360;
	60	61	0.00264	0.0135	0.01456	9900	0	0	0	0	1	-360	360;
	60	62	0.0123	0.0561	0.01468	9900	0	0	0	0	1	-360	360;
	61	62	0.00824	0.0376	0.0098	9900	0	0	0	0	1	-360	360;
	63	59	0	0.0386	0	9900	0	0	0.96	0	1	-360	360;
	63	64	0.00172	0.02	0.216	9900	0	0	0	0	1	-360	360;
	64	61	0	0.0268	0	9900	0	0	0.985	0	1	-360	360;
	38	65	0.00901	0.0986	1.046	9900	0	0	0	0	1	-360	360;
	64	65	0.00269	0.0302	0.38	9900	0	0	0	0	1	-360	360;
	49	66	0.018	0.0919	0.0248	9900	0	0	0	0	1	-360	360;
	49	66	0.018	0.0919	0.0248	9900	0	0	0	0	1	-360	360;
	62	66	0.0482	0.218	0.0578	9900	0	0	0	0	1	-360	360;
	62	67	0.0258	0.117	0.031	9900	0	0	0	0	1	-360	360;
	65	66	0	0.037	0	9900	0	0	0.935	0	1	-360	360;
	66	67	0.0224	0.1015	0.02682	9900	0	0	0	0	1	-360	360;
	65	68	0.00138	0.016	0.638	9900	0	0	0	0	1	-360	360;
	47	69	0.0844	0.2778	0.07092	9900	0	0	0	0	1	-360	360;
	49	69	0.0985	0.324	0.0828	9900	0	0	0	0	1	-360	360;
	68	69	0	0.037	0	9900	0	0	0.935	0	1	-360	360;
	69	70	0.03	0.127	0.122	9900	0	0	0	0	1	-360	360;
	24	70	0.00221	0.4115	0.10198	9900	0	0	0	0	1	-360	360;
	70	71	0.00882	0.0355	0.00878	9900	0	0	0	0	1	-360	360;
	24	72	0.0488	0.196	0.0488	9900	0	0	0	0	1	-360	360;
	71	72	0.0446	0.18	0.04444	9900	0	0	0	0	1	-360	360;
	71	73	0.00866	0.0454	0.01178	9900	0	0	0	0	1	-360	360;
	70	74	0.0401	0.1323	0.03368	9900	0	0	0	0	1	-360	360;
	70	75	0.0428	0.141	0.036	9900	0	0	0	0	1	-360	360;
	69	75	0.0405	0.122	0.124	9900	0	0	0	0	1	-360	360;
	74	75	0.0123	0.0406	0.01034	9900	0	0	0	0	1	-360	360;
	76	77	0.0444	0.148	0.0368	9900	0	0	0	0	1	-360	360;
	69	77	0.0309	0.101	0.1038	9900	0	0	0	0	1	-360	360;
	75	77	0.0601	0.1999	0.04978	9900	0	0	0	0	1	-360	360;
	77	78	0.00376	0.0124	0.01264	9900	0	0	0	0	1	-360	360;
	78	79	0.00546	0.0244	0.00648	9900	0	0	0	0	1	-360	360;
	77	80	0.017	0.0485	0.0472	9900	0	0	0	0	1	-360	360;
	77	80	0.0294	0.105	0.0228	9900	0	0	0	0	1	-360	360;
	79	80	0.0156	0.0704	0.0187	9900	0	0	0	0	1	-360	360;
	68	81	0.00175	0.0202	0.808	9900	0	0	0	0	1	-360	360;
	81	80	0	0.037	0	9900	0	0	0.935	0	1	-360	360;
	77	82	0.0298	0.0853	0.08174	9900	0	0	0	0	1	-360	360;
	82	83	0.0112	0.03665	0.03796	9900	0	0	0	0	1	-360	360;
	83	84	0.0625	0.132	0.0258	9900	0	0	0	0	1	-360	360;
	83	85	0.043	0.148	0.0348	9900	0	0	0	0	1	-360	360;
	84	85	0.0302	0.0641	0.01234	9900	0	0	0	0	1	-360	360;
	85	86	0.035	0.123	0.0276	9900	0	0	0	0	1	-360	360;
	86	87	0.02828	0.2074	0.0445	9900	0	0	0	0	1	-360	360;
	85	88	0.02	0.102	0.0276	9900	0	0	0	0	1	-360	360;
	85	89	0.0239	0.173	0.047	9900	0	0	0	0	1	-360	360;
	88	89	0.0139	0.0712	0.01934	9900	0	0	0	0	1	-360	360;
	89	90	0.0518	0.188	0.0528	9900	0	0	0	0	1	-360	360;
	89	90	0.0238	0.0997	0.106	9900	0	0	0	0	1	-360	360;
	90	91	0.0254	0.0836	0.0214	9900	0	0	0	0	1	-360	360;
	89	92	0.0099	0.0505	0.0548	9900	0	0	0	0	1	-360	360;
	89	92	0.0393	0.1581	0.0414	9900	0	0	0	0	1	-360	360;
	91	92	0.0387	0.1272	0.03268	9900	0	0	0	0	1	-360	360;
	92	93	0.0258	0.0848	0.0218	9900	0	0	0	0	1	-360	360;
	92	94	0.0481	0.158	0.0406	9900	0	0	0	0	1	-360	360;
	93	94	0.0223	0.0732	0.01876	9900	0	0	0	0	1	-360	360;
	94	95	0.0132	0.0434	0.0111	9900	0	0	0	0	1	-360	360;
	80	96	0.0356	0.182	0.0494	9900	0	0	0	0	1	-360	360;
	82	96	0.0162	0.053	0.0544	9900	0	0	0	0	1	-360	360;
	94	96	0.0269	0.0869	0.023	9900	0	0	0	0	1	-360	360;
	80	97	0.0183	0.0934	0.0254	9900	0	0	0	0	1	-360	360;
	80	98	0.0238	0.108	0.0286	9900	0	0	0	0	1	-360	360;
	80	99	0.0454	0.206	0.0546	9900	0	0	0	0	1	-360	360;
	92	100	0.0648	0.295	0.0472	9900	0	0	0	0	1	-360	360;
	94	100	0.0178	0.058	0.0604	9900	0	0	0	0	1	-360	360;
	95	96	0.0171	0.0547	0.01474	9900	0	0	0	0	1	-360	360;
	96	97	0.0173	0.0885	0.024	9900	0	0	0	0	1	-360	360;
	98	100	0.0397	0.179	0.0476	9900	0	0	0	0	1	-360	360;
	99	100	0.018	0.0813	0.0216	9900	0	0	0	0	1	-360	360;
	100	101	0.0277	0.1262	0.0328	9900	0	0	0	0	1	-360	360;
	92	102	0.0123	0.0559	0.01464	9900	0	0	0	0	1	-360	360;
	101	102	0.0246	0.112	0.0294	9900	0	0	0	0	1	-360	360;
	100	103	0.016	0.0525	0.0536	9900	0	0	0	0	1	-360	360;
	100	104	0.0451	0.204	0.0541	9900	0	0	0	0	1	-360	360;
	103	104	0.0466	0.1584	0.0407	9900	0	0	0	0	1	-360	360;
	103	105	0.0535	0.1625	0.0408	9900	0	0	0	0	1	-360	360;
	100	106	0.0605	0.229	0.062	9900	0	0	0	0	1	-360	360;
	104	105	0.00994	0.0378	0.00986	9900	0	0	0	0	1	-360	360;
	105	106	0.014	0.0547	0.01434	9900	0	0	0	0	1	-360	360;
	105	107	0.053	0.183	0.0472	9900	0	0	0	0	1	-360	360;
	105	108	0.0261	0.0703	0.01844	9900	0	0	0	0	1	-360	360;
	106	107	0.053	0.183	0.0472	9900	0	0	0	0	1	-360	360;
	108	109	0.0105	0.0288	0.0076	9900	0	0	0	0	1	-360	360;
	103	110	0.03906	0.1813	0.0461	9900	0	0	0	0	1	-360	360;
	109	110	0.0278	0.0762	0.0202	9900	0	0	0	0	1	-360	360;
	110	111	0.022	0.0755	0.02	9900	0	0	0	0	1	-360	360;
	110	112	0.0247	0.064	0.062	9900	0	0	0	0	1	-360	360;
	17	113	0.00913	0.0301	0.00768	9900	0	0	0	0	1	-360	360;
	32	113	0.0615	0.203	0.0518	9900	0	0	0	0	1	-360	360;
	32	114	0.0135	0.0612	0.01628	9900	0	0	0	0	1	-360	360;
	27	115	0.0164	0.0741	0.01972	9900	0	0	0	0	1	-360	360;
	114	115	0.0023	0.0104	0.00276	9900	0	0	0	0	1	-360	360;
	68	116	0.00034	0.00405	0.164	9900	0	0	0	0	1	-360	360;
	12	117	0.0329	0.14	0.0358	9900	0	0	0	0	1	-360	360;
	75	118	0.0145	0.0481	0.01198	9900	0	0	0	0	1	-360	360;
	76	118	0.0164	0.0544	0.01356	9900	0	0	0	0	1	-360	360;
];

bus_gen=gen118(:,1);
type = bus118(:,2);
Vsp= bus118(:,8);
theta= bus118(:,9); %in degrees
PG= zeros(size(bus118,1), 1);
PG(bus_gen)=gen118(:,2)/100;
QG= zeros(size(bus118,1), 1);
QG(bus_gen)=gen118(:,3)/100;
PL= bus118(:,3)/100;
QL= bus118(:,4)/100;
Qmin= -ones(size(bus118,1),1)*3;
Qmax= ones(size(bus118,1),1)*3;

line118 = [line118(:, 1:5) ones(size(line118,1),1) zeros(size(line118,1),1)];
bus118 = [bus118(:,1), Vsp, theta, PG, QG, PL, QL , zeros(size(bus118,1),2), type, Qmax, Qmin];

switch num
    case 3
        bus=bus3;
        line = line3;
    case 14
        bus = bus14;
        line = line14;
    case 30
        bus = bus30;
        line = line30;
    case 9
        bus = bus9;
        line = line9;
    case 118
        bus = bus118;
        line = line118;
    case 301
        bus = bus30g1;
        line = line30g1;
    case 302
        bus = bus30g2;
        line = line30g2;
    case 303
        bus = bus30g3;
        line = line30g3;
end        
end