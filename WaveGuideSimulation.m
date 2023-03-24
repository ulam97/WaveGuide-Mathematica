(* WaveGuideSimulation.m *)
(*This is basically a comment, any line of code within asterisk-bracket pair is simply suppressed within  running session*) 
ClearAll["Global`*"] 

(*Initialization all in SI*) 
E0 = 1; 
wavlen = 1550 10^-9(*wavelength*); 
Ncr = 1.464 (*n core*); 
Ncl = 1.459(*n cladding*); 
\[Beta]0 = 2 Pi Ncr/wavlen (*initial beta0*); 
\[Beta]1 = \[Beta]0 (*Initializing beta1*); 
lex = 37 (*scaled domain length horizontal*); 
ley = 37(*scaled domain length vertical*); 
lcx = 10(*scaled core horizontal length*); 
lcy = 7 (*scaled core vertical length*); 

domx = lex 10^-6(*domain horizontal*); 
domy = ley 10^-6(*domain vertical*); 
gnumx = lex*5(*Grid numbers horizontal*); 
gnumy = ley*5(*Grid numbers vertical*); 

(*Core Coordinate*) 
coreXS = (20 gnumx)/lex + 2 (*x start*); 
coreXE = coreXS + (lcx gnumx)/lex(*x end*); 
coreYS = (20 gnumy)/ley + 2(*y start*); 
coreYE = coreYS + (lcy gnumy)/ley(*y end*);
dx = domx/gnumx (*dx*); 
dy = domy/gnumy (*dy*); 
\[Gamma] = 10^-3 (*beta error tolerance*); 
Efield0 = Array[If[gnumx + 3 > #1 > 2 && gnumy + 3 > #2 > 2, E0, 0] &,  {gnumy + 2, gnumx + 2}] (*initial Field*); 
Efield1 = Efield0 (*Initialize Field1*); 

(*Definitions*) 
\[Alpha][n_] := (((2 \[Pi] n)/ 
 wavlen)^2 - \[Beta]1^2 )(*long coefficient*); 
PxE1[p_, q_] := (  
 Efield1[[p, q]] - 2 Efield1[[p - 1, q]] + Efield1[[p - 2, q]])/  dx^2 (*Discrete backward 2nd derivative with respect to horizontal \ at pq for field1*); 
PyE1[p_, q_] := (  
 Efield1[[p, q]] - 2 Efield1[[p, q - 1]] + Efield1[[p, q - 2]])/  dy^2 (*Discrete backward 2nd derivative with respect to vertical at \ pq for field1*); 
update[p_, q_, n_] := ((( 
 2 Efield0[[p - 1, q]] - Efield0[[p - 2, q]])/(dx^2)) + (( 
 2 Efield0[[p, q - 1]] - Efield0[[p, q - 2]])/( 
 dy^2)) )/(\[Alpha][n] + 1/dx^2 + 1/ 
 dy^2) (*update field value at pq*); 
\[Beta]u := \[Sqrt]((\!\( 
\*UnderoverscriptBox[\(\[Sum]\), \(p = 3\), \(gnumy + 2\)]\(
\*UnderoverscriptBox[\(\[Sum]\), \(q = 3\), \(gnumx + 2\)]If[  coreXS - 1 < p < coreXE + 1 &&  
 coreYS - 1 < q <  
 coreYE +  
 1, \((Efield1[[p, q]]\ \((PxE1[p, q] + PyE1[p, q])\) +  \*SuperscriptBox[\(( 
\*FractionBox[\(2\ \[Pi]\ Ncr\), \(wavlen\)])\), \(2\)])\)\ dx\ dy, \ \((Efield1[[p, q]]\ \((PxE1[p, q] + PyE1[p, q])\) +  
\*SuperscriptBox[\(( 
\*FractionBox[\(2\ \[Pi]\ Ncl\), \(wavlen\)])\), \(2\)])\)\ dx\ dy]\)\ \))/(\!\( 
\*UnderoverscriptBox[\(\[Sum]\), \(p = 3\), \(gnumy\)]\( \*UnderoverscriptBox[\(\[Sum]\), \(q = 3\), \(gnumx\)] \*SuperscriptBox[\(Efield1[[p, q]]\), \(2\)]\)\))); (*beta updater*); 

(*Initial Computation*) 
For[j = 3, j < gnumy + 3, j++, 
 For[i = 3, i < gnumx + 3, i++, 
 Efield1[[i, j]] =  
 If[coreXS - 1 < i < coreXE + 1 && coreYS - 1 < j < coreYE + 1,   update[i, j, Ncr], update[i, j, Ncl]]]; 
 ]; 
Efield0 = Efield1; 
\[Beta]1 = \[Beta]u; 
Print[\[Beta]1];
Print[Abs[\[Beta]1 - \[Beta]0]/\[Beta]0]; 

(*Print[ListDensityPlot[Efield1]];*) 

(*Looping*) 
While[Abs[\[Beta]1 - \[Beta]0]/\[Beta]0 > \[Gamma] , \[Beta]0 = \[Beta]1; 
For[j = 3, j < gnumy + 3, j++, 
 For[i = 3, i < gnumx + 3, i++, 
 Efield1[[i, j]] =  
 If[coreXS - 1 < i < coreXE + 1 && coreYS - 1 < j < coreYE + 1,   update[i, j, Ncr], update[i, j, Ncl]]; 
 ]; 
 ]; 
Efield0 = Efield1; 
\[Beta]1 = \[Beta]u; 
Print[\[Beta]1] (*beta*); 
Print[Abs[\[Beta]1 - \[Beta]0]/\[Beta]0] (*beta error*); (*Print[ListDensityPlot[Efield1]] (*Plot*);*) 
]
