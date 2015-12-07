(*: Mathematica Version: 6.0 *)

(*: Package Version: 1.0 *)

(*: Name: "KnoxProb6`Utilities`" *)

(*: Context: "KnoxProb6`Utilities`" *)

(*: Title: Commands for Introduction to Probability 
with Mathematica, v. 6 *)

(*: Author: Kevin J. Hastings (Knox College) *)

(*: Copyright:  *)

(*: History:
	Version 1.0 begun Feb. 20, 2000, a subset of the 
KnoxStat packages. Histogram, ProbabilityHistogram, and 
DotPlot from the KnoxStat`StatGraphics` package, and 
DrawIntegerSample, SimulateMysteryX, KPermutations,
RandomKPermutation from KnoxStat`Sampling`, and also
new commands DensityHistogram, PlotStepFunction, and
PlotContsProb added as manuscript was rewritten.  Updated to
version 2.0 June 9, 2000 to incorporate changes in the larger KnoxStat
packages, particularly options for the graphics commands. 
Updated to Mathematica 6 compatible version beginning Sept. 5, 2007.

 *)

(*: Reference: Usage messages only *)

(*: Summary:
	This package contains routines to draw histograms 
and dot plots of one-dimensional lists of data, and discrete
histograms of probability mass functions. It also contains
commands dealing with sampling, and some special plot
commands as described in the History above.
*)

(*: Keywords: histogram, dot plot, 
    permutations, plotting *)

(*: Requirements: Requires Mathematica v. 6 or above. *)

(*: Warning: None *)

(*: Sources: Standard statistics texts and stat package manuals *)

BeginPackage["KnoxProb6`Utilities`",
{"BarCharts`",
"Combinatorica`",
"Utilities`FilterOptions`"}]


(* Usage messages *)

Histogram::usage=
"Histogram[datalist,numrectangles] plots a histogram of a list of data, 
with a desired number of rectangles. It inherits some of the options of 
GeneralizedBarChart and has four of its own.  The option Type has any
of the values Relative (default), Absolute, or Scaled depending on whether 
you want bars to have heights which are relative frequencies, absolute
frequencies, or relative frequencies divided by interval length.  The option
Endpoints may be set to a list {a,b} of real numbers with a<b to force
the histogram to be plotted between these endpoints. Otherwise the command
uses the min and max of the datalist as endpoints. The
option NumDigits (initialized to 2) can be used to set the number of decimal
digits used in the tick marks on the x-axis. The option
Distribution->Continuous may be reset to Discrete in order to force
a histogram whose boxes are at the integer values between the lowest
and highest integer data value. The user cannot override the PlotRange option,
nor AxesOrigin, nor Ticks, nor BarOrientation, in the interest of
having a well-formed graph."	

Type::usage=
"Type is an option for Histogram that can be given values Relative (default),
Absolute, or Scaled for a histogram that uses relative frequencies,
absolute frequencies, or frequencies scaled so that the total bar
area is 1."

Absolute::usage=
"Absolute is a possible value for option Type in Histogram, which
means that the histogram bars are absolute frequencies."

Scaled::usage=
"Scaled is a possible value for option Type in Histogram, which 
means that the histogram bars are relative frequencies divided by
interval lengths."

Relative::usage=
"Relative is a possible value for option Type in Histogram, which
means that the histogram bars are relative frequencies."

Endpoints::usage=
"Endpoints is an option for Histogram that can be set to a list {a,b}
with a<b to force the histogram to begin at a and end at b."

NumDigits::usage=
"NumDigits is an option for Histogram and 
DensityHistogram that controls the number of decimal 
digits used in x-axis category tickmarks."

Distribution::usage=
"Distribution is an option for Histogram which is given
the default value Continuous which is appropriate for
histograms of general data sets, but if it is set to 
Discrete, a histogram is created in which all integers
from the lowest integer data value to the highest are
made into their own categories."

Continuous::usage=
"Continuous is a designation for the kind of distribution
sampled from in the Distribution option of the Histogram
command."

DensityHistogram::usage=
"DensityHistogram[densityfn,n,a,b] takes a named density 
function, a number n indicating the number of subintervals to use
in a discrete approximation of the density function, and
the endpoints a and b of the interval over which the
density is defined, and returns a combined plot of the
density and the discrete approximation.  It accepts an
option NumDigits->2 to control how many digits appear in the
x-axis tick marks, and it accepts all options of
GeneralizedBarChart."

ProbabilityHistogram::usage=
"ProbabilityHistogram[statelist,problist] takes a list 
of states, assumed to be in order, an associated list of 
probabilities for those states, assumed to have the same 
length, and graphic options for bars as in Histogram, 
and displays a bar chart for the probability distribution.
The user cannot override the PlotRange option,
nor AxesOrigin, nor Ticks, nor BarOrientation, in the interest of
having a well-formed graph."

DotPlot::usage=
"DotPlot[datalist] draws a dot plot of the given list of data. It takes
some of the options of ListPlot, although in the interest of getting
a good dotplot it suppresses AspectRatio, AxesOrigin, Axes,
PlotRange, and AxesLabel. It has four options of its own; VariableName,
set to be a string with which to label the horizontal axis,
NumCategories, initialized to 30, which gives the number
of categories for stacking dots, DotSize, initially .02,
and DotColor, initially RGBColor[0,0,0], that is black."

DotSize::usage=
"DotSize is an option for PlotStepFunction and DotPlot. It controls 
the size of the dots on the graph."

VariableName::usage=
"VariableName is an option used by DotPlot to give the horizontal axis 
only a variable name. It suppresses the use of option AxesLabel."

DotColor::usage=
"DotColor is an option for DotPlot which applies a style to the dots.  Its default is black RGBColor[0,0,0]."

NumCategories::usage=
"NumCategories->30 is an option for DotPlot which
controls the number of categories into which the axis is divided;
Dots falling in the same category are stacked."

KPermutations::usage=
"KPermutations[fromlist,k] returns a list of all permutations of k objects 
taken from the given list."

RandomKPermutation::usage=
"RandomKPermutation[fromlist,k] returns a randomly selected permutation of k 
objects from the given list."	

DrawIntegerSample::usage=
"DrawIntegerSample[a, b, n, options] selects a sample of size n from the 
set of integers a,..., b with optional Boolean arguments Ordered->True,
and Replacement->False determining the four possible sampling scenarios."

SimulateMysteryX::usage=
"SimulateMysteryX[m] returns a list of m values 
simulated from a discrete distribution. It is used
for a problem in Section 1.1 of Hastings."

Replacement::usage=
"Replacement is an option for DrawIntegerSample which
is boolean valued and specifies whether the sample is
to be chosen with replacement. Default is False"

Ordered::usage=
"Ordered is an option for DrawIntegerSample which
is boolean valued and specifies whether the sample is
to be chosen with regard to order. Default is True"

PlotStepFunction::usage=
"PlotStepFunction[fn,domain,jumplist] plots a step
function on domain specified, with jumps at the points in 
jumplist, which is a list of sorted numbers.  The step function 
is assumed to be right continuous, as a c.d.f. is. It accepts option 
DotSize->.017 to change the size of the dots, StepStyle->RGBColor[0,0,0]
to assign a style to the steps, and it inherits any options that
make sense for Show."

StepStyle::usage=
"StepStyle is an option for PlotStepFunction which gives a style
to the steps. Its default is RGBColor[0,0,0], or black."

PlotContsProb::usage=
"PlotContsProb[density,domain,between] plots the area
under the given function on the given domain between
the points in the list between, which is assumed to consist 
of two points in increasing order.  Options are the
options that make sense for Show, and ShadingStyle->RGBColor[1,0,0] 
which can be used to give a style to the shaded area region."

ShadingStyle::usage=
"ShadingStyle is an option for PlotContsProb that sets the
style of the shaded region. Default is RGBColor[1,0,0]."

Begin["Private`"]

(* Trunc is a private function used by Histogram
to control the number of decimal digits printed
on the horizontal axis categories *)

Trunc[number_, digits_] := 
   Round[10^digits*number]/((1.0)*10^digits)

(* Histogram plots a histogram of a list of data, 
  with a desired number of rectangles. The Type option can be set
 to Absolute or Scaled according to whether the user wants rectangle
heights to be absolute frequencies, or relative freqs divided by 
interval length so as to scale the area to 1.  By default Type is 
Relative, which makes rectangle heights the relative freqs. Can also set
option Endpoints to a list {a,b} where a is less than b to force
histogram endpoints, otherwise by default it uses the min and max
of the data list. Option Distribution is Continuous by default,
in which case the interval from low to high is split up into
numrectangles equally sized subintervals. If set to anything else,
Histogram assumes a discrete integer valued distribution taking
values between low and high, and uses as categories all integers
in that range. Uses BinCounts to count observations in categories. Note that because barstend to be close together on the x-axis by default we choose the
small font TimesNewRoman 8 pt. *)

Options[Histogram] =
  Join[{Type -> Relative, Distribution -> Continuous, Endpoints -> {}, 
        NumDigits -> 2, BaseStyle -> {FontFamily->"Times", FontSize -> 8},
         DisplayFunction->$DisplayFunction}, Options[GeneralizedBarChart]];

 Histogram[datalist_, numrectangles_, opts___] :=
	Module[{numdigs, htype, ends, dist, baropts, sampsize,
      low, high, dx, cutoffs, hts, midpoints, widths, states, 
      numstates,dispfn,bstyle},
		numdigs = NumDigits /. {opts} /. Options[Histogram];
		htype = Type /. {opts} /. Options[Histogram];
		ends = Endpoints /. {opts} /. Options[Histogram];
    	      dist = Distribution /. {opts} /. Options[Histogram];
		dispfn = DisplayFunction/. {opts} /. Options[Histogram];
            bstyle = BaseStyle/. {opts} /. Options[Histogram];
		baropts = FilterOptions[GeneralizedBarChart, opts];
		sampsize = Length[datalist];
            If[ToString[dist] == ToString[Continuous],
		If[Length[ends] == 2, low = ends[[1]]; high = ends[[2]],
			 low = Min[datalist]; high = Max[datalist]];
		dx = N[(high - low)/numrectangles];
		cutoffs = Append[Table[low + i*dx, {i, 0, numrectangles - 1}],Infinity]; 
		If[ToString[htype] == ToString[Scaled],
      		hts = N[BinCounts[datalist, {cutoffs}]/(sampsize*dx)],
		If[ToString[htype] == ToString[Absolute],  
			hts = N[BinCounts[datalist, {cutoffs}]],
			hts = N[BinCounts[datalist, {cutoffs}]/sampsize]]];
		midpoints = Table[low + (i + .5)*dx, {i, 0, numrectangles - 1}];
		widths = Table[dx, {i, 1, numrectangles}];
		GeneralizedBarChart[Transpose[{midpoints, hts, widths}],
		PlotRange -> All, 
		AxesOrigin -> {low, 0},BarOrientation->Vertical,
		Ticks -> {Trunc[midpoints, numdigs], Automatic}, 
                DisplayFunction->dispfn,BaseStyle->bstyle,baropts], 
      	low = Min[datalist]; high = Max[datalist];
      	states = Table[i, {i, low, high}];
		cutoffs = Append[states,Infinity];
      	numstates = Length[states];
      	widths = Table[1, {i, 1, numstates}];
      	If[ToString[htype] == ToString[Absolute],
               hts = N[BinCounts[datalist, {cutoffs}]],
      		hts = N[BinCounts[datalist, {cutoffs}]/(sampsize*1)]];
      	GeneralizedBarChart[Transpose[{states, hts, widths}],
		PlotRange -> All, BarOrientation->Vertical,
		AxesOrigin -> {low - .5, 0},
		Ticks -> {states, Automatic}, 
                DisplayFunction->dispfn,BaseStyle->bstyle,baropts]]]


(* ProbabilityHistogram takes a list of states, assumed 
to be in order, an associated list of probabilities for
those states, assumed to have the same length, and
graphic options for bars as in Histogram, and displays
a bar chart for the probability distribution. The code is
straightforward, except for the part where the minimum gap
between states is found.  The bar widths are all made equal
to that minimum gap.*)

Options[ProbabilityHistogram]=Join[{BaseStyle -> {FontFamily->"Times", FontSize -> 8}, DisplayFunction->$DisplayFunction}, Options[GeneralizedBarChart]];

ProbabilityHistogram[statelist_,problist_,opts___]:=
	Module[{numstates,mingap,widths},
	  baropts = FilterOptions[GeneralizedBarChart, opts];
	  numstates=Length[statelist];
	  mingap=Min[Table[statelist[[i]]-statelist[[i-1]],{i,2,numstates}]];
	  widths=Table[mingap,{i,1,numstates}];
	  GeneralizedBarChart[Transpose[{statelist,problist,widths}],
		     PlotRange->All,
		AxesOrigin->{statelist[[1]]-.5*mingap,0},
		Ticks->{statelist,Automatic}, BarOrientation->Vertical,
                baropts]]

(* DensityHistogram takes a named density function,
a number n indicating the number of subintervals to use
in a discrete approximation of the density function, and
the endpoints a and b of the interval over which the
density is defined, and returns a combined plot of the
density and the discrete approximation.  It accepts an
option NumDigits->2 to control how many digits appear in the
x-axis tick marks, and it accepts all options of
GeneralizedBarChart *)


Options[DensityHistogram]=
  Join[{NumDigits->2},Options[GeneralizedBarChart]]

DensityHistogram[fn_, {x_,a_,b_}, n_, opts___] :=
	Module[{numdigs,hstyle,hvalues,hedges,hedgestyle,
		horientation,Deltax, ptlist, densitylist, 
		widthlist, g1, g2,F},
		hstyle=BarStyle/.{opts}/.Options[DensityHistogram];
		hvalues=BarValues/.{opts}/.Options[DensityHistogram];
		hedges=BarEdges/.{opts}/.Options[DensityHistogram];
		hedgestyle=BarEdgeStyle/.{opts}/.Options[DensityHistogram];
		horientation=BarOrientation/.{opts}/.Options[DensityHistogram];
		F:= Function[x,fn];
		Deltax = (b - a)/n;
		ptlist = N[Table[a + k*Deltax, {k, 0, n}]];
		densitylist = Map[F[#]&,ptlist];
		widthlist = Table[Deltax, {k, 0, n}];
		g1 = GeneralizedBarChart[Transpose[{ptlist, densitylist, widthlist}],  
			BaseStyle -> {FontFamily->"Times", FontSize -> 8},
			BarStyle->hstyle,
		     	BarValues->hvalues,
			BarEdges->hedges,
			BarEdgeStyle->hedgestyle,
			BarOrientation->horientation];
		g2 = Plot[F[x], {x, a, b}, 
        		BaseStyle -> {FontFamily->"Times", FontSize -> 8}];
		Show[g1, g2,
			AxesOrigin -> {a - Deltax, 0}, 
      			Ticks -> {Trunc[ptlist, 2], Automatic}]]

(* DotPlot displays a dot plot of a list of data. It takes
the options of ListPlot, although in the interest of getting
a good dotplot it suppresses most of them, such as keeping a y-axis
in the picture. It has three options of its own; VariableName,
set to be a string with which to label the horizontal axis,
NumCategories, initialized to 30, which gives the number
of categories for stacking dots, Dot Color, set to black initially,
and DotSize, initially .02. After initializing the options, it finds 
the largest and smallest data points and splits up the interval 
between into equally sized subintervals.  The number of observations 
in each is counted, and a data list of dot positions is built in 
which the horizontal dot positions are at the left endpoints of each 
category, and the vertical coordinates are the positive integers up 
to the frequency for the category.  ListPlot is then called, with some
options that by trial and error seem to produce a nice graph.*)

Options[DotPlot]=Join[{VariableName->" ", NumCategories->30,
			DotSize->.02,DotColor->RGBColor[0,0,0]}, 
                     Options[ListPlot]];

DotPlot[datalist_,opts___]:=
	Module[{low,high,dx,cutoffs,freqs,midpoints,dotlist,
		maxcount,numcats,listplotopts,dsize,varname,dcolor},
		dsize=DotSize/.{opts}/.Options[DotPlot];
		dcolor=DotColor/.{opts}/.Options[DotPlot];
		numcats=NumCategories/.{opts}/.Options[DotPlot];
		varname=VariableName/.{opts}/.Options[DotPlot];
		listplotopts=FilterOptions[ListPlot,opts];
		low=Min[datalist];
		high=Max[datalist];
		dx=(high-low)/numcats;
		cutoffs = Append[Table[low+i*dx,{i,0,numcats-1}],Infinity];
		freqs=BinCounts[datalist,{cutoffs}];
		midpoints=Table[low+i*dx,{i,0,numcats-1}];
		dotlist={};
		Do[Do[AppendTo[dotlist,{midpoints[[i]],j}],{j,1,freqs[[i]]}],{i,1,
        numcats}];
		maxcount=Max[freqs];
		ListPlot[dotlist,PlotStyle->{PointSize[dsize],dcolor},
            AspectRatio->maxcount/(.5*numcats),
      AxesOrigin->{low-dx,0},Axes->{True,False},
            PlotRange->{{low-dx,high+dx},{0,1.1*maxcount}},
		AxesLabel->{varname," "},
		listplotopts]]


(* KPermutations returns all permutations of k objects from the given list, 
by first using KSubsets to create the subsets, and then joining up the 
permutations into the complete list. *)


KPermutations[fromlist_,k_]:=
	Module[{currlist,subsetlist,len},
			subsetlist=KSubsets[fromlist,k];
			len = Length[subsetlist];
		   currlist={};
			Do[currlist=Join[currlist,Permutations[subsetlist[[i]]]],{i,1,len}];
				currlist]


(* RandomKPermutation returns a list of k objects sampled randomly in 
sequence without replacement. First the RandomKSubset function is used 
to get a random set, then the selected set is permuted randomly with 
RandomPermutation.*)

RandomKPermutation[fromlist_,k_]:=
	Module[{sample,intperm},
		   sample = RandomKSubset[fromlist,k];
		intperm=RandomPermutation[k];
		Table[sample[[intperm[[i]]]],{i,1,k}]]


Options[DrawIntegerSample]={Ordered->True,Replacement->False}

(* DrawIntegerSample selects a sample of size n from the set of integers 
a,..., b with optional Boolean arguments Ordered, and Replacement 
determining the four possible sampling scenarios. Default values for the 
optional arguments are in order and without replacement. It needs the 
discrete uniform distribution, and also 
uses the RandomKPermutation function in this package and RandomKSubset in 
the Combinatorica package. The four cases are considered in turn
by nested If statements; first in the with replacement case, if 
the sample is to be drawn in sequence it suffices to call on RandomArray,
and if the order does not matter we can still do that but we
eliminate the apparent sequencing in the output by sorting it in
ascending order. Else in the no replacement case, the order and
no order subcases are covered respectively by RandomKPermutation
and RandomKSubset.*)

DrawIntegerSample[a_,b_,n_,opts___]:=
	Module[{inorder,withrepl,thesample,population,i},
		inorder=Ordered/.{opts}/.Options[DrawIntegerSample];
		withrepl=
			Replacement/.{opts}/.Options[DrawIntegerSample];
		If[withrepl,
				If[inorder, thesample=RandomInteger[DiscreteUniformDistribution[{a, b}],n],
				    thesample=Sort[RandomInteger[DiscreteUniformDistribution[{a,b}],n]]],
			population=Table[i,{i,a,b}];
			 If[inorder,thesample=RandomKPermutation[population,n],
				thesample=RandomKSubset[population,n]]]]



(* SimulateMysteryX is a command for a Section 1.1 exercise
 which simulates a list of m values of the random variable
 X which equals 1 with probability 1/6, 2 with probability
 1/3, and 3 with probability 1/2 *)

SimulateMysteryX[m_]:=
	Module[{thelist,newrand,i},
		thelist={};
		Do[newrand=Random[];
			If[newrand<1/6,AppendTo[thelist,1],
				  If[newrand<1/2, AppendTo[thelist,2],
					          AppendTo[thelist,3]]],{i,1,m}];
		 thelist]

(* PlotStepFunction[fn,domain,jumplist] plots a step
function on the specified domain, with jumps at the points 
in jumplist, which is a list of sorted numbers.  The step function 
is assumed to be right continuous, as a c.d.f. is. It accepts option 
DotSize->.017 to change the size of the dots,  and it inherits the 
options of Plot. After reconstructing the input function as F[x],
a list of points (x,F(x)) at the jumps is created, then a list
of horizontal line segments between jump points.  By default the
x-axis is labeled with tick marks at the jumps. Later I added
StepStyle in order to set a style for the steps. *)

Options[PlotStepFunction] = Join[{DotSize -> .017,StepStyle->RGBColor[0,0,0]},
                                Options[Plot]];

PlotStepFunction[fn_, {x_,a_,b_}, jumplist_, opts___] :=
  Module[{plotopts, dsize, len, xticklist,F,stepsty,x},
    plotopts=FilterOptions[Plot,opts];
    dsize = DotSize /. {opts} /. Options[PlotStepFunction];
    stepsty=StepStyle/.{opts}/.Options[PlotStepFunction];
    len = Length[jumplist];
    F:=Function[x,fn];
    g1 = Graphics[
        Join[{stepsty},{PointSize[dsize]}, {Point[{a, F[a]}]},
           Table[Point[{jumplist[[i]], F[jumplist[[i]]]}], {i, 1, len}]]];
    g2 = Graphics[
        Join[{stepsty},{Line[{{a, F[a]}, {jumplist[[1]], F[a]}}]}, 
          Table[Line[{{jumplist[[i]], F[jumplist[[i]]]}, {jumplist[[i + 1]], 
                  F[jumplist[[i]]]}}], {i, 1, 
              len - 1}], {Line[{{jumplist[[len]], 
                  F[jumplist[[len]]]}, {b, F[jumplist[[len]]]}}]}]];
    xticklist = Join[{a}, jumplist, {b}];
    Show[{g1, g2}, plotopts, Ticks -> {xticklist, Automatic},
          Axes -> True]]

(* PlotContsProb[density,domain,between] plots the area
under the given function on the given domain, between
the points in the list between, which is assumed to consist 
of two points in increasing order.  Options are the
options for Plot, and one of its own, ShadingStyle, which can
assign a style to the shaded region. The plot interval bounds 
a,b and shaded interval bounds c and d are extracted from the input, 
and the shaded interval is split into 30 equal subintervals. 
Then the density function is reconstituted with name f and
variable x, and the polygon beginning with the leftmost
shaded point on the x-axis, and proceding through the polygonal
curve approximation, then back down to the x-axis is made. 
The plots of the shaded polygon and the density function are
then superimposed. *)

Options[PlotContsProb] = Join[{ShadingStyle->RGBColor[1, 0, 0],PlotRange->All},
Options[Plot]];

PlotContsProb[density_, {x_,a_,b_}, between_, opts___] :=
  Module[{c, d, f, xlist, n, dx, polylist, g1, g2,dispfn,
shade,pltrng,pstyle,tks,aor,plotopts},
    dispfn=DisplayFunction/.{opts}/.Options[PlotContsProb];
    shade=ShadingStyle/.{opts}/.Options[PlotContsProb];
    pltrng=PlotRange/.{opts}/.Options[PlotContsProb];
    pstyle=PlotStyle/.{opts}/.Options[PlotContsProb];
    tks=Ticks/.{opts}/.Options[PlotContsProb];
    aor=AxesOrigin/.{opts}/.Options[PlotContsProb];
    plotopts=FilterOptions[Plot,opts];
    n = 30;
    c = between[[1]];
    d = between[[2]];
    dx = (d - c)/n;
    f:= Function[x,density];
    xlist = Table[c + k*dx, {k, 0, n}];
    polylist = 
      Join[{{xlist[[1]], 0}}, 
        Table[{xlist[[i]], f[xlist[[i]]]}, {i, 1, 
            n + 1}], {{xlist[[n + 1]], 0}}];
    g1 = Plot[f[x], {x, a, b}, PlotRange->pltrng, PlotStyle->pstyle,Ticks->tks,AxesOrigin->aor (*plotopts*)];
    g2 = Graphics[{shade, Polygon[polylist]}];
    Show[{g1, g2}]]

End[]

SetAttributes[Trunc,{Listable,ReadProtected}];
SetAttributes[Histogram,ReadProtected];
SetAttributes[ProbabilityHistogram,ReadProtected];
SetAttributes[DotPlot,ReadProtected];
SetAttributes[Type,ReadProtected];
SetAttributes[Absolute,ReadProtected];
SetAttributes[Scaled,ReadProtected];
SetAttributes[Relative,ReadProtected];
SetAttributes[Endpoints,ReadProtected];
SetAttributes[DensityHistogram,ReadProtected];
SetAttributes[KPermutations,ReadProtected];
SetAttributes[RandomKPermutation,ReadProtected];
SetAttributes[DrawIntegerSample,ReadProtected];
SetAttributes[SimulateMysteryX,ReadProtected];
SetAttributes[Replacement,ReadProtected];
SetAttributes[Ordered,ReadProtected];
SetAttributes[PlotStepFunction,ReadProtected];
SetAttributes[DotSize,ReadProtected];
SetAttributes[PlotContsProb,ReadProtected];
SetAttributes[Distribution,ReadProtected];
SetAttributes[Continuous,ReadProtected];
SetAttributes[DotSize,ReadProtected];
SetAttributes[DotColor,ReadProtected];
SetAttributes[NumCategories,ReadProtected];
SetAttributes[VariableName,ReadProtected];
SetAttributes[StepStyle,ReadProtected];
SetAttributes[ShadingStyle,ReadProtected];

Protect[Trunc,Histogram,ProbabilityHistogram,DotPlot,
Type,Absolute,Scaled,Relative,Endpoints,KPermutations,
RandomKPermutation,DrawIntegerSample,SimulateMysteryX,
Replacement,Ordered,PlotStepFunction,DotSize,
PlotContsProb,DensityHistogram,Distribution,
Continuous,DotSize,DotColor,VariableName,NumCategories,
StepStyle,ShadingStyle];

EndPackage[]
