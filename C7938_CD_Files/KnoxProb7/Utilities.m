(*: Mathematica Version: 7.0 *)

(*: Package Version: 1.0 *)

(*: Name: "KnoxProb7`Utilities`" *)

(*: Context: "KnoxProb7`Utilities`" *)

(*: Title: Commands for Introduction to Probability 
with Mathematica, v. 7 *)

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
Updated to Mathematica 7 compatible version beginning Jan. 5, 2009
eliminating the local version of Histogram and the loading
of the obsolete BarCharts package in view of the corresponding
changes in version 7. Since the GeneralizedBarChart function
had been eliminated, the ProbabilityHistogram command had to be
rewritten to use primitive graphics objects instead of that
function.
 *)

(*: Reference: Usage messages only *)

(*: Summary:
	This package contains routines to draw histograms 
and dot plots of one-dimensional lists of data, and discrete
histograms of probability mass functions. It also contains
commands dealing with sampling, and some special plot
commands as described in the History above.
*)

(*: Keywords: dot plot, 
    permutations, plotting *)

(*: Requirements: Requires Mathematica v. 7 or above. *)

(*: Warning: None *)

(*: Sources: Standard statistics texts and stat package manuals *)

BeginPackage["KnoxProb7`Utilities`",
{"Combinatorica`",
"Utilities`FilterOptions`"}]


(* Usage messages *)

ProbabilityHistogram::usage=
"ProbabilityHistogram[statelist,problist] takes a list 
of states, not necessarily in order, an associated list of 
probabilities for those states, assumed to have the same 
length, and displays a bar chart for the probability distribution.
Options are those of Graphics, together with BarColor to
set the color of the bars."

BarColor::usage=
"BarColor is an option for ProbabilityHistogram that sets the 
color of the histogram bars. Default value is Blue."

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


(* ProbabilityHistogram takes a list of states, assumed 
to be in order, an associated list of probabilities for
those states, assumed to have the same length, and
graphic options for bars as in Histogram, and displays
a bar chart for the probability distribution. The code is
straightforward, except for the part where the minimum gap
between states is found.  The bar widths are all made equal
to the smaller of that minimum gap and 1.*)

 Options[ProbabilityHistogram]=Join[{BarColor->Blue,BaseStyle->{FontFamily->"Times",FontSize->8},Axes->True,AspectRatio->.5},Options[Graphics]];

 ProbabilityHistogram[statelist_,problist_,opts___]:=Module[{numstates,mingap,width,bcolor,baropts,slist,plist,tlist,
graphdirect,g1,g2},
  baropts=FilterOptions[Graphics,opts];
  bcolor=BarColor/.{opts}/.Options[ProbabilityHistogram];
  tlist=Ticks/.{opts}/.Options[Graphics];
  numstates=Length[statelist];
  {slist,plist}=Transpose[Sort[Transpose[{statelist,problist}],#1[[1]]<#2[[1]]&]];
  mingap=Min[Table[slist[[i]]-slist[[i-1]],{i,2,numstates}]];
  width=Min[{mingap,1}];
  rectlist=Table[Polygon[{{slist[[i]]-.5*width,0},{slist[[i]]-.5*width,plist[[i]]},{slist[[i]]+.5*width,plist[[i]]},{slist[[i]]+.5*width,0}}],{i,1,numstates}];
 graphdirect=Join[{bcolor,EdgeForm[Thin]}, rectlist];
 g1=Graphics[graphdirect, PlotRange->All, Axes->True, AxesOrigin->{slist[[1]]-.5*width,0}, Ticks->{slist,Automatic}, AspectRatio->.5];
g2=Graphics[graphdirect,PlotRange->All, Axes->True, AxesOrigin->{slist[[1]]-.5*width,0},Ticks->tlist, AspectRatio->.5];
If[ToString[tlist]==ToString[Automatic],g1,g2]]




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

Options[PlotStepFunction] = Join[{DotSize -> .017,StepStyle->RGBColor[0,0,0]}, Options[Plot]];

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

SetAttributes[ProbabilityHistogram,ReadProtected];
SetAttributes[BarColor,ReadProtected];
SetAttributes[DotPlot,ReadProtected];
SetAttributes[KPermutations,ReadProtected];
SetAttributes[RandomKPermutation,ReadProtected];
SetAttributes[DrawIntegerSample,ReadProtected];
SetAttributes[SimulateMysteryX,ReadProtected];
SetAttributes[Replacement,ReadProtected];
SetAttributes[Ordered,ReadProtected];
SetAttributes[PlotStepFunction,ReadProtected];
SetAttributes[DotSize,ReadProtected];
SetAttributes[PlotContsProb,ReadProtected];
SetAttributes[DotSize,ReadProtected];
SetAttributes[DotColor,ReadProtected];
SetAttributes[NumCategories,ReadProtected];
SetAttributes[VariableName,ReadProtected];
SetAttributes[StepStyle,ReadProtected];
SetAttributes[ShadingStyle,ReadProtected];

Protect[ProbabilityHistogram,BarColor,DotPlot,KPermutations,
RandomKPermutation,DrawIntegerSample,SimulateMysteryX,
Replacement,Ordered,PlotStepFunction,DotSize,
PlotContsProb,DotSize,DotColor,VariableName,NumCategories,
StepStyle,ShadingStyle];

EndPackage[]
