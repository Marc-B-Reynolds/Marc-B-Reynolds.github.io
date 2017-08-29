---
layout: post
title:  Pigeonhole principle bias
tagline: and mapping integer sets
categories: [math]
tags : [bias, distribution, random]
description : a note on bias introduced by the pigeonhold principle
plotly: true
revision: 20170809
---

This is a quick note about how the Pigeonhole Principle applies to mapping integers to a smaller set. One usage is generating unbiased uniform random integers.

Imagine we have a machine which has a large mixing bin.  Into this machine we put $b$ physically identical balls except one is blue and the remainder are red.  After sufficient mixing the balls drop into a long narrow transparent tube in single file.  Connected to the tube is a chute that can distribute the balls to $n$ bins and the machine is set-up to use some method to distribute the balls equally into the bins (of course numbered from 0 to $n-1$). After the distribution is complete we call the bin with the blue ball the result.

Let's say we have 16 balls and 4 bins, then after distribution each has 4 balls.  All is good.  Each bin has a 4 in 16 chance of getting the ball each run which is a probability of $\frac{1}{n}$.

Now we add an extra bin and things go abit wrong.  Since 16 isn't a multiple of 5 we end up with 4 bins of three balls and one with four:

$$ \left\{\frac{4}{16}, \frac{3}{16}, \frac{3}{16}, \frac{3}{16}, \frac{3}{16} \right\} $$

If we were the repeated apply this process then we'd get a zero result about 25% and 1-4 each about 18.75% of the time instead of the desired 20%.  If we were to up the number of bins to 11, then we'd get 5 bins recieving 2 and 6 receiving 1.  So instead of the desire ~9.09% result for each we'd get 12.5% and 6.25% for 2 and 1 ball bins respectively.  One possible measure of "wrongness" is to sum up the too much/little amounts:

$$ 5\left(\frac{2}{16}-\frac{1}{11}\right)+6\left(\frac{1}{11}-\frac{1}{16}\right) = \frac{15}{44} \approx .341 $$

\\
Which says we're getting an incorrect answer ~34.1% of the time.  We can visualize this as the sum of areas: chop off the too much (the orange) and fill in the too little (the green) and that area is the wrongness or bias.

<div id="ex1" style="width:100%"></div>

If you're thinking to yourself: "Why don't we just distribute an equal number of balls to each bin?  And if the blue ball is still in the tube...re-run the process until it does make into one?"  Yes, yes, you're very clever.  That's the rejection method and why I'm bothering to write this to up.

The machine parts break down as follows:

* The mixer is a generating function.
* The position of the ball in the transparent tube is the integer result of that function.
* The chute to bin is the map from the larger to smaller set. My examples are equivalent to modulus.

<br>

------

Math brushstrokes
------

\\
Given $b$ balls and $n$ bins with $b \geq n$, then each bin receives at least $a$ balls:

$$ \begin{equation} \label{a}
a = \left\lfloor \frac{b}{n} \right\rfloor 
\end{equation} $$ 

\\
We have $e$ excess balls:

$$ \begin{equation} \label{excess}
e = b \bmod n = b-n \left\lfloor \frac{b}{n} \right\rfloor 
\end{equation}
$$

\\
which is zero if $b$ is equally divisible by $n$.  So we have $\left(n-e\right)$ bins with $a$ balls and $e$ bins with $\left(a+1\right)$.  We can express the bias (like the 11 bin example above) as:

$$
\left(n-e\right)\left(\frac{1}{n}-\frac{a}{b}\right) + e\left(\frac{a+1}{b}-\frac{1}{n}\right) 
$$

\\
Plugging in $a$ and $e$ and then reducing gives:

$$ \begin{equation} \label{dbias}
\frac{2 \left(b-n \left\lfloor \frac{b}{n}\right\rfloor \right) \left(n \left\lfloor \frac{b}{n}\right\rfloor
   -b+n\right)}{b~n} 
\end{equation} $$


\\
We can convert $\eqref{dbias}$ into a univariant continuous function replacing $b$ by 1 and $n$ with $x \in \left\[0,1\right\] $:

$$ \begin{equation} \label{cbias}
\frac{2 \left(1-x \left\lfloor \frac{1}{x}\right\rfloor \right) \left(x \left\lfloor
   \frac{1}{x}\right\rfloor +x-1\right)}{x}
\end{equation} $$

\\
Equation $\eqref{cbias}$ is effectively superimposing all choices of $b$ and $n$ into a single graph (with normalized range) where $x=\frac{n}{b}$.  The following plots $\eqref{cbias}$ and the sample points that occur for $b=32$ and $b=64$.


<div id="phole" style="width:100%"></div>

\\
The function is zero when $x=\frac{1}{i}$ where $i$ is a positive integer since it is a factor of some $b$. We can rewrite $\eqref{cbias}$ as an infinite piecewise function, one per range of an integer $i$ and it's successor.  Each of these produces the $i^{th}$ lobe from the right ($i$ or $i+1$ balls/bin).  Example $i=2$ is second lobe from the right: $[\frac{1}{i+1}, \frac{1}{i}] = [\frac{1}{3}, \frac{1}{2}]$


$$  \begin{equation} \label{ibias}
f_i\left(x\right) = 2\left(2i+1-\frac{1}{x}-i\left(i+1\right)x  \right) 
\end{equation} $$

\\
The derivative of $\eqref{ibias}$:

$$ \begin{equation} \label{dibias}
{f_i}^\prime\left(x\right) = 2\left(\frac{1}{x^2}-i\left(i+1\right) \right) 
\end{equation} $$

\\
Solving ${f_i}^\prime\left(x\right) = k$ for $x$:

$$ \begin{equation} \label{kdibias}
  \frac{1}{\sqrt{i+i^2+\frac{k}{2}}}
\end{equation} $$

\\
Eyeballing the continous plot it appears that the bias is bound by the line through the origin $\frac{1}{2}x$.

<div id="bound" style="width:100%"></div>

\\
Since the slope of the proposed bounding line is $\frac{1}{2}$, we can plug $k=\frac{1}{2}$ into $\eqref{kdibias}$ to find the $x$ value within a lobe that has that slope:

$$ \begin{equation} \label{slope}
  \frac{2}{1+2i}
\end{equation} $$

\\
Plugging $\eqref{slope}$ into $\eqref{ibias}$ gives the bias at these points: 

$$ \frac{1}{1+2i} $$

\\
Which is $\frac{x}{2}$, so each lobe kisses the line once and the line is an upper bound of the bias.  Okay I'm cheating by not showing the bias never goes over the line.  And we could formulate a much tigher bound. My goal in this part was to come up with rule of thumb that can be quickly performed in your head.

You have fixed $b$ and variable $n \leq n_{\text{max}}$ 

$$ \text{bias} \leq \frac{n_{\text{max}}}{2b} $$

\\
Another converative measure is a "stair-step" function that bounds the bias. Determine the furthest to right lobe of the range: $i = 1+w-\{LeadZeroCount}\left(n_{max}\right)$, where $w$ is the number of bits in working word.

The slope is zero when $\eqref{dibias}$ is zero:

$$ x = \left(i+i^2\right)^{-\frac{1}{2}} $$

\\
Plugging this back into the piecewise bias function yields:

$$ 2\left(1+2\left(i-\sqrt{i\left(i+1\right)}\right) \right) $$

\\
then the first 128 lobe peaks are:


<div id="plot" style="width:100%"></div>

<br>

<script>

function pholeF(x)
{
  var a = 1./x;
  var b = 1.-x*Math.floor(a) 
  return 2.*b*(x-b)*a;
}

var phole;

{
  var cX=[];
  var cY=[];
  var limit = 0.5;
  var dx    = 1.0/512.0;
  var x = 1.0;
  for (x=0.0; x<=1.0; x+=1.0/2048.0) { cX.push(x); cY.push(pholeF(x)); }
  phole = {x:cX, y:cY, type:'scatter', mode:'line', name:'bias'};
}

function pholeB(x) { return .5*x; }
var phole2;

{
  var cX =[];
  var cY =[];
  var dx = .708/64.0;
  var x;

  for (x=0.0; x<=.708; x+=dx) { cX.push(x); cY.push(pholeB(x)); }
  
  phole2 = {x:cX, y:cY, type:'scatter', mode:'line', name:'bound'};
}


var phole;

function dphole(b) {
  var cX=[];
  var cY=[];
  var dx=1/b;
  for (var x=0.0; x<=1.0; x+=dx) { cX.push(x); cY.push(pholeF(x)); }
  return {x:cX, y:cY, type:'scatter', mode:'markers', name:'b='+b, marker:{size:4}};
}

var boundData = [phole, phole2];
var pholeData = [phole,dphole(64), dphole(32)];

var pholeLayout = {title: 'pigeonhole bias'};
Plotly.newPlot('phole', pholeData, pholeLayout, {showLink:false });

var pholeLayout = {title: 'bound'};
Plotly.newPlot('bound', boundData, pholeLayout, {showLink:false });


var x1a = [0,1,2,3,4,5,6,7,8,9,10];
var y1a = [1/11,1/11,1/11,1/11,1/11,1/16,1/16,1/16,1/16,1/16,1/16];
var x1b = [0,1,2,3,4];
var x1c = [5,6,7,8,9,10];
var y1b = [3/88,3/88,3/88,3/88,3/88];
var y1c = [5/176,5/176,5/176,5/176,5/176,5/176];
var ex1a    = { x: x1a, y: y1a, name: 'base', opacity: 0.75,  type: "bar" };
var ex1b    = { x: x1b, y: y1b, name: 'over', opacity: 0.75,  type: "bar" };
var ex1c    = { x: x1c, y: y1c, name: 'under', opacity: 0.75,  type: "bar" };
var ex1data = [ex1a,ex1b,ex1c];


var ex1layout = {
  bargap: 0.5, 
  bargroupgap: 0.2, 
  barmode: "stack", 
};

Plotly.newPlot('ex1', ex1data, ex1layout);

var plotData;

{
  var ptwo=1;
  var cX=[];
  var cY=[];
  for (var x=1; x<=128; x++) { 
    cX.push(x); 
	cY.push(2*(1+2*(x-Math.sqrt(x+x*x))));
	ptwo *= 2;
  }
  plotData = [{x:cX, y:cY, type:'scatter', mode:'markers', marker:{size:4}}];
}

var layout = {title: 'max bias', yaxis:{type: 'log', hoverformat:'e' } };
Plotly.newPlot('plot', plotData, layout, {showLink:false });

</script>
