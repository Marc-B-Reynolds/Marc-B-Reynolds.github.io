---
layout:       post
title:        Square/Disc mappings
categories:   [math]
tags:         [quantization, distribution, sampling, map]
description:  a visualization of area/shape distortions and point distributions
plotly:       true
---

This note is for some auxilary information on some maps between the disc and square.  A fair amount of the math is convered by Fong[^fong2015] [^fong2015a] and the survey by Lambers[^lambers].  Specifically this contains:

* A pair of animated plots between the various maps 
* Show the Jacobian matrices and determinates and visualize with interactive plots
* Reduced complexity radial stretch in both directions
* Reduced complexity concentric disc to square
* Show a (perhaps) new radial approximate area preserving map

I will use the convention that $ \left(x,~y \right)$ is a coordinate on the square and $ \left(u,~v \right)$ a coordinate on the disc.

Define a signum like function:

$$
\text{sgn}\left(x\right) =
\begin{cases}
1  & x \geq 0 \\[2ex]
-1 & x < 0
\end{cases}
$$

<br>

------

Area distortion
------

\\
This shadertoy is intended to visualize how shapes and angles are transformed under the various maps.  Clicking through to the site allows interactivity.

<iframe width="600" height="400" src="https://www.shadertoy.com/embed/MtySRw?gui=true&t=10&paused=false&muted=false" frameborder="0" allowfullscreen></iframe>{: .center-image }

<br>

------

Point set
------

\\
As an alternate visualization this plot take uniform points on the square and shows how the various square to disc mappings transform the point set.

<div id="graph" style="width:100%"></div><br>

<br>

------

Radial Stretching
------

\\
Simplest square to disc map is to simply stretch (scale factor of $L_1$ over $L_2$ norm):

$$ \left(u,~v \right) = \frac{\text{max}\left( \abs{x},~\abs{y} \right)}{\sqrt{x^2+y^2}}\left(x,~y\right) $$

\\
Fong followed by Lambers provide an alternate formulation intended to be computationally friendly.  Note that we can reformulate the scale factor $s$ applied to the input coordinate as the following pseudo-code:

{% highlight c %}
float x2 = x*x;
float y2 = y*y;
float m  = x2 >= y2 ? x : y;
float s  = abs(m)*inversesqrt(x2+y2+epsilon);
// return s*(x,y)
{% endhighlight %}

\\
where $\text{epsilon}$ is some sufficiently small constant[^eps].  Formulated in this manner only needs a single select/cmov like operation and no branching to handle all cases including degenerate.

\\
A method to measure area distortion (growth/shrinkage) at infinitesimal around a point is to compute the Jacobian determinate[^jacobian] of the function. 

First we need the Jacobian, for $x^2 > y^2$:

$$
\left(x^2+y^2\right)^{-\frac{3}{2}} \left[
\begin{array}{cc}
 x^3+2 y^2 x & -x^2 y \\
 y^3         & x^3    \\
\end{array}
\right]
$$

\\
and for $x^2 \leq y^2$

$$
\left(x^2+y^2\right)^{-\frac{3}{2}}
\left[
\begin{array}{cc}
  y^3   & x^3 \\
 -x y^2 & y^3+2 x^2 y \\
\end{array}
\right]
$$

\\
The determinates of these are respectively:

$$
\frac{x^2}{x^2+y^2} \\
\frac{y^2}{x^2+y^2} 
$$

\\
Combining these into a single expression gives:

$$
\frac{ \max \left(x^2, y^2\right)}{x^2+y^2}
$$

\\
The area of the square is 4 and that of the disc is $\pi$, so the value where no area distortion occurs is $\frac{\pi}{4} \approx 0.785398$. Larger/smaller values are the local stretch/compressing respectively. The plot of the area distortion (determinate):

<div id="sheat" style="width:100%"></div>

\\
The disc to square transform is simply inverting the scale factor:

$$ \left(x,~y \right) = \frac{\sqrt{u^2+v^2}}{\text{max}\left( \abs{u},~\abs{v} \right)}\left(u,~v\right) $$

\\
and like above all cases can be handled with no branches and a single select.

<br>

------

Concentric
------

\\
Peter Shirley and Kenneth Chiu in 1997[^shirley1997] derived an area preserving map between square to disc.

\\
The square to disc map as given in Shirley's blog post[^shirley2011] (modifications noted by Dave Cline):

$$
\left(u,~v \right) =
\begin{cases}
\left(x \cos\left( \frac{\pi}{4} \frac{y}{x} \right), x \sin\left( \frac{\pi}{4} \frac{y}{x} \right) \right)  & \abs{x} \geq \abs{y} \\[2ex]
\left(y \cos\left( \frac{\pi}{2} - \frac{\pi}{4} \frac{x}{y} \right), y \sin\left( \frac{\pi}{2} - \frac{\pi}{4} \frac{x}{y} \right) \right)  & \abs{x} < \abs{y}
\end{cases}
$$

\\
Phase shift and pulling out signs:

$$
\left(u,~v \right) =
\begin{cases}
\left(x \cos\left( \frac{\pi}{4} \frac{y}{x} \right), x \sin\left( \frac{\pi}{4} \frac{y}{x} \right) \right)  & \abs{x} \geq \abs{y} \\[2ex]
\left(y \sin\left( \frac{\pi}{4} \frac{x}{y} \right), y \cos\left( \frac{\pi}{4} \frac{x}{y} \right) \right)  & \abs{x} < \abs{y}
\end{cases}
$$

\\
The more interesting form will be architecture/code specific. The second form trig inputs are on $\pm\frac{\pi}{4}$, so the cosine terms are always positive and removable without sign fixup and a narrow enough range for light weight approximation.

\\
The Jacobian for $\abs{x} \geq \abs{y}$:

$$
\left[
\begin{array}{cc}
 \cos \left(  \frac{\pi y}{4 x} \right) + \frac{\pi y}{4 x} \sin \left(\frac{\pi y}{4 x}\right) &
 -\frac{\pi}{4} \sin \left(\frac{\pi y}{4 x}\right) \\
 \sin \left(\frac{\pi y}{4 x}\right)-\frac{\pi y}{4 x} \cos \left(\frac{\pi y}{4 x}\right) &
 \frac{\pi}{4} \cos \left(\frac{\pi y}{4 x}\right) \\
\end{array}
\right]
$$

\\
The Jacobian for $\abs{x} < \abs{y}$:

$$
\left[
\begin{array}{cc}
 \frac{\pi}{4}  \cos \left(\frac{\pi}{4 y} (x+y) \right) & 
  \sin \left(\frac{\pi}{4 y} (x+y)\right)-\frac{\pi x}{4 y} \cos \left(\frac{\pi}{4 y}(x+y)\right)\\
 -\frac{\pi}{4}  \sin \left(\frac{\pi}{4 y} (x+y)\right) & 
  \cos \left(\frac{\pi}{4 y}(x+y)\right)+\frac{\pi x}{4 y} \sin \left(\frac{\pi}{4 y}(x+y)\right) \\
\end{array}
\right]
$$


\\
The determinant of both are constant $\left(\frac{\pi}{4} \right)$ as expected for an area preserving map.

The disc to square map:

$$
\left(x,~y \right) =
\begin{cases}
\text{sgn}(u)\sqrt{u^2+v^2}  \left(1, \frac{4}{\pi}\text{atan}\left(\frac{v}{u}\right) \right)   & u^2 > v^2 \\[2ex]
\text{sgn}(v)\sqrt{u^2+v^2}  \left(\frac{4}{\pi}\text{atan}\left(\frac{u}{v}\right) , 1\right)  & u^2 \leq v^2 
\end{cases}
$$

\\
Where $\text{atan}$ is single parameter.  This is degenerate when $v$ is approaching zero. This could be corrected by a select or since $\text{atan}$ is an odd function one possible rewrite choice is:

$$
\left(x,~y \right) =
\begin{cases}
\sqrt{u^2+v^2}  \left(\text{sgn}(u), \frac{4}{\pi}\text{atan}\left(\frac{v}{\abs{u}}\right) \right)  & u^2 > v^2 \\[2ex]
\sqrt{u^2+v^2}  \left(\frac{4}{\pi}\text{atan}\left(\frac{u}{\abs{v}+\epsilon}\right) , \text{sgn}(v)\right)  & u^2 \leq v^2 
\end{cases}
$$

\\
In either case the range of $\text{atan}$ in on $\pm 1$ so can be lightweight (for it) approximated.  ($\epsilon$ is some small constant[^eps] as above)

<br>

------

F.G. Squircle
------

\\
Fong[^fong2015] derivated mappings based on Fernandez-Guasti's squircle.

$$
\left(u,~v \right) = \frac{\sqrt{x^2+y^2 - x^2 y^2}}{\sqrt{x^2+y^2}} \left( x,y \right)
$$


\\
The Jacobian:

$$
\frac{1}{\left(x^2+y^2\right)^{3/2} \sqrt{y^2-x^2 \left(y^2-1\right)}} 

\left[
\begin{array}{cc}
 y^4-\left(y^2-1\right) x^4-2 y^2 \left(y^2-1\right) x^2 & -x^5 y \\
 -x y^5 & y^4+\left(1-2 y^2\right) x^4-y^2 \left(y^2-2\right) x^2 \\
\end{array}
\right]

$$

\\
The Jacobian determinate:

$$ 
1-\frac{2 x^2 y^2}{x^2+y^2}
$$

\\
the plot of the area distortion (determinate):

<div id="fheat" style="width:100%"></div>

\\
The disc to square map:

$$
n = u^2+v^2 \\
  \\
s = \frac{\text{sgn}(uv)}{\sqrt{2}}\sqrt{n-\sqrt{n\left(n-4u^2v^2\right)}}
$$

$$
\left(x,~y \right) = s \left(\frac{1}{u}, ~\frac{1}{v} \right)
$$


<br>

------

Elliptical
------

\\
The square to disc mapping was derived by Nowell in a blog post[^nowell]:

$$
\left(u,~v \right) = \left( x \sqrt{1-\frac{y^2}{2}}  , ~y \sqrt{1-\frac{x^2}{2}} \right)
$$

\\
The Jacobian:

$$
\left(
\begin{array}{cc}
 \sqrt{1-\frac{y^2}{2}} & -\frac{x y}{\sqrt{4-2 y^2}} \\
 -\frac{x y}{\sqrt{4-2 x^2}} & \sqrt{1-\frac{x^2}{2}} \\
\end{array}
\right)
$$

\\
The Jacobian determinate:

$$
\frac{2-\left(x^2+y^2\right)}{\sqrt{2-x^2} \sqrt{2-y^2}}
$$

\\
the plot of the area distortion (determinate):

<div id="eheat" style="width:100%"></div>

\\
Fong[^fong2015] provides two derivations of the inverse map, one of which:

$$
\left(x,~y \right) = \frac{1}{2}\left(\sqrt{2+t+2\sqrt{2}u}-\sqrt{2+t-2\sqrt{2}u} ,\sqrt{2-t+2\sqrt{2}v}-\sqrt{2-t-2\sqrt{2}v} \right)
$$

where:

$$ 
t = u^2-v^2
$$

<br>

------

Approximate equal area
------

\\
I haven't been able to find this method referenced anywhere and a twitter query came up negative as well.

The basic idea here is to form a low complexity middle ground between *radial stretch* and *concentric*. See the point set animation above. As such it is only of interest if computationally cheaper than concentric and higher performance is a priority.

\\
The disc to square map:

$$
\left(x,~y \right) =
\begin{cases}
\left( \text{sgn}(u) \sqrt{u^2+v^2} , ~\sqrt{2}v       \right)  & u^2 > v^2 \\[2ex]
\left( \sqrt{2}u      , ~\text{sgn}(v) \sqrt{u^2+v^2}  \right)  & u^2 \le v^2
\end{cases}
$$

\\
The square to disc map:

$$
\left(u,~v \right) =
\begin{cases}
\left( x \sqrt{1-\frac{y^2}{2x^2}} , ~\frac{y}{\sqrt{2}} \right)           & x^2 > y^2 \\[2ex]
\left( \frac{x}{\sqrt{2}}          , ~y \sqrt{1-\frac{x^2}{2y^2}} \right)  & x^2 \le y^2
\end{cases}
$$

\\
Which is degenerate when $y \le x $ and $y$ approaching zero.  Since the division term is positive we can simply add a bias as before. More likely to be interesting is we can multiply through by the denominator and get:

$$
\left(u,~v \right) =
\begin{cases}
\frac{1}{\sqrt{2}} \left(\text{sgn}(x)\sqrt{2x^2-y^2} , ~y  \right)  & x^2 > y^2 \\[2ex]
\frac{1}{\sqrt{2}} \left(x  , ~\text{sgn}(y)\sqrt{2y^2-x^2} \right)  & x^2 \le y^2
\end{cases}
$$

\\
or leave the half inside:

$$
\left(u,~v \right) =
\begin{cases}
 \left(\text{sgn}(x)\sqrt{x^2-\frac{y^2}{2}} , ~\frac{y}{\sqrt{2}}  \right)  & x^2 > y^2 \\[2ex]
 \left(\frac{x}{\sqrt{2}} , ~\text{sgn}(y)\sqrt{y^2-\frac{x^2}{2}} \right)  & x^2 \le y^2
\end{cases}
$$



\\
The Jacobian for $\abs{x} \geq \abs{y}$:

$$
\left(
\begin{array}{cc}
 \frac{1}{\sqrt{1-\frac{y^2}{2 x^2}}} & -\frac{y}{x \sqrt{4-\frac{2 y^2}{x^2}}} \\
 0 & \frac{1}{\sqrt{2}} \\
\end{array}
\right)
$$

\\
and the determinate:

$$ \frac{1}{\sqrt{2-\frac{y^2}{x^2}}} $$

\\
The Jacobian for $\abs{x} < \abs{y}$:

$$
\left(
\begin{array}{cc}
 \frac{1}{\sqrt{2}} & 0 \\
 -\frac{x}{\sqrt{4-\frac{2 x^2}{y^2}} y} & \frac{1}{\sqrt{1-\frac{x^2}{2 y^2}}} \\
\end{array}
\right)
$$

\\
and the determinate:

$$ \frac{1}{\sqrt{2-\frac{x^2}{y^2}}} $$

\\
and the heat map of area distortion. Recall equal area value is $\frac{\pi}{4} \approx 0.785398$

<div id="aheat" style="width:100%"></div>


<br>

------

Approximate area preserving (variant 2)
------

\\
In another post[^ballcube] we have a volume preserving map between cylinder and sphere which breaks the space into two parts: point inside and outside of the embedded conic. We can take a planar slice through the cylinder's axis to form a square to disc map.  Reducing the 3D equation to the slice gives:

$$
\left(u,~v \right) =
\begin{cases}
\left(x \sqrt{1-\frac{4}{9}\frac{y^2}{x^2}},~\frac{2}{3}y \right)  & x^2 \geq y^2 \\[2ex]
\left(x \sqrt{\frac{2}{3}-\frac{1}{9}\frac{x^2}{y^2}},~ y-\frac{1}{3}\frac{x^2}{y}\right)  & x^2  < y^2
\end{cases}
$$

\\
We can rewrite the $x^2 \geq y^2$ case to eliminate division (which can be zero) as:

$$
t = \frac{2}{3}y \\
\left( \text{sgn}\left(x\right) \sqrt{x^2-t^2}, ~t\right)
$$

\\
and re-express the $x^2<y^2$ case as:

$$
t = \frac{x}{3y} \\
\left( x\sqrt{\frac{2}{3}-t^2}, ~y-xt\right)
$$

\\
The Jacobian for for the cases are:

$$
\left(
\begin{array}{cc}
 \frac{3}{\sqrt{9-\frac{4 y^2}{x^2}}} &
 -\frac{4y}{3 x \sqrt{9-\frac{4 y^2}{x^2}}} \\
 0 &
 \frac{2}{3} \\
\end{array}
\right)
$$

$$
\left(
\begin{array}{cc}
 \frac{6y^2-2x^2}{3\sqrt{6-y^2\frac{x^2}{y^2}}} & \frac{x^3}{3y^3\sqrt{6-\frac{x^2}{y^2}}} \\
 -\frac{2x}{3y} & 1+\frac{x^2}{3y^2} \\
\end{array}
\right)
$$

\\
with the determinates:

$$
\frac{2}{\sqrt{9-\frac{4y^2}{x^2}}} \\
\frac{2}{\sqrt{6-\frac{x^2}{y^2}}} 
$$

\\
The area distortion:

<div id="bheat" style="width:100%"></div>

\\
For the disc to square map, start with a common term:

$$
t = \sqrt{u^2+v^2}
$$

\\
then the map is:

$$
\left(x,~y \right) =
\begin{cases}
\left(u \sqrt{\frac{3t}{t+\abs{v}}}, ~\text{sgn}\left(v\right)t \right)  & \frac{5}{4}v^2 \geq u^2 \\[2ex]
\left(\text{sgn}\left(u\right)t,~ \frac{3}{2}v\right)  & \frac{5}{4}v^2 < u^2
\end{cases}
$$


<br>

------

References and Footnotes
------

[^eps]:      A reasonable(ish) epsilon for the funcs here under generic usage is the smallest normal fp value (single: $2^{-126}$).

[^jacobian]: **"Jacobian matrix and determinant"**, Wikipedia. ([link](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant))

[^shirley1997]: **"A Low Distortion Map Between Disk and Square"**, Peter Shirley, Kenneth Chiu, 1997. ([PDF](http://pdfs.semanticscholar.org/4322/6a3916a85025acbb3a58c17f6dc0756b35ac.pdf))

[^shirley2011]: **"Improved code for concentric map"**, Peter Shirley's blog, 2011. ([link](http://psgraphics.blogspot.fr/2011/01/improved-code-for-concentric-map.html))

[^fong2015]: **"Analytical Methods for Squaring the Disc"**, Chamberlain Fong, 2015. ([PDF](http://arxiv.org/abs/1509.06344))

[^fong2015a]: **"Analytical Methods for Squaring the Disc (presentation slides)"**, Chamberlain Fong, 2014. ([link](http://www.slideshare.net/chamb3rlain/analytical-methods))

[^lambers]:  **"Mappings between Sphere, Disc, and Square"**, Martin Lambers, 2016. ([link](http://jcgt.org/published/0005/02/01/))

[^nowell]:  **"Mapping a Square to a Circle"**, Philip Nowell, 2005. ([link](http://mathproofs.blogspot.fr/2005/07/mapping-square-to-circle.html))

[^ballcube]: **Cube/cylinder/ball mappings**, ([post]({{site.base}}/math/2017/01/27/CubeBall.html))


<script>

'use strict';

// single: ulp(1.0)
var ulp1    = (1.0/68719476736.0);
var sqrt2o2 = 0.5*Math.sqrt(2.0);

var options = {displaylogo: false, autosizable: true};

var size = 1024;
var dSize = 2.0/size;
var sheatZ = new Array(size);
var eheatZ = new Array(size);
var fheatZ = new Array(size);
var aheatZ = new Array(size);
var bheatZ = new Array(size);

for(var i=0; i<size; i++) {
 sheatZ[i] = new Array(size);
 eheatZ[i] = new Array(size);
 fheatZ[i] = new Array(size);
 aheatZ[i] = new Array(size);
 bheatZ[i] = new Array(size);
}

for(var j=0; j<size; j++) {
  var y  = j*dSize-1.0;
  var y2 = y*y;
  for(var i=0; i<size; i++) {
    var x    = i*dSize-1.0;
    var x2   = x*x;
    var norm = x2+y2;
    sheatZ[j][i] = Math.max(x2,y2)/(norm+ulp1);
    eheatZ[j][i] = (2.0-norm)/(Math.sqrt(2.0-x2)*Math.sqrt(2.0-y2)+ulp1);
    fheatZ[j][i] = 1.0-(2.0*x2*y2)/(norm+ulp1);
	
	if (x2 > y2) {
      aheatZ[j][i] = 1.0/Math.sqrt(2.0-y2/x2);
      bheatZ[j][i] = 2.0/Math.sqrt(9.0-4*y2/x2);
    } else {
      aheatZ[j][i] = 1.0/Math.sqrt(2.0-x2/(y2+ulp1));
      bheatZ[j][i] = 2.0/Math.sqrt(6.0-x2/(y2+ulp1));
    }
  }
}

var heatColor;

var Blackbody =
[
  [0,   'rgb(0,0,0)'], 
  [0.2, 'rgb(230,0,0)'],
  [0.4, 'rgb(230,210,0)'],
  [0.78539816, 'rgb(255,255,255)'],
  [1,   'rgb(160,200,255)']
];

var BlackbodyH =
[
  [0,   'rgb(0,0,0)'], 
  [0.2, 'rgb(230,0,0)'],
  [0.4, 'rgb(230,210,0)'],
  [0.78539816, 'rgb(255,255,255)'],
  [1,   'rgb(160,200,255)']
];

var Bar = [
  [0.,        'rgb(128,  0, 38)'],
  [0.0981748, 'rgb(189,  0, 38)'],
  [0.19635,   'rgb(227, 26, 28)'],
  [0.294524,  'rgb(252, 78, 42)'],
  [0.392699,  'rgb(253,141, 60)'],
  [0.490874,  'rgb(254,178, 76)'],
  [0.589049,  'rgb(254,217,118)'],
  [0.78539816,'rgb(255,255,255)'],
  [0.812223,  'rgb(199,233,180)'],
  [0.839049,  'rgb(127,205,187)'],
  [0.865874,  'rgb( 65,182,196)'],
  [0.892699,  'rgb( 29,145,192)'],
  [0.919524,  'rgb( 34, 94,168)'],
  [0.94635,   'rgb( 37, 52,148)'],
  [1,         'rgb(  8, 29, 88)'],
];

//heatColor = 'YIGnBu';
heatColor = BlackbodyH;
//heatColor = Bar;

function makeDataUS(set)
{
  return [{
  z: set,
  colorscale: heatColor,
  type:  'heatmap',
  xtype: 'scaled',
  x0:    -1,
  dx:     2/size,
  ytype: 'scaled',
  y0:    -1,
  dy:     2/size,
  }
]; 
}

function makeData(set)
{
  return [{
  z: set,
  colorscale: heatColor,
  type:  'heatmap',
  zmax:   1,
  zmin:   0,
  xtype: 'scaled',
  x0:    -1,
  dx:     2/size,
  ytype: 'scaled',
  y0:    -1,
  dy:     2/size,
  }
]; 
}

var sheatData = makeData(sheatZ);
var eheatData = makeData(eheatZ);
var fheatData = makeData(fheatZ);
var aheatData = makeData(aheatZ);
var bheatData = makeData(bheatZ);

var heatLayout = { height: 600, width: 620 };

Plotly.newPlot('sheat', sheatData, heatLayout, options);
Plotly.newPlot('eheat', eheatData, heatLayout, options);
Plotly.newPlot('fheat', fheatData, heatLayout, options);
Plotly.newPlot('aheat', aheatData, heatLayout, options);
Plotly.newPlot('bheat', bheatData, heatLayout, options);

var axisDef = { showgrid: false, showline: false, zeroline: false, range:[-1.2, 1.2], fixedrange: true };
var markDef = { size:4, color:'6565FF'};
var plotType = 'scatter';

var options = {displaylogo: false, autosizable: true};

var frames = [
  {name: 'square',     data: [{x: [], y: []}]},
  {name: 'stretch',    data: [{x: [], y: []}]},
  {name: 'squircle',   data: [{x: [], y: []}]},
  {name: 'concentric', data: [{x: [], y: []}]},
  {name: 'elliptic',   data: [{x: [], y: []}]},
  {name: 'aea',        data: [{x: [], y: []}]},
  {name: 'aea2',       data: [{x: [], y: []}]},
];


// per dim
var pointSetInc  = (1.0/20.0);


function sgn(x) { if (x>=0.0) { return 1.0; } else { return -1.0; } }

function buildPointSets()
{
  var S = -1.0+pointSetInc/2.0;
  var i = 0;
  for(var x=S; x<=1.0; x += pointSetInc) {
    for(var y=S; y<=1.0; y += pointSetInc) {
      // square
      frames[0].data[0].x[i] = x;
      frames[0].data[0].y[i] = y;
        
       // common stuff
       var x2 = x*x;
       var y2 = y*y;
       var d  = x2+y2;
       var t  = Math.sqrt(d);
       var xy = x*y;
       var rx,ry,s;
    
       if (x2 >= y2) {
         s  = sgn(x)/t;
         rx = x2;
         ry = xy;
       }
       else {
         s  = sgn(y)/t;
         rx = xy;
         ry = y2;
       }
    
       frames[1].data[0].x[i] = s*rx;
       frames[1].data[0].y[i] = s*ry;

       s  = Math.sqrt(1-x2*y2/d);
       frames[2].data[0].x[i] = s*x;
       frames[2].data[0].y[i] = s*y;

       // concentric
       if (x2 > y2) {
        s = x;
        t = (Math.PI/4.)*(y/x);
       } else {
        s = y;
        t = (Math.PI/2.0) - (Math.PI/4.0)*(x/y);
       }
         
      frames[3].data[0].x[i] = s*Math.cos(t);
      frames[3].data[0].y[i] = s*Math.sin(t);

      // elliptic
      var xs = Math.sqrt(1.0-0.5*y*y);
      var ys = Math.sqrt(1.0-0.5*x*x);
         
      frames[4].data[0].x[i] = x*xs;
      frames[4].data[0].y[i] = y*ys;
	  
      // approx equal area
      if (x2 > y2) {
	    var t = 2/3*y;
        frames[5].data[0].x[i] = x*Math.sqrt(1.0-.5*(y2/x2));
        frames[5].data[0].y[i] = y*sqrt2o2;
        frames[6].data[0].x[i] = sgn(x)*Math.sqrt(x2-t*t);
        frames[6].data[0].y[i] = t;
       } else {
	    var t = x/(3*y);
        frames[5].data[0].x[i] = x*sqrt2o2;
        frames[5].data[0].y[i] = y*Math.sqrt(1.0-.5*(x2/(y2+ulp1)));
        frames[6].data[0].x[i] = x*Math.sqrt(2/3-t*t);
        frames[6].data[0].y[i] = y-x*t;
       }

      i++;
    }
  }
}

buildPointSets();


Plotly.plot('graph', [{
  x: frames[0].data[0].x,
  y: frames[0].data[0].y,
  mode: 'markers',
  marker: markDef,
  line: {simplify: false},
}], {
  xaxis: axisDef,
  yaxis: axisDef,
  height: 600,
  width:  630,
  displaylogo: false,
  updatemenus: [{
    buttons: [
      {method: 'animate', args: [['square']],     label: 'square'},
      {method: 'animate', args: [['stretch']],    label: 'stretch'},
      {method: 'animate', args: [['concentric']], label: 'concentric'},
      {method: 'animate', args: [['aea']],        label: 'aea'},
      {method: 'animate', args: [['aea2']],       label: 'aea2'},
      {method: 'animate', args: [['squircle']],   label: 'squircle'},
      {method: 'animate', args: [['elliptic']],   label: 'elliptic'}
    ]
  }]
}).then(function() {
  Plotly.addFrames('graph', frames);
});


</script>
