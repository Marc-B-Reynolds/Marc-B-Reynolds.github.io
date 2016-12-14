---
layout:       post
title:        Viz of square to disc mappings
categories:   [math]
tags:         [distribution, sampling, map]
description:  a visualization of area/shape distortions and point distributions
plotly:       true
---

foo foo foo foo

<iframe width="600" height="400" src="https://www.shadertoy.com/embed/XldSWS?gui=true&t=10&paused=false&muted=false" frameborder="0" allowfullscreen></iframe>{: .center-image }

foo foo foo foo[^fong2015] [^fong2015a].

<br>

<script>

'use strict';

var stretchData, ellipticalData, squircularData, concentricData;
</script>

------

Radial Stretching <small></small>
------

ba

<div id="stretch" style="width:100%"></div>


<br>

------

Concentric <small>Shirley/Chiu modified by Cline and 'Franz'</small>
------

<div id="concentric" style="width:100%"></div>

<br>

------

Squircle <small></small>
------

<div id="squircular" style="width:100%"></div><br>

<br>

------

Elliptical <small></small>
------

<div id="elliptical" style="width:100%"></div><br>

<br>

------

foo <small></small>
------

<div id="graph" style="width:100%"></div><br>

<br>

------

References and Footnotes
------

[^fong2015]:*"Analytical Methods for Squaring the Disc"*, Chamberlain Fong, 2015. ([arXiv](https://arxiv.org/abs/1509.06344))


[^fong2015a]:*"Analytical Methods for Squaring the Disc (presentation slides)"*, Chamberlain Fong, 2015. ([arXiv](https://www.slideshare.net/chamb3rlain/analytical-methods))

[^lambers]:*"Mappings between Sphere, Disc, and Square"*, Martin Lambers, 2016. ([jcgt](https://jcgt.org/published/0005/02/01/))




<script>

var axisDef = { showgrid: false, showline: false, zeroline: false, range:[-1.2, 1.2], fixedrange: true };
var markDef = { size:4, color:'6565FF'};
var plotType = 'scatter';

// single: ulp(1.0)
var ulp1 = (1.0/68719476736.0);

function sgn(x) { if (x>=0.0) { return 1.0; } else { return -1.0; } }

// build a reference unit circle (100 points)
var unitCircle;

{
  var cX=[];
  var cY=[];

  for(var t=-Math.PI; t<=Math.PI; t+=(Math.PI/100)) {
    var x = Math.cos(t);
    var y = Math.sin(t);
    cX.push(x);
    cY.push(y);
  }
  
  unitCircle = {x:cX, y:cY, type:plotType, mode:'line', name:'unit circle', line: {color: 'C0C0FF'}};
}

// per dim
var pointSetInc  = (1.0/32.0);

// build a reference unit square
var unitSquare, squareX, squareY;
{
  squareX=[];squareY=[];

  var s = -1.0+pointSetInc;

  for(var i=s; i<=1.0; i += pointSetInc) {
    for(var j=s; j<=1.0; j += pointSetInc) {
      squareX.push(i);
      squareY.push(j);
	}
  }
  
  unitSquare = {x:squareX, y:squareY, type:plotType, mode:'markers', marker:markDef, name:'points'};
}


var stretch, stretchX, stretchY;
{
  var e = squareX.length;
  stretchX=[]; stretchY=[];
  
  for(var i=0; i<e; i++) {
    var x  = squareX[i]+ulp1;
    var y  = squareY[i];
	var x2 = x*x;
	var y2 = y*y;
	var n  = x2+y2;
	var t  = Math.sqrt(n);
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
	
    stretchX.push(s*rx);
    stretchY.push(s*ry);
  }
  stretch = {x:stretchX, y:stretchY, type:plotType, mode:'markers', marker:markDef, name:'points'};
}

var squircular, squircularX,squircularY;
{
  var e = squareX.length;
  squircularX=[]; squircularY=[];
  
  for(var i=0; i<e; i++) {
    var x  = squareX[i]+ulp1;
    var y  = squareY[i];
	var x2 = x*x;
	var y2 = y*y;
	var d  = x2+y2;
    var s  = Math.sqrt(1-x2*y2/d);
    squircularX.push(s*x);
    squircularY.push(s*y);
  }
  squircular = {x:squircularX, y:squircularY, type:plotType, mode:'markers', marker:markDef, name:'points'};
}

var eliptic, elipticX,elipticY;
{
  var e = squareX.length;
  elipticX=[]; elipticY=[];
  for(var i=0; i<e; i++) {
    var x  = squareX[i];
    var y  = squareY[i];
	var xs = Math.sqrt(1.0-0.5*y*y);
	var ys = Math.sqrt(1.0-0.5*x*x);
    elipticX.push(x*xs);
    elipticY.push(y*ys);
  }

  eliptic = {x:elipticX, y:elipticY, type:plotType, mode:'markers', marker:markDef, name:'points'};
}


var concentric, concentricX,concentricY;
{
  var e = squareX.length;
  concentricX=[]; concentricY=[];
  for(var i=0; i<e; i++) {
    var x  = squareX[i];
    var y  = squareY[i];
    var t,r;

    if (x*x> y*y) {
      r = x;
      t = (Math.PI/4.)*(y/x);
    } else {
      r = y;
      t = (Math.PI/2.0) - (Math.PI/4.0)*(x/y);
    }
    concentricX.push(r*Math.cos(t));
    concentricY.push(r*Math.sin(t));
  }

  concentric = {x:concentricX, y:concentricY, type:plotType, mode:'markers', marker:markDef, name:'points'};
}

// some hackiness to approximate square aspect ratio, since I couldn't figure it out.
var discLayout = {
  xaxis:  axisDef,
  yaxis:  axisDef,
  height: 600,
  width:  630,
  //aspectmode: 'data',
  hovermode: !1,
  updatemenus: [{
  buttons: [
    {method: 'animate', args: [['disc']],   label: 'disc'},
    {method: 'animate', args: [['square']], label: 'square'},
  ]
  }]
};

//var ease 

function to_square(g) { 
  Plotly.animate(g, {data: [unitCircle, unitSquare] }, 
  { transition: { duration: 2000, ease: 'cubic-in-out' }});
}

function to_disc(g,m) { 
  console.log('there');
  Plotly.animate(g, {data: [unitCircle, m] }, 
  { transition: { duration: 2000, ease: 'cubic-in-out'  }});
  console.log(concentricData);
}

function stretch_disc() { to_disc('squircular', squircularData); }


var options = {displaylogo: false, autosizable: true};

var frames0 = [{name: 'disc',   data: stretchData}, {name: 'square', data: unitSquare}];


var axisDef = { showgrid: false, showline: false, zeroline: false, range:[-1.2, 1.2], fixedrange: true };
var markDef = { size:4, color:'6565FF'};

var frames = [
  {name: 'square',     data: [{x: [], y: []}]},
  {name: 'stretch',    data: [{x: [], y: []}]},
  {name: 'squircle',   data: [{x: [], y: []}]},
  {name: 'concentric', data: [{x: [], y: []}]},
  {name: 'eliptic',    data: [{x: [], y: []}]},
];


// per dim
var pointSetInc  = (1.0/16.0);

var ulp1 = (1.0/68719476736.0);

function sgn(x) { if (x>=0.0) { return 1.0; } else { return -1.0; } }


// build a reference unit square
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

      var s  = Math.sqrt(1-x2*y2/d);
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

	 // eliptic
	 var xs = Math.sqrt(1.0-0.5*y*y);
	 var ys = Math.sqrt(1.0-0.5*x*x);
		 
    frames[4].data[0].x[i] = x*xs;
    frames[4].data[0].y[i] = y*ys;
		 
		i++;
	}
  }
  
  //unitSquare = {x:squareX, y:squareY, type:plotType, mode:'markers', marker:markDef, name:'points'};
}

stretchData      = [unitCircle, stretch ];
ellipticalData   = [unitCircle, eliptic ];
squircularData   = [unitCircle, squircular ];
concentricData   = [unitCircle, concentric ];

Plotly.newPlot('stretch',    stretchData,    discLayout, options);
Plotly.newPlot('elliptical', ellipticalData, discLayout, options);
Plotly.newPlot('squircular', squircularData, discLayout, options);
Plotly.newPlot('concentric', concentricData, discLayout, options);

Plotly.plot('graph', [{
  x: frames[0].data[0].x,
  y: frames[0].data[0].y,
  mode: 'markers',
  line: {simplify: false},
}], {
  xaxis: axisDef,
  yaxis: axisDef,
  height: 600,
  width:  630,
  updatemenus: [{
    buttons: [
	  {method: 'animate', args: [['square']],     label: 'square'},
      {method: 'animate', args: [['stretch']],    label: 'stretch'},
      {method: 'animate', args: [['squircle']],   label:  'squircle'},
	  {method: 'animate', args: [['concentric']], label:  'concentric'},
	  {method: 'animate', args: [['eliptic']],    label:  'eliptic'}
    ]
  }]
}).then(function() {
  Plotly.addFrames('graph', frames);
});


</script>
