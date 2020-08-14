---
layout:       post
title:        "Error growth of composing rotations"
tagline:      "it ain't much"
categories:   [quaternions]
tags:         [numeric,rotation]
description:  'Some quick-n-dirty empirical results to spitball the error growth on composing rotations'
plotly: true
---

\\
**DISCLAIMER:**  This is a half written blog post I started ages ago and got bored with. I'm tossing out as is since the topic has shown up on my timeline a couple of time recently. The TL;DR version is: worry about quat to matrix conversion and you don't *need* to normalize unless feeding to some method that expects the inputs to be normalized. shrug.

<br>

------

Problem statement and the set-up
------

\\
Some *set-up* bullet points:

* All quaternions (other than zero) represent a rotation.
* All quaternions that fall on a line through the origin represent the same rotation and the space of quaternion is that of infinitesimal rotations.
* The quaternion product is composition of rotations and the similarity transform $QPQ^{-1}$ applies a rotation $Q$ to a coordinate $P$.
* Assuming $Q$ is a unit quaternion we can reduce $QPQ^{-1}$ to $QPQ^{*}$ and if that assumption is false then the resulting transform is still the correct rotation composed with a uniform scaling of $Q \cdot Q$ (the magitude squared).
* The previous holds if converted into a matrix before application **provided** the diagonal components are computed without assuming unit magnitude. This is worth considering in all cases (SEE: [On quaternion/rotation matrix conversions and errors]({{site.base}}/quaternions/2017/08/08/QuatRotMatrix.html)).

What I am going to do is a little empirical demo with some graphs. Toy code can be found: [HERE](http://github.com/Marc-B-Reynolds/Stand-alone-junk/blob/master/src/Posts/quat_compose_error.c)

<br>

------

The test <small>quaternion style</small>
------

\\
The simple manufactured-solution test I'm going use is as follows. Note this is not attempting to pound on worse case input situations but to instead give a rough idea of what might actually occur in common cases.

1. Initialize quaternion $Q_0 = 1$
2. iteratively compose $n$ times: $ Q_{i+1} = R_i~Q_i $, where $R_i$ is a [uniform random rotation]({{site.base}}/distribution/2017/01/27/UniformRot.html).
3. iteratively compose $n$ times: $ Q_{i+1} = Q_i~\left(R_{i-n}\right)^* $ (the $R_j$ sequence is same as previous step) 

\\
As an example (in painful detail) if $n=4$ we have the following composition sequence:

$$ 
\begin{array}{llll}
Q_0 & = 1 \\
Q_1 & = R_0 ~ Q_0   & = R_0 \\
Q_2 & = R_1 ~ Q_1   & = R_1 ~ R_0 \\
Q_3 & = R_2 ~ Q_2   & = R_2 ~ R_1 ~ R_0 \\
Q_4 & = R_3 ~ Q_3   & = R_3 ~ R_2 ~ R_1 ~ R_0 \\
Q_5 & = Q_4 ~ R_0^* & = R_3 ~ R_2 ~ R_1 ~ R_0 ~ R_0^*          & = R_3 ~ R_2 ~ R_1\\
Q_6 & = Q_5 ~ R_1^* & = R_3 ~ R_2 ~ R_1 ~ R_0 ~ R_0^* ~ R_1^*  &=  R_3 ~ R_2\\
Q_7 & = Q_6 ~ R_2^* & = R_3 ~ R_2 ~ R_1 ~ R_0 ~ R_0^* ~ R_1^* ~ R_2^*  &= R_3 \\
Q_8 & = Q_7 ~ R_3^* & = R_3 ~ R_2 ~ R_1 ~ R_0 ~ R_0^* ~ R_1^* ~ R_2^* ~ R_3^* &= 1\\
\end{array}
$$

\\
so we've composed $2n-1$ products and the exact result is $Q_n = 1$.  We can measure how far we've deviated from unity by:

$$ \left|~ 1 - \sqrt{Q_n \cdot Q_n} ~ \right| $$

\\
Setting $n = 2^x$, performing computation in single precision (so $2^{x+1}+1$ products):


<div id="mag" style="width:100%"></div>


\\
So at the right end we have $x=16$ so we've composed $2\cdot2^{16}-1 = 131071$ products and the magnitude of the results is wavering around $ 1 \pm 8\cdot10^{-5} $ which would result in wanky uniform scale factor on $\left[0.999841, 1.00016\right]$. Given the $x$ axis is a log scale let's look again as a log plot:

<div id="lmag" style="width:100%"></div>

\\
So the error in magnitude roughly linearly increases (wrt the number of product). Now let's look at how the angle error (measured in degrees) behaves.

$$ \frac{360}{\pi} \text{atan}\left(\frac{\sqrt{U \cdot U}}{w} \right)$$

\\
where $U$ and $w$ are the bivector and scalar parts respectively. (recall an *extra* factor of 2 is needed to convert the angle measure from quaternion to 3D space)

<div id="angle"  style="width:100%"></div>
<div id="langle" style="width:100%"></div>


\\
First let me mention that using the standard methods for quaternion-to-matrix and matrix-to-quaternion introduces a peek round-trip error of ~0.000121 degrees. This is more than what we're seeing at $x=5$ which corresponds to a chain of 63 products. OK, like magnitude the angle error is roughly linearly increasing wrt the number of products. These plots also have some addition traces: 

| *unity*    |is as described for magnitude                                          |
| *tiny*     |starts with $Q_0 = 2^{-100}$ 											 |
| *huge*     |starts with $Q_0 = 2^{100}$											 |
| *norm*     |performs a normalization step after each product 						 |
| *no-fma*   |fp-contraction (auto FMA usage) is disabled    						 |
| *fma*      |the product is rewritten into eight $ab \pm cd$ terms in the style of Kahan's determinate (uses a *TwoProduct*) |
| *promote*  |same as *unity* except promotes to doubles, computes product and lowers to singles |
| *binary64* |same as *unity* except double precision quaternions (for later) |

\\
The point of *tiny*, *huge* and *norm* is to show that the magnitude of the inputs (and post normalizing) has almost no effect on the angle error. Notice that allowing auto & explict FMAs and promoting to double computation *do* significantly decrease the error.


<br>

------

The matrix test <small>I love LA</small>
------

\\
For matrices can just repeat the same process with a couple of tweaks. Swap quaternion to matrix product, the random quaternion is converted to matrix just prior to that (any angle error introduced by this step doesn't matter since the two terms which should cancel will be the same). To measure the angle error I should be using the fact that the trace is $1+2~\cos \theta$, but I'm being lazy and just converting back to a quaternion (we're spillballing and doing matrices as well came as an afterthought).

<div id="mangle"  style="width:100%"></div>
<div id="lmangle" style="width:100%"></div>

<br>

------

Revised test
------

\\
An obvious question would be: "But are uniform rotations a good choice?" especially if you've thought about what that means. Quickly: as a space very little of the total *volume* of rotations are those with a small rotation angle. Visualize rotating an object about all possible rotations that are approximately the indentity (yeah..it's pretty static).  Now attempt to visualize rotating about all possiable rotations that are approximately maximum angle (3D rotation of $\pi$). The [average uniform random rotation angle](http://marc-b-reynolds.github.io/quaternions/2017/11/10/AveRandomRot.html) works out to be $\frac{\pi}{2}+\frac{2}{\pi} \approx 126.476$ degrees. 


<div id="range_angle"  style="width:100%"></div>

<div id="lrange_angle"  style="width:100%"></div>

\\
So just about the same.

<br>

------

Some formal treatment
------

\\
I'm aware of a couple of papers with formal proofs WRT quaternion error bounds. Most recent is [*"Algorithms for manipulating quaternions in floating-point arithmetic"*](https://hal.archives-ouvertes.fr/hal-02470766) which covers some bounds on the norm, product and a specific quaternion-to-matrix conversion method. The other is [*"Verification of the Functional Behavior of a Floating-Point Program: an Industrial Case Study"*](https://hal.inria.fr/hal-00967124) which (seems to be more of a tutorial oriented paper of using some tools ([Frama-C](http://frama-c.com/), [Coq](https://coq.inria.fr/) & [Gappa](http://gappa.gforge.inria.fr/))) proves a (not tight) upper bound on the magnitude error. The example result is for in double precision, composing $10^{10} \left(\approx 2^{33.22}\right)$ rotations then the magnitude error is bound by: $\|q\| \leq 1.00003 $.

<br>

<script>

const options = {displaylogo: false, autosizable: true};
const gtDef   = {size:5, opacity:1.0};

function title(n)
{
 return { title: n, yaxis:{zeroline:false, hoverformat: 'g', type: 'linear', exponentformat: 'power'}, height: 500, width: 800 };
}

function title_log(n)
{
 return { title: n, yaxis:{zeroline:false, hoverformat: 'g', type: 'log', exponentformat: 'power'}, height: 500, width: 800 };
}

function make_trace(d,n) 
{
  return {x0:1, dx:1, y:d, name:n, type: 'scatter', mode:'lines+markers', marker: gtDef, opacity:0.7};
}

// standard dataset
const m_error   = [7.152557e-07,8.344650e-07,1.192093e-06,1.549721e-06,2.145767e-06,2.920628e-06,3.933907e-06,5.722046e-06,8.344651e-06,1.078844e-05,1.591444e-05,2.062321e-05,2.455713e-05,3.534554e-05,5.882976e-05,7.957213e-05];

const a_error = [2.670319e-05,3.258070e-05,5.133700e-05,6.712071e-05,1.021746e-04,1.430897e-04,1.738748e-04,2.458539e-04,3.137333e-04,4.493568e-04,6.069020e-04,9.099382e-04,1.179033e-03,1.642148e-03,2.564920e-03,3.403522e-03];


const a_error_x = [2.557049e-05,3.875041e-05,5.520555e-05,8.658650e-05,1.105201e-04,1.510674e-04,2.180266e-04,2.914399e-04,4.052511e-04,5.631373e-04,8.275389e-04,1.104928e-03,1.511751e-03,1.999411e-03,2.887302e-03,3.960702e-03];
const a_error_s = [2.670319e-05,3.258070e-05,5.133700e-05,6.712071e-05,1.021746e-04,1.430897e-04,1.738748e-04,2.458539e-04,3.137333e-04,4.866216e-04,6.369500e-04,9.663118e-04,1.172949e-03,1.716473e-03,2.337954e-03,3.493425e-03];
const a_error_h = [2.670319e-05,3.258070e-05,5.133700e-05,6.712071e-05,1.021746e-04,1.430897e-04,1.738748e-04,2.458539e-04,3.137333e-04,4.432504e-04,5.815096e-04,8.195559e-04,1.130657e-03,1.752940e-03,2.381121e-03,3.446459e-03];
const a_error_n = [2.617911e-05,3.923363e-05,4.916207e-05,6.690536e-05,9.958950e-05,1.314751e-04,1.730813e-04,2.506232e-04,3.439723e-04,4.750497e-04,6.148450e-04,9.225699e-04,1.159890e-03,1.740209e-03,2.400539e-03,3.869248e-03];
const a_error_f = [2.107959e-05,3.242128e-05,4.912591e-05,6.660270e-05,8.899992e-05,1.277100e-04,1.801008e-04,2.693142e-04,3.361507e-04,4.276641e-04,6.505829e-04,8.717225e-04,1.251805e-03,1.736884e-03,2.383090e-03,3.050364e-03];
const a_error_d = [3.935361e-14,5.986838e-14,9.617267e-14,1.265967e-13,1.846165e-13,2.397264e-13,3.218126e-13,4.352896e-13,6.512923e-13,8.322679e-13,1.114998e-12,1.587263e-12,2.331831e-12,3.198696e-12,4.259375e-12,6.848263e-12];
const a_error_p = [1.023900e-05,1.740074e-05,2.686347e-05,3.934514e-05,5.313407e-05,7.395189e-05,1.061775e-04,1.661962e-04,1.902090e-04,2.678548e-04,3.383249e-04,5.159365e-04,6.652062e-04,9.970340e-04,1.235722e-03,1.962197e-03];

const a_error_m0 = [8.913655e-06,1.304626e-05,2.036568e-05,2.756809e-05,3.611169e-05,5.151726e-05,7.135360e-05,9.814979e-05,1.423859e-04,1.957317e-04,2.551091e-04,3.445398e-04,5.320411e-04,7.677650e-04,1.034589e-03,1.277229e-03];
const a_error_m1 = [7.521756e-06,1.061353e-05,1.615059e-05,2.173843e-05,3.040226e-05,4.346609e-05,6.013318e-05,8.887199e-05,1.251394e-04,1.549613e-04,2.144815e-04,3.264943e-04,4.646432e-04,5.295614e-04,7.324356e-04,1.125635e-03];
const a_error_m2 = [8.517089e-06,1.193790e-05,1.779396e-05,2.410881e-05,3.567651e-05,4.731459e-05,6.351174e-05,9.036445e-05,1.169462e-04,1.782237e-04,2.400608e-04,3.358640e-04,4.291050e-04,6.725352e-04,8.421604e-04,1.150298e-03];

const a_error_a0 = [2.282713e-05,3.377459e-05,4.810045e-05,6.845824e-05,9.837508e-05,1.254243e-04,1.895656e-04,2.322295e-04,3.305645e-04,4.300595e-04,6.424422e-04,9.550615e-04,1.104680e-03,1.559405e-03,2.221777e-03,3.390861e-03];
const a_error_a1 = [2.282713e-05,3.377459e-05,4.810045e-05,6.845824e-05,9.837508e-05,1.254243e-04,1.895656e-04,2.322295e-04,3.305645e-04,4.300595e-04,6.424422e-04,9.550615e-04,1.104680e-03,1.559405e-03,2.221777e-03,3.390861e-03];
const a_error_a2 = [2.550368e-05,3.554277e-05,4.743884e-05,6.544314e-05,1.048058e-04,1.223550e-04,1.749512e-04,2.445619e-04,4.054647e-04,4.348241e-04,6.343069e-04,8.513650e-04,1.233221e-03,1.613091e-03,2.458394e-03,3.202742e-03];
const a_error_a3 = [2.488439e-05,3.468507e-05,4.652171e-05,6.769763e-05,9.487974e-05,1.264632e-04,1.758952e-04,2.602807e-04,3.484199e-04,4.523323e-04,6.043485e-04,8.642467e-04,1.186876e-03,1.714051e-03,2.418753e-03,3.274364e-03];
const a_error_a4 = [2.297942e-05,3.609799e-05,4.963255e-05,6.724680e-05,9.201753e-05,1.359564e-04,1.680106e-04,2.571109e-04,3.378848e-04,4.613654e-04,6.223112e-04,8.170242e-04,1.166652e-03,1.672521e-03,2.360407e-03,3.616338e-03];


Plotly.newPlot('mag',  [make_trace(m_error, 'm')], title("magnitude error"),        options);
Plotly.newPlot('lmag', [make_trace(m_error, 'm')], title_log("magnitude error (log)"), options);

Plotly.newPlot('angle', 
  [make_trace(a_error,   'unity'),
   make_trace(a_error_s, 'tiny'),
   make_trace(a_error_h, 'huge'),
   make_trace(a_error_n, 'norm'),
   make_trace(a_error_f, 'fma'),
   make_trace(a_error_x, 'no-fma'),
   make_trace(a_error_p, 'promote'),
   make_trace(a_error_d, 'binary64')
  ], 
  title("angle error"), options);

Plotly.newPlot('langle', 
  [make_trace(a_error,   'unity'),
   make_trace(a_error_s, 'tiny'),
   make_trace(a_error_h, 'huge'),
   make_trace(a_error_n, 'norm'),
   make_trace(a_error_f, 'fma'),
   make_trace(a_error_x, 'no-fma'),
   make_trace(a_error_p, 'promote'),
   make_trace(a_error_d, 'binary64')
  ], 
  title_log("angle error (log)"), options);


Plotly.newPlot('mangle', 
  [make_trace(a_error,    'quat (unity)'),
   make_trace(a_error_m0, 'matrix no-fma'),
   make_trace(a_error_m1, 'matrix'),
   make_trace(a_error_m2, 'matrix fma')
  ], 
  title("angle error"), options);

Plotly.newPlot('lmangle', 
  [make_trace(a_error,   'quat (unity)'),
   make_trace(a_error_m0, 'matrix no-fma'),
   make_trace(a_error_m1, 'matrix'),
   make_trace(a_error_m2, 'matrix fma')
  ], 
  title_log("angle error (log)"), options);

Plotly.newPlot('range_angle', 
  [make_trace(a_error_a0, 'pi/2'),
   make_trace(a_error_a1, 'pi/4'),
   make_trace(a_error_a2, 'pi/8'),
   make_trace(a_error_a3, 'pi/16'),
   make_trace(a_error_a4, 'pi/2048')
  ], 
  title("angle error"), options);

Plotly.newPlot('lrange_angle', 
  [make_trace(a_error_a0, 'pi/2'),
   make_trace(a_error_a1, 'pi/4'),
   make_trace(a_error_a2, 'pi/8'),
   make_trace(a_error_a3, 'pi/16'),
   make_trace(a_error_a4, 'pi/2048')
  ], 
  title_log("angle error (log)"), options);
