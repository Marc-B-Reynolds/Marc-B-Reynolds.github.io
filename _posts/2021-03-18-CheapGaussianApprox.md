---
layout:       post
title:        A cheap normal distribution approximation
categories:   [distribution]
tags:         [random]
description:  A brief explaination and implementation of the standard normal distribution approximation "on the cheap".
plotly:       true
---

Just about a year ago I tossed out a [standard normal distribution](https://en.wikipedia.org/wiki/Normal_distribution) approximation in a tweet just for the fun factor. I ran across it recently and noticed I could make a couple of improvements. Per Vognsen thoughfully broke down the "how it works" in a pair of tweets:

<blockquote class="twitter-tweet"><p lang="en" dir="ltr">BTW, an intuitive way to think about this is to remember that addition of random variables corresponds to convolution of their probability densities. So adding the two uniform variates is equivalent to box filtering the stair-stepped distribution twice.</p>&mdash; Per Vognsen (@pervognsen) <a href="https://twitter.com/pervognsen/status/1240900754052763649?ref_src=twsrc%5Etfw">March 20, 2020</a></blockquote> <script async src="https://platform.twitter.com/widgets.js" charset="utf-8"></script> 

<br>

------


Baseline version
------

\\
Before jumping into the improvements let's take a look at an alternate approximation for comparison purposes. The super old skool method of approximating a normal distriubtion is to simply sum up a bunch of uniform distributions (so also [central limit theorem](https://en.wikipedia.org/wiki/Central_limit_theorem) based). Since our *population count* version draws two 64-bit uniform numbers let's stick with that to keep the uniform draw cost the same.

{% highlight c %}
float dist_normal_approx_sum()
{
  // Generate two 64-bit uniform random integers. These could
  // be passed in as parameters and/or the two values could
  // be drawn from independent generators to break the serial
  // dependency.
  uint64_t u0 = rng_u64();
  uint64_t u1 = rng_u64();
  
  // partition both into two 32-bit uniform integers each
  int64_t  a  = (int64_t)(u0 & 0xffffffff);
  int64_t  b  = (int64_t)(u0 >> 32);
  int64_t  c  = (int64_t)(u1 & 0xffffffff);
  int64_t  d  = (int64_t)(u1 >> 32);

  // magic constant (not really)
  const float K = 0x1.b566e2p-32f;  // ~1.70860*2^-32

  // perform the sum and scale. The pairings can be in any
  // if that helps with throughput. Flipping the two adds
  // to subs and the one sub to add is equivalent.
  return K*((a+b)-(c+d));
}
{% endhighlight %}

\\
As for the scaling value $K$ we need to include the factor to take our fixed-point integer to floating point which is $2^{-32}$.   We can look up the variance for the sum of $n$ uniform values which is $\frac{n}{12}$ and the standard deviation is the square root of the variance.  Since the standard normal distribution we need an additional factor of the reciprocal of our sum's standard deviation: $\sqrt{3} \approx 1.732051$. The value of $K$ in the code is different from this because I ran an optimization to *minimize the maximum error* (for absolute error).

<br>

------


Improved popcount version
------

\\
The modified population count version:


{% highlight c %}
float dist_normal_approx()
{
  // Generate two 64-bit uniform random integers. These could
  // be passed in as parameters and/or the two values could
  // be drawn from independent generators. This would allow
  // to break the serial dep on u0 & u1.
  uint64_t u0 = rng_u64();

  // We compute population count (aka hamming weight) of u0.
  // this gives us a binomial distribution (p=1\2, cut 64).
  // The subtraction centers the distribution on [-32,32].
  int64_t  bd = (int64_t)__builtin_popcountl(u0)-32;
  
  // draw & partition u1 into two 32-bit integers
  uint64_t u1 = rng_u64();
  int64_t  a  = (int64_t)(u1 & 0xffffffff);
  int64_t  b  = (int64_t)(u1 >> 32);

  // triangle distribution
  int64_t  td = a-b;                      // <-- changed add to sub

  // sum the binomial and triangle distributions
  float    r  = (float)((bd<<32) + td);   // <-- nuked a shift here

  // scale the result 
  return r * 0x1.fb760cp-35f;             // <-- nuked a constant add here & magic constant!
}

{% endhighlight %}

\\
There are three changes (marked as side comments):

1. The triangle distribution uses a subtraction instead of addition. This makes it centered around zero which removes the final subtraction (-0.25f) in the original version. This reduces the dependency chain length by one.
2. The original was (foolishly) designed to have no rounding when converting the integer to floating point. That's a non feature here so the left shift has been adjusted up and the right shift is eliminated.
3. In addition to nuking the subtraction (as per item 1) the scaling constant has been *minimax* optimized like the *summation* method above.

So the cost of the transform is now ~9 primitive ops (dropping two) on 64-bit hardware with native popcount.

<br>

------


32-bit popcount version
------

\\
In the case of hardware that either has 32-bit popcount (but not 64-bit) or the population count must be performed in software it seems worthwhile to toss out a variant for those cases. Adjust the centering of the triangle distribution and create new search for a new scaling constant...done:

{% highlight c %}
// in graphs: pop32
float dist_normal_approx()
{
  uint64_t u0 = rng_u64();
  int64_t  bd = (int64_t)__builtin_popcount((uint32_t)u0)-16;
  uint64_t u1 = rng_u64();
  int64_t  a  = (int64_t)(u1 & 0xffffffff);
  int64_t  b  = (int64_t)(u1 >> 32);
  int64_t  td = a-b;
  float    r  = (float)((bd<<31) + td);

  return r *  0x1.59db68p-33f;
}
{% endhighlight %}

\\
A sadface about the previous version is we're tossing away 32 perfectly good random bits. It'll cost us a couple of operations but we can instead use them as another convolution step:

{% highlight c %}
// in graphs: pop32x 
float dist_normal_approx()
{
  uint64_t u0 = rng_u64();
  int64_t  bd = (int64_t)__builtin_popcount((uint32_t)u0)-16;
  uint64_t u1 = rng_u64();
  int64_t  a  = (int64_t)(u1 & 0xffffffff);
  int64_t  b  = (int64_t)(u1 >> 32);
  int64_t  td = a-b+((int64_t)u0>>UINT64_C(32));
  float    r  = (float)((bd<<31) + td);

  return r *  0x1.540aep-33f;
}
{% endhighlight %}


<br>

------


Comparisons and comments
------

\\
Let me first state that the summation method is probably sufficient for most "horseshoes and handgrenades" usages.  The *"[Encyclopedia of Mathematics](https://encyclopediaofmath.org/wiki/Uniform_distribution)"* goes so far to say: "the approximation for $n=3$ is already satisfactory for many practical purposes".

The most commonly used method to generate a normal distribution is to use the [Box-Muller transform](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform) which generates two uniform random numbers and transforms with a logarithm, square root and a sin/cos pair to produce two results.

Although the scaling constants have been optimized they are not optimal. My search/testing method was a fast hack. They are, however, nearly optimal (famous last words) because the peak errors of each swing nearly to the same peak values.

{: .center }
<div id="sum" style="width:100%"></div><br>


{: .center }
<div id="err" style="width:100%"></div><br>


{: .center }
|method     | ops | range | error |
|:---       |:---:| :---: |:---: |
|sum        |  8  | $\pm3.4172022$ | $8.898866 \times 10^{-3}$  |
|pop        |  9  | $\pm8.1768640$ | $9.249441 \times 10^{-4}$  |
|pop32      | 10  | $\pm6.079518$  | $ 2.213490\times 10^{-3}$  |
|pop32x     | 12  | $\pm4.611270$ | $ 1.391753\times 10^{-3}$  |
|box-muller | --  | $\pm5.768108 / \pm8.571674$ | --  |

<br>

* ops is the number of CPU operations ops needed for the transform on x64
* range is the range of the output values. The Box-Muller numbers are for single/double precision respectively.
* error is peak absolute error

<script>


var add = function(a, b){ return a + b; };
var sub = function(a, b){ return a - b; };
var mul = function(a, b){ return a * b; };
var div = function(a, b){ return a / b; };
var rer = function(a, b){ return a / b - 1.0; };


function array_op(A,B,F) { return B.map(function(B,i) { return F(A[i],B); }); }

const ref = [1.516364e-04,1.958286e-04,2.480530e-04,3.136446e-04,3.964780e-04,4.985649e-04,6.257096e-04,7.851323e-04,9.742326e-04,1.206519e-03,1.491380e-03,1.831372e-03,2.239637e-03,2.744640e-03,3.327692e-03,4.033016e-03,4.876170e-03,5.871608e-03,7.002380e-03,8.349732e-03,9.907851e-03,1.173731e-02,1.381549e-02,1.624693e-02,1.897255e-02,2.207990e-02,2.560267e-02,2.960182e-02,3.407681e-02,3.903184e-02,4.459711e-02,5.073967e-02,5.747823e-02,6.485610e-02,7.296069e-02,8.169752e-02,9.108779e-02,1.012288e-01,1.120722e-01,1.235783e-01,1.357357e-01,1.485096e-01,1.618004e-01,1.756016e-01,1.898665e-01,2.045063e-01,2.193840e-01,2.344249e-01,2.495362e-01,2.645291e-01,2.794713e-01,2.939928e-01,3.081636e-01,3.217034e-01,3.344638e-01,3.464527e-01,3.573340e-01,3.673223e-01,3.759834e-01,3.834218e-01,3.893615e-01,3.940416e-01,3.971954e-01,3.986207e-01,3.988255e-01,3.972333e-01,3.940350e-01,3.895744e-01,3.833662e-01,3.760379e-01,3.672841e-01,3.574198e-01,3.464841e-01,3.344293e-01,3.216952e-01,3.081482e-01,2.939254e-01,2.794463e-01,2.645187e-01,2.496151e-01,2.344244e-01,2.193671e-01,2.044715e-01,1.899161e-01,1.755403e-01,1.617916e-01,1.484702e-01,1.356634e-01,1.235113e-01,1.120898e-01,1.012232e-01,9.104664e-02,8.165530e-02,7.293026e-02,6.486913e-02,5.748405e-02,5.072847e-02,4.462532e-02,3.903323e-02,3.408638e-02,2.959889e-02,2.557674e-02,2.208696e-02,1.895728e-02,1.623120e-02,1.381721e-02,1.173703e-02,9.907635e-03,8.350238e-03,7.014831e-03,5.845463e-03,4.862021e-03,4.053447e-03,3.336112e-03,2.750124e-03,2.235271e-03,1.829360e-03,1.490471e-03,1.207160e-03,9.781221e-04,7.849759e-04,6.245324e-04,5.008076e-04,3.976850e-04,3.136297e-04,2.459965e-04,1.930493e-04,1.521207e-04,1.521207e-04];

const sum = [0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,2.443790e-07,9.167194e-06,5.193949e-05,1.581371e-04,3.548324e-04,6.699026e-04,1.131409e-03,1.775354e-03,2.612495e-03,3.702295e-03,5.046606e-03,6.670487e-03,8.616418e-03,1.093525e-02,1.361244e-02,1.669770e-02,2.023021e-02,2.420340e-02,2.869787e-02,3.367322e-02,3.925070e-02,4.539006e-02,5.211847e-02,5.943418e-02,6.749890e-02,7.624607e-02,8.567292e-02,9.580989e-02,1.067893e-01,1.184696e-01,1.307372e-01,1.435354e-01,1.569161e-01,1.705699e-01,1.845005e-01,1.986378e-01,2.128828e-01,2.271220e-01,2.412823e-01,2.553942e-01,2.691811e-01,2.826668e-01,2.958069e-01,3.084136e-01,3.203961e-01,3.317557e-01,3.424185e-01,3.520955e-01,3.609921e-01,3.687957e-01,3.756003e-01,3.810536e-01,3.854111e-01,3.884972e-01,3.898449e-01,3.898425e-01,3.884443e-01,3.854532e-01,3.812521e-01,3.755914e-01,3.687241e-01,3.610601e-01,3.520663e-01,3.423627e-01,3.318208e-01,3.204743e-01,3.083590e-01,2.957803e-01,2.826315e-01,2.691858e-01,2.553413e-01,2.412929e-01,2.271073e-01,2.128339e-01,1.986499e-01,1.844268e-01,1.704987e-01,1.568915e-01,1.436317e-01,1.307877e-01,1.184602e-01,1.068185e-01,9.585915e-02,8.567011e-02,7.620355e-02,6.750749e-02,5.944643e-02,5.210643e-02,4.532149e-02,3.924521e-02,3.364754e-02,2.866204e-02,2.421386e-02,2.020907e-02,1.669668e-02,1.362515e-02,1.092229e-02,8.628482e-03,6.684065e-03,5.055350e-03,3.698462e-03,2.620113e-03,1.777339e-03,1.136845e-03,6.678939e-04,3.520846e-04,1.568019e-04,5.155802e-05,8.922815e-06,2.503395e-07,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00];

const pop = [1.102734e-04,1.529224e-04,2.028494e-04,2.516022e-04,3.078774e-04,4.223829e-04,5.390939e-04,6.699974e-04,8.091505e-04,1.074015e-03,1.348882e-03,1.628339e-03,1.969257e-03,2.536688e-03,3.136146e-03,3.728577e-03,4.436675e-03,5.578196e-03,6.755349e-03,7.938458e-03,9.299673e-03,1.140002e-02,1.356340e-02,1.572692e-02,1.821188e-02,2.179932e-02,2.542800e-02,2.911052e-02,3.321725e-02,3.887953e-02,4.459991e-02,5.035033e-02,5.673773e-02,6.491404e-02,7.322310e-02,8.152615e-02,9.065641e-02,1.016296e-01,1.126776e-01,1.237251e-01,1.355591e-01,1.490582e-01,1.625288e-01,1.760126e-01,1.900386e-01,2.050912e-01,2.201806e-01,2.351240e-01,2.500982e-01,2.650253e-01,2.799241e-01,2.948060e-01,3.089172e-01,3.216683e-01,3.344630e-01,3.472480e-01,3.582297e-01,3.668678e-01,3.755961e-01,3.841132e-01,3.901758e-01,3.932781e-01,3.962864e-01,3.995017e-01,3.993343e-01,3.963439e-01,3.933280e-01,3.901930e-01,3.842332e-01,3.754738e-01,3.669582e-01,3.582016e-01,3.472324e-01,3.344691e-01,3.217110e-01,3.088181e-01,2.948521e-01,2.799161e-01,2.650281e-01,2.501043e-01,2.351843e-01,2.201167e-01,2.050654e-01,1.900767e-01,1.760737e-01,1.625338e-01,1.489952e-01,1.355742e-01,1.237203e-01,1.126560e-01,1.016368e-01,9.061501e-02,8.152077e-02,7.322668e-02,6.488979e-02,5.673100e-02,5.032578e-02,4.462168e-02,3.888133e-02,3.323714e-02,2.909544e-02,2.542226e-02,2.176047e-02,1.821135e-02,1.572185e-02,1.357310e-02,1.140256e-02,9.298147e-03,7.933803e-03,6.753817e-03,5.566101e-03,4.438237e-03,3.729906e-03,3.130138e-03,2.526310e-03,1.964899e-03,1.637948e-03,1.355982e-03,1.069193e-03,8.084710e-04,6.660753e-04,5.427359e-04,4.226154e-04,3.089205e-04,2.505352e-04,2.018242e-04,1.538821e-04,1.087236e-04,1.087236e-04];

const pop32 = [9.188919e-05,1.286961e-04,1.612651e-04,1.967607e-04,2.382169e-04,3.514630e-04,4.782518e-04,6.100057e-04,7.351672e-04,8.709922e-04,1.145784e-03,1.556162e-03,1.965395e-03,2.371684e-03,2.778837e-03,3.326900e-03,4.378100e-03,5.460695e-03,6.556833e-03,7.655903e-03,8.775740e-03,1.087372e-02,1.339125e-02,1.588618e-02,1.839611e-02,2.090991e-02,2.427716e-02,2.915327e-02,3.410633e-02,3.902810e-02,4.397682e-02,4.923154e-02,5.702515e-02,6.538387e-02,7.376509e-02,8.211781e-02,9.048099e-02,1.007451e-01,1.129034e-01,1.250285e-01,1.371634e-01,1.492104e-01,1.618241e-01,1.764657e-01,1.912690e-01,2.060341e-01,2.208038e-01,2.356308e-01,2.504148e-01,2.653031e-01,2.800264e-01,2.947995e-01,3.096315e-01,3.232020e-01,3.344351e-01,3.454579e-01,3.565710e-01,3.677475e-01,3.781543e-01,3.835545e-01,3.876902e-01,3.919179e-01,3.959604e-01,4.000780e-01,4.001331e-01,3.960177e-01,3.918215e-01,3.877378e-01,3.834797e-01,3.782292e-01,3.676171e-01,3.566315e-01,3.454869e-01,3.343887e-01,3.232564e-01,3.096992e-01,2.948591e-01,2.800719e-01,2.651887e-01,2.504195e-01,2.356457e-01,2.208796e-01,2.060161e-01,1.912469e-01,1.764076e-01,1.617874e-01,1.492519e-01,1.371487e-01,1.250150e-01,1.129428e-01,1.008542e-01,9.049097e-02,8.211626e-02,7.377213e-02,6.538110e-02,5.704622e-02,4.922940e-02,4.400359e-02,3.901587e-02,3.411379e-02,2.915129e-02,2.426239e-02,2.092220e-02,1.840457e-02,1.589090e-02,1.338277e-02,1.087750e-02,8.796012e-03,7.658812e-03,6.554323e-03,5.464957e-03,4.371651e-03,3.319879e-03,2.776626e-03,2.370391e-03,1.956669e-03,1.557705e-03,1.147560e-03,8.711770e-04,7.372474e-04,6.081937e-04,4.826448e-04,3.537937e-04,2.390395e-04,1.939472e-04,1.630473e-04,1.276292e-04,9.371911e-05,9.371911e-05];

const pop32x = [9.157955e-05,1.228611e-04,1.646393e-04,2.178680e-04,2.824279e-04,3.643750e-04,4.610986e-04,5.971983e-04,7.731033e-04,9.784002e-04,1.228819e-03,1.526477e-03,1.910629e-03,2.383524e-03,2.961833e-03,3.627829e-03,4.400116e-03,5.335784e-03,6.477959e-03,7.834784e-03,9.397937e-03,1.115095e-02,1.321140e-02,1.565896e-02,1.844841e-02,2.162804e-02,2.518167e-02,2.917022e-02,3.370022e-02,3.879573e-02,4.445690e-02,5.072892e-02,5.755398e-02,6.502591e-02,7.321177e-02,8.215358e-02,9.169969e-02,1.020341e-01,1.128786e-01,1.244548e-01,1.367191e-01,1.495320e-01,1.629561e-01,1.769504e-01,1.911957e-01,2.056707e-01,2.204244e-01,2.354560e-01,2.506282e-01,2.656687e-01,2.802975e-01,2.944952e-01,3.084002e-01,3.219038e-01,3.346854e-01,3.463770e-01,3.570553e-01,3.666397e-01,3.752281e-01,3.825732e-01,3.885933e-01,3.930099e-01,3.959838e-01,3.974781e-01,3.974910e-01,3.959047e-01,3.930560e-01,3.885337e-01,3.826339e-01,3.752287e-01,3.666564e-01,3.570688e-01,3.463482e-01,3.346649e-01,3.219127e-01,3.083443e-01,2.944862e-01,2.801934e-01,2.655960e-01,2.506602e-01,2.354719e-01,2.203714e-01,2.056244e-01,1.911906e-01,1.769321e-01,1.629989e-01,1.495549e-01,1.366781e-01,1.244824e-01,1.128593e-01,1.019536e-01,9.172483e-02,8.207774e-02,7.323801e-02,6.505485e-02,5.756207e-02,5.067689e-02,4.448271e-02,3.879800e-02,3.369653e-02,2.917349e-02,2.518469e-02,2.161650e-02,1.844626e-02,1.566009e-02,1.321347e-02,1.115670e-02,9.396583e-03,7.833228e-03,6.487419e-03,5.333328e-03,4.392766e-03,3.623043e-03,2.954889e-03,2.383911e-03,1.905139e-03,1.525130e-03,1.230828e-03,9.852073e-04,7.762922e-04,6.009952e-04,4.620463e-04,3.633200e-04,2.830597e-04,2.200913e-04,1.652294e-04,1.216868e-04,9.118018e-05,9.118018e-05];

const sum_diff    = array_op(ref,sum,   sub);
const pop_diff    = array_op(ref,pop,   sub);
const pop32_diff  = array_op(ref,pop32, sub);
const pop32x_diff = array_op(ref,pop32x,sub);

const ref_data      = {y: ref,      x0:-4, dx: 8./128., mode: 'lines', name: 'reference'};
const sum_data      = {y: sum,      x0:-4, dx: 8./128., mode: 'lines', name: 'sum'};
const pop_data      = {y: pop,      x0:-4, dx: 8./128., mode: 'lines', name: 'pop'};
const pop32_data    = {y: pop32,    x0:-4, dx: 8./128., mode: 'lines', name: 'pop32'};
const pop32x_data   = {y: pop32x,   x0:-4, dx: 8./128., mode: 'lines', name: 'pop32x'};

const sum_diff_data = {y: sum_diff, x0:-4, dx: 8./128., mode: 'lines', name: 'sum'};
const pop_diff_data = {y: pop_diff, x0:-4, dx: 8./128., mode: 'lines', name: 'pop'};
const pop32_diff_data  = {y: pop32_diff,  x0:-4, dx: 8./128., mode: 'lines', name: 'pop32'};
const pop32x_diff_data = {y: pop32x_diff, x0:-4, dx: 8./128., mode: 'lines', name: 'pop32x'};

const options = {displaylogo: false};

const elayout = {
  title:  'absolute error',
  yaxis:  {showline:false, hoverformat: 'g', exponentformat: 'power' },
  xaxis:  {range:[-4.,4.], nticks:9, zeroline:false },
  height: 400,
  width:  630,
};

Plotly.newPlot('sum', [ref_data,sum_data,pop_data,pop32_data,pop32x_data], elayout, options);
Plotly.newPlot('err', [sum_diff_data,pop_diff_data,pop32_diff_data,pop32x_diff_data], elayout, options);

</script>
