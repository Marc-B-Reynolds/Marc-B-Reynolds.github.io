---
layout:       post
title:        'Popcount walks: next, previous, toward and nearest'
categories:   [math]
tags:         [integer,bithack]
description:  'Informal (and by-example) derivation of walking adjacent integers with the same population count as the input.'
plotly:       false
---

\\
Harold Aptroot posted a method to compute the nearest integer with the same [population count](https://en.wikipedia.org/wiki/Hamming_weight) (`popcount`) as the input (see below). This got me thinking about the various enumeration methods so I'll toss out some comments. I'll show 32-bit example code (64-bit follows directly) and use 16-bit for visual examples.

Implementations in 32 & 64 bit with visual examples (like below) can be found [HERE](https://gcc.godbolt.org/z/8e46EaWYY).

<br>

------

Next popcount <small></small>
------

\\
Visting the next greater integer with the same popcount is the commonly seen method. Let's walk through to a method (found in *"Hacker's Delight"*, second edition in exercises. Warren credits David de Kloet) which has the interesting properly that if the input is the maximum integer (of the given popcount) then all ones is returned.

First let's look at what happens if we were to just increment until we find the result:

{% highlight c %}
  p = popcount(x);
  a = 1;
  
  while(popcount(x+a) != p) a++;
  
  return x+a;
{% endhighlight %}

\\
with the following example input and after a "few" iterations (.=0 x=doesn't matter)

$$ \begin{array}{lll}
\texttt{x}   & = & \texttt{xxxxxx.1.111....} \\
\texttt{a}   & = & \texttt{............1111} \\
\texttt{x+a} & = & \texttt{xxxxxx.1.1111111}
\end{array} $$


\\
all values of $a$ up to and including the above (as many 1s as there are trailing zeroes in the input) result in $(x+a)$ having a greater popcount than $x$ so we speed up the increment method by initializing $a$ to the value of the lowest set bit in $x$ which can be expressed as `(x & -x)`.

Our example then looks like:

$$ \begin{array}{lll}
\texttt{x} & &                 &=& \texttt{xxxxxx.1.111.... : } \text{input} \\
\texttt{a} &=& \texttt{x & -x} &=& \texttt{...........1.... : } \text{isolate lowest set bit} \\
\texttt{b} &=& \texttt{x + a}  &=& \texttt{xxxxxx.11....... : } \text{closer!} \\
\end{array} $$

\\
which generally has a smaller popcount than $x$. Specifically if we call $L$ the length of the lowest run of 1s (here 3) and target popcount $p$ then $\text{popcount}\left(b\right) = p-(L-1)$ because the $L$ 1s of the run flip to zero and the first zero after flips to 1. So to get to our desired result we'd need to keep adding until the bottom $(L-1)$ bits are set.

Now let's note that we can isolate the lowest run of 1s in a integer using our temporary value $b$ as `(x & ~b)`. We simply need to right shift this value by one more than it's *"trailing zero count"* to have the needed number of 1s in place:

$$ \begin{array}{lll}
\texttt{x} &   &                      & = & \texttt{xxxxxx.1.111.... : } \text{input}\\
\texttt{a} & = & \texttt{x & -x}      & = & \texttt{...........1.... : } \text{isolate lowest set bit}    \\
\texttt{b} & = & \texttt{x + a}       & = & \texttt{xxxxxx.11....... : } \text{cascade the run}           \\
\texttt{c} & = & \texttt{x & ~b}      & = & \texttt{.........111.... : } \text{isolate lowest run of 1s}  \\
\texttt{c} & = & \texttt{c >> ctz(c)} & = & \texttt{.............111 : } \text{  * move to bottom}        \\
\texttt{c} & = & \texttt{c >> 1}      & = & \texttt{..............11 : } \text{  * lose one more bit}     \\
\texttt{r} & = & \texttt{b ^ c}       & = & \texttt{xxxxxx.11.....11 : } \text{result (add/or/xor all work)} 
\end{array} $$


\\
Note that this shifting twice thing is to sidestep *undefined behavior* (UB) and we'd have to add one anyway to flip it to a single shift which is a zero sum gain.  A direct translation into code then looks like:

{% highlight c %}
static inline uint32_t pop_next_32_verbose(uint32_t x)
{
  uint32_t a = x & -x;   // 1) isolate lowest set bit (lsb) of 'x'
  uint32_t b = x +  a;   // 2) add lsb to 'x'
  uint32_t c = x & ~b;   // 3) isolate lowest bit string (all 1's) of 'x'
  uint32_t d;
  
  // if the length of the lowest bit string of 'x' (popcount(c)) is 'L'
  // then popcount(b) = popcount(x)-(L-1). We need (L-1) set bits at the
  // bottom to "add" to 'b'.
  d = c >> (ctz_32(c) & 0x1f); // 4) shift c to bottom. avoid UB for x=0
  d = d >> 1;                  // 5)   and 1 more to have L-1

  return b^d;                  // 6) result
}
{% endhighlight %}


\\
This currently returns a *funky* result when the input is the maximum value of a given popcount. The max value is when the lowest bit string ($b$ above) includes the sign bit and our pair of shifts then result in the smallest number with a popcount of $\left(p-1\right)$ mod the bitwidth. 

If instead we use signed (arithmetic) shifts then we get the above mentioned "feature" of max value input returning all 1s. 

{% highlight c %}

#if defined(__x86_64__)
static inline uint32_t ctz_32(uint32_t x) { return (uint32_t)_tzcnt_u32(x); }
#else
// this works for Intel as well and compilers lower to one opcode (YMMV)
static inline uint32_t ctz_32(uint32_t x) { return (x!=0) ? (uint32_t)__builtin_ctz(x)  : 32; }
#endif

// signed shift of unsigned + no UB on shift
static inline uint32_t shr_32(uint32_t x, uint32_t n)
{
  return (uint32_t)((int32_t)x >> (n & 0x1f));
}

// next larger integer with same popcount as 'x'
// * input of zero returns zero
// * max value of a given popcount returns all 1s (-1)
static inline uint32_t pop_next_32(uint32_t x)
{
  uint32_t t = x + (x & -x);
  x = x & ~t;
  x = shr_32(x, ctz_32(x));
  x = shr_32(x, 1);
  return t^x;
}
{% endhighlight %}

\\
Example lowering into Intel and ARM (which can probably be improved by an ARM smart person):

<div class="container">
  <div class="row">
    <div class="col-sm">

{% highlight nasm %}
; pop_next_32 (intel)
blsi    ecx, edi
add     ecx, edi
andn    eax, ecx, edi
tzcnt   edx, eax
sarx    eax, eax, edx
sar     eax
xor     eax, ecx
ret
{% endhighlight %}

</div>
<div class="col-sm" markdown="1">

{% highlight nasm %}
; pop_next_32 (ARM)
neg     w8, w0
and     w8, w8, w0
add     w8, w8, w0
bic     w9, w0, w8
rbit    w10, w9       ; lack of trailing zero count xforms to
clz     w10, w10      ;   bit-reverse/clz. input is single bit
asr     w9, w9, w10   ;   or zero. exercise for reader! 
eor     w0, w8, w9, asr #1
ret
{% endhighlight %}
</div>
</div>
</div>


<br>

------

Previous popcount  <small></small>
------

\\
I couldn't find any efficient methods to walk to the next smaller. The only constant time version I found was a mod of R.W. Gosper HAKMEM method which looks like this:

{% highlight c %}
uint32_t pop_prev_DO_NOT_USE_32(uint32_t x)
{
  uint32_t a = x + 1;
  uint32_t b = (a ^ x) + 1;
  uint32_t c = (a & x);
  uint32_t d = c & -c;

  // had to add this since two inputs produce b=0 and BOOM!
  if (b != 0) 
    return c - d/b;

  // bottom (bit_width-1) bits are set
  return (uint32_t)((int32_t)x >> 31);
}
{% endhighlight %}

\\
Needing to perform an integer division is very unfortunate. However since we already have a version that walks to the next we can simply apply a linear transform to map/unmap the input/result.  The ones complement `(~x)` reverses the notion of unsigned integer ordering so the follow function walks to the previous:

{% highlight c %}
static inline uint32_t pop_prev_32(uint32_t x) { return ~pop_next_32(~x); }
{% endhighlight %}

\\
Which isn't too bad at +2 basic logical ops. We could *"carry through the derivation"*: pull the complements inside, apply identities and win!  But's lets just re-run through the *pop_next* example:


$$ \begin{array}{lll}
\texttt{x} &   &                      & = & \texttt{xxxxxx1.1...1111 : } \text{input: ~x of 'next' example}\\
\texttt{a} & = & \texttt{~x & (x+1)}  & = & \texttt{...........1.... : } \text{isolate lowest clear bit}  \\
\texttt{b} & = & \texttt{ x -  a}     & = & \texttt{xxxxxx.11....... : } \text{cascade the run}           \\
\texttt{c} & = & \texttt{~x &  b}     & = & \texttt{.........111.... : } \text{isolate lowest run of 0s}  \\
\texttt{c} & = & \texttt{c >> ctz(c)} & = & \texttt{.............111 : } \text{  * move to bottom}        \\
\texttt{c} & = & \texttt{c >> 1}      & = & \texttt{..............11 : } \text{  * lose one more bit}     \\
\texttt{r} & = & \texttt{b ^ c}       & = & \texttt{xxxxxx1..11111.. : } \text{zero out bottom bits} 
\end{array} $$

\\
List of modifications needed:
* `x` is changed to `~x` of *pop_next* case because we're interested in odd inputs here and even there (more below)
* `a` goes from lowest set `(x & -x)` to lowest clear `~x & (x+1)`
* `b` goes from `x+a` to `x-a`. Both flip the run and the next higher bit
* `c` remains the same
* `r` remains the same but requires using XOR to flip the needed bits


{% highlight c %}
// next smaller integer with same popcount as 'x'
// * if the input is zero or min value of a given
//   popcount then returns zero.
uint32_t pop_prev_32(uint32_t x)
{
  uint32_t a = ~x & (x+1);  // 1) lowest clear bit of 'x'
  uint32_t b =  x - a;      // 2) cascade run of 0s to 1s
  uint32_t c = ~x & b;      // 3) isolate lowest run of zeroes of 'x'

  // 4) move isloated bits to bottom and remove one
  uint32_t d = shr_32(c, ctz_32(c));
  d = shr_32(d, 1);

  return b^d;               // 5) flip the bits
}
{% endhighlight %}


\\
Example lowering into Intel for transformed next vs. reworked:

<div class="container">
  <div class="row">
    <div class="col-sm">

{% highlight nasm %}
; pop_prev_32 (intel, via transform)
lea     eax, [rdi + 1] ; compiler transforms: (~x) & -(~x)
andn    eax, edi, eax  ;   to ~x & (x+1)
mov     ecx, edi
sub     ecx, eax
andn    ecx, edi, ecx
mov     edx, edi
not     edx
add     edx, eax
tzcnt   eax, ecx
sarx    eax, ecx, eax
sar     eax
or      eax, edx
not     eax
ret
{% endhighlight %}

</div>
<div class="col-sm" markdown="1">

{% highlight nasm %}
; pop_prev_32 (intel, reworked)
lea     eax, [rdi + 1]
andn    eax, edi, eax
mov     ecx, edi
sub     ecx, eax
andn    eax, edi, ecx
tzcnt   edx, eax
sarx    eax, eax, edx
sar     eax
xor     eax, ecx
ret
{% endhighlight %}
</div>
</div>
</div>



<br>

------

Popcount toward some target direction
------

\\
Given the transform version of `pop_prev` above it's straightword to create runtime selection of the walk direction.

{% highlight c %}
// m =  0 : next
// m = ~0 : prev
// otherwise: fun?
static inline uint32_t pop_walk_directed_32(uint32_t x, uint32_t m)
{
  uint32_t t;

  x ^= m;
  
  // pop_next core
  t  = x + (x & -x);
  x  = x & ~t;
  x  = shr_32(x, ctz_32(x));
  x  = shr_32(x, 1);
  t ^= x;

  return m^t;
}

// closest integer to 'x' with the same popcount in the direction of 'y'
static inline uint32_t pop_toward_32(uint32_t x, uint32_t y)
{
  uint32_t dir = (y>x) ? 0 : UINT32_C(~0);
  pop_walk_directed_32(x, dir);
}
{% endhighlight %}

<br>

------

Nearest popcount
------

\\
Now back to Harold's post:

<p align="center" vertical-align="middle">
<iframe src="https://mastodon.gamedev.place/@harold/111281691234096070/embed" class="mastodon-embed" style="max-width: 100%; border: 0" width="400" allowfullscreen="allowfullscreen"></iframe><script src="https://mastodon.gamedev.place/embed.js" async="async"></script>
</p>

\\
Let's write out a 32-bit version in C:

{% highlight c %}
static inline uint32_t pop_nearest_32(uint32_t x)
{
  uint32_t a = -x & (x + 1);  // lowest change bit (see below)
//uint32_t b = a ^ (a >> 1);  // not a Gray code so I'd want to write it
  uint32_t b = a | (a >> 1);  // like this instead (see below)
  return x ^ b;
}
{% endhighlight %}


\\
The first temp value $a$ gives the "first bit that differs from bit-0". So for even $x$ it isolates the lowest set bit and for odd the lowest clear bit. If $x$ is either all zeros or all ones then $a$ is zero. Otherwise $a$ is always a single set bit with the lowest clear.

The second temp value $b$ looks like it's computing a [Gray code](https://en.wikipedia.org/wiki/Gray_code) but given the possible values of $a$ it's simply creating a number with 2 adjacent bits set (so the above code mod).

If in the by-example derivation of `pop_next` we'd used an odd integer then we would have came up with the `pop_nearest` above except for the computation of `a` would have been 'isolate the lowest clear bit' `(~x & (x+1))` and if we'd carried through with an even integer for `pop_prev` then `a` would have been 'isolate the lowest set bit' `(x & -x)`. So this neat little thing just marries together the simple cases of next and prev. Nice!


Let's look at way happens with 5 trailing zeroes (x0) and 5 trailing ones (x1):

    xxxxxxxxxx100000 : x0 = even x
    xxxxxxxxxx011111 : x1 = odd  x
    0000000000100000 : both yield this as 'a'
	0000000000110000 : both yield this as 'b'
    
    xxxxxxxxxx010000 : x0 ^ b
    xxxxxxxxxx101111 : x1 ^ b

			  
	
\\
Notice that both decrease the length of the low run by 1 so we can percolate down:

    xxxxxxxxxx000010 : x0 = even x
    xxxxxxxxxx111101 : x1 = odd  x
    0000000000000010 : both yield this as 'a'
    0000000000000011 : both yield this as 'b'
    
    xxxxxxxxxx000001 : x0 ^ b
    xxxxxxxxxx111110 : x1 ^ b
              -----  : each leave behind the original bit-0
                       and lowest 2 bits toggle each step from now on

