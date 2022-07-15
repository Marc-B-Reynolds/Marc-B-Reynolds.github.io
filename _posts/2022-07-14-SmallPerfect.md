---
layout: post
title:  "A fast perfect shuffle for n&#8804;64"
categories: [distribution]
tags : [random]
description : "Describes a perfect shuffle using a bitset and the original pencil-and-paper Fisher-Yates method."
plotly: false
---

The problem is to perform a perfect shuffle when $n \le 64$ and we're going to attack it with the standard method: The [Fisher-Yates shuffle](https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle).  This is sufficient to shuffle a standard deck of card or to generate a $64\times64$ permutation matrix. 

The important assumption is that the target hardware has a *bit scatter* opcode (and not microcoded: AMD prior to Zen3 I'm looking at you) and is exposed in the target language (`PDEP` in Intel land).  Let's jump straight in:

{% highlight c %}

// remember: this is strawman/example code only
typedef struct { uint64_t cards; uint8_t pop; } foo_t;

// set up for a permutation of 'n' with n<=64
void foo_init(foo_t* data, uint32_t n)
{
  data->cards = UINT64_C(-1) >> (64-n);
  data->pop   = n;
}

{% endhighlight %}

\\
The first step of *Fisher-Yates* is to initialize an array to: $0$ through $n-1$. Instead of an explicit array we're initializing a *bit-set* (`cards`) to represent our (logically) unshuffled deck.  We're storing $n$ in `pop` to mentally note that this is redundant data which is maintained to be the [population count](https://en.wikipedia.org/wiki/Hamming_weight) of `cards`.

For the rest we're ignoring modern variants and are running with the original *pencil-and-paper* method (except for the counting from zero instead of one part):


{% highlight c %}

// wrap the Intel family specific bit scatter intrinsic
static inline uint64_t bit_scatter_64(uint64_t x, uint64_t m) 
{ 
  return _pdep_u64(x, m); 
} 

// wrap the Intel family trailing zero count. 
// WARNING: __builtin_ctzl is treated as UB with input of zero
// by at least clang (on intel targets) and intel compilers. So
// avoiding that pit.
static inline uint32_t ctz_64(uint64_t x) { return (uint32_t)_tzcnt_u64(x); }


// return 'x' with the n^th set bit cleared (counting from the LSB). So
// 0 clears lowest, 1=second lowest, etc. If n >= popcount(x) then
// 'x' is return unmodified.
static inline uint64_t bit_clear_nth_set_64(uint64_t x, uint32_t n)
{
  const uint64_t m = ~UINT64_C(0);     // all bits set
  const uint64_t o =  UINT64_C(1);     // 64-bit '1'
  const uint64_t t = m ^ (o<<n);       // all bits set except the nth

  return bit_scatter_64(t, x);
}

// randomly draw a remaining card. if 'empty' it returns 64
uint32_t foo_next(foo_t* gen)
{
  uint32_t i = rng_range(gen->pop);
  uint64_t c = bit_clear_nth_set_64(gen->cards, i);
  uint64_t r = gen->cards ^ c;
  
  gen->cards = c;
  gen->pop--;
  
  return ctz_64(r);
}
{% endhighlight %}

\\
So the first two helper functions are just intrinsic wrappers for the scatter op (`bist_scatter_64`) and *trailing zero count* (`ctz_64`). The next function (`bit_clear_nth_set_64`) is lynchpin of the method that performs the *strike* the n<sup>th</sup> card operation. The `m` constant is just all ones, `o` is just one and `t` is all ones except the n<sup>th</sup> bit is clear.  The result of the function is to apply the bit-scatter op. Notice the ordering of the parameters: the input `x` is being used as the *selector* and the `t` is value that is being scattered. This transforms our n<sup>th</sup> bit clear `t` into clearing the n<sup>th</sup> set bit of `x`.

Finally `foo_next` just snaps all the pieces together:

* `i` uniformly selects a remain card in our unshuffled deck
* `c` is the updated `cards` (selected one striked)
* `r` isolates which card was striken
*  update the data structure with new cards and decrease the popcount
* return the number of the drawn card which is the trailing zero count of its isolated position

Gluing the initialization code and looping over the next function to fill an array gives a standardish kind of thing:

{% highlight c %}
// fill array 'a' with a random permutation of {0..n-1}
void foo_array_gen(uint32_t n, uint32_t a[n])
{
  foo_t gen;
  foo_init(&gen, n);
  
  for(int i=0; i<n; i++)
    a[i] = foo_next(&gen);
}
{% endhighlight %}


\\
And that's it! "Obviously" this can be extended to supporting larger $n$. As examples for some $64 \cdot m$ max value with $m$ very small then a linear search to find the correct 64-bit window of the bitset to modifiy is pretty straight forward and one could push to $64 \cdot 2^p$ which can be performed in $p$ steps to find the window and update the hierarchical popcounts. Generally using a standard method instead would probably be the better choice however YMMV.

Closing note: For the presented toy code to be "correct" requires that `rng_range` returns a uniform result. There are recent solutions which are quite performant like the one by Daniel Lemire which is talked about in this [blog post](https://lemire.me/blog/2016/06/30/fast-random-shuffling) (paper linked at end).

On review the example code is rather ugly so here's generating a binary $64\times64$ permutation matrix. The real version doesn't use a global PRNG...but you get the idea. Otherwise this is the example code minus the zero-counting which isn't needed here (row result is the choosen bit position).

{% highlight c %}
void bmat_random_permutation_64(uint64_t m[static 64])
{
  // initialize the bitset to contain [0..64]
  uint64_t cards = UINT64_C(-1);
  uint32_t row   = 64;

  // for each row perform a pencil & paper Fisher-Yates
  do {
    uint32_t u = rng_range(row--);                // uniform 'u' on [0,row)
    uint64_t t = UINT64_C(-1) ^ (UINT64_C(1)<<u); // all set except bit 'u'
    uint64_t c = bit_scatter_64(t, cards);        // clear u^th element of bitset

    // the row is the isolated element we just cleared
    m[row] = cards ^ c;
    cards  = c;                                   // update the bitset
  } while(row);
}
{% endhighlight %}



