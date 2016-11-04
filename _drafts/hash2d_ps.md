---
layout: post
title:  2D hash power spectrum
tags : [random, low-discrepency, sampling]
description : visualizing
---

The following is a math block:

$$a\times b$$

But next comes a paragraph with an inline math statement:

Foo\$$ 5 + 5 $$

{% highlight c++ %}
  static inline uint32_t pack32(uint32_t x, uint32_t y)
  {
    return (y<<16) ^ x;
  }
{% endhighlight %}


<iframe width="420" height="315" src="http://www.shadertoy.com/embed/MddGDB" frameborder="0" allowfullscreen></iframe>


{% highlight c %}
uint32_t hash(uint32_t x, uint32_t y)
{
  x *= 0x3504f333;   //
  y *= 0xf1bbcdcb;   //
  x ^= y;
  x *= 741103597;    // 
  return x;
}
{% endhighlight %}



$$
\begin{align*}
  & \phi(x,y) = \phi \left(\sum_{i=1}^n x_ie_i, \sum_{j=1}^n y_je_j \right)
  = \sum_{i=1}^n \sum_{j=1}^n x_i y_j \phi(e_i, e_j) = \\
  & (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right)
\end{align*}
$$