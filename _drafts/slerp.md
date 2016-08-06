---
layout:       post
title:        On certain results of quaternionic interpolation
tagline:      (or let's walk a 2D path)
categories:   [quaternions]
tags:         [interpolation]
description:  foo
---

This is a quick note to walk through slerp and how the same notions can be applied to other basic forms of quaternion interpolation.

* foo
* bar
* baz


------

Real power
------

\\
We start with some unit quaternion:

$$
\begin{eqnarray}
  Q & = & w + \left(x,~y,~z\right) \label{qxform} \\ 
    & = & \cos\left(\theta\right)+\sin\left(\theta\right)~\mathbf{u} \label{qpform}
\end{eqnarray}
$$

where we explictly have the four real values of $ \eqref{qxform} $ which implictly represents a rotation of $2\theta$ about $\mathbf{u}$ as per $ \eqref{qpform} $.  


\\
To linearly parameterize the angle between zero and $\theta$ :

$$ \begin{equation} \label{eq:ppower}
Q^{t} = \cos\left(t~\theta\right)+\sin\left(t~\theta\right)~\mathbf{u}
\end{equation} $$

\\
To linearly parameterize the angle between zero and $\theta$ :

<br>

------

SLERP <small>the dreaded circular arc</small>
------

\\



$$ \begin{equation} \label{eq:slerp}
s\left(t\right) = \left(BA^{*}\right)^{t}A = A\left(A^{*}B\right)^{t}
\end{equation} $$

$$
\begin{align}
  s\left(t\right) = \left(BA^{*}\right)^{t}A  \\
                  = A\left(A^{*}B\right)^{t}
		  \label{eq:slerpl}
\end{align}
$$


------

Basic analysis tools <small></small>
------

\\
The behavior of slerp is defined by $ \eqref{eq:ppower} $ which operates in the complex plane spanning $Q$ and the line of reals.

$$ \begin{equation} \label{eq:fpoint}
p\left(t\right) = \left(p_x\!\left(t\right),~p_y\!\left(t\right)\right)
\end{equation} $$

$$ \begin{equation} \label{eq:fangle}
a\left(t\right) =\arctan \left(\frac{p_y\!\left(t\right)}{p_x\!\left(t\right)}\right)
\end{equation} $$

\\
To place our 2D coordinate back into the original space we need to directions in that space to  $A~(f_x\left(t\right)$


$$
\left(A\cdot A\right)B - \left(A\cdot B\right)A
$$

and since we are limiting ourselves to unit magnitudes this reduces to:

$$
B - \left(A\cdot B\right)A
$$


$$
\frac{B - \left(A\cdot B\right)A}{\sqrt{1 -\left(A\cdot B\right)^{2}}}
$$


$$
$$

<br>

------

Reference formulation of slerp
------

\\
The algebraic form of slerp $ \eqref{eq:slerp} $

$$ \begin{equation} \label{eq:davis}
s\left(t\right) = \frac{\sin\left(\left(1-t\right)\theta\right)}{\sin\left(\theta\right)}A+\frac{\sin\left(t\theta\right)}{\sin\left(\theta\right)}B
\end{equation} $$

------

xxx
------

#### LERP

\\
Returns the current position in the sequence.

$$
\frac{\left(A\cdot A\right)B - \left(A\cdot B\right)A}{\sqrt{\left(A\cdot A\right)\left(B\cdot B\right)-\left(A\cdot B\right)^{2}}}
$$



$$ 
 p_x\left(t\right) = 1+t\left(\cos\left(\theta\right)-1\right)
$$

$$
p_y\left(t\right) = t~\sin\left(\theta\right)
$$



------

xx
------


------

References and Footnotes
------
[^rotation]: "Animating Rotation with Quaternion Curves", Ken Shoemake, 1985. ([PDF](http://run.usc.edu/cs520-s12/assign2/p245-shoemake.pdf))
[^shoemake85]: "Animating Rotation with Quaternion Curves", Ken Shoemake, 1985. ([PDF](http://run.usc.edu/cs520-s12/assign2/p245-shoemake.pdf))