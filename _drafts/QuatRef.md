---
layout:       post
title:        Quaternion reference definations and notation
tagline:      
categories:   [quaternions]
tags:         [interpolation, quantization]
description:  foo bar baz
---


------

Sets
------

### Set builder

{% highlight c %}
SetDefine: SetName '=' '{' SetVariable ':' SetContraints '}'

SetConstraints: SetContraint
                SetContraint ';' SetContraints


SetName: SetName '{' variable ':' contraints '}'

{% endhighlight %}


$$ \begin{equation*}
\mathbb{X}=\left\{ \text{u} : \text{abc} \right\} 
\end{equation*} $$


### Set builder

$$ \begin{equation} \label{set:E3}
\mathbb{E}^{3}=\left\{ \mathbf{u}\in\mathbb{R}^{3}:\left|\mathbf{u}\right|=1\right\} 
\end{equation} $$

$$ \begin{equation} \label{set:S3}
\mathbb{S}^{3}=\left\{ \mathbf{u}\in\mathbb{R}^{4}:\left|\mathbf{u}\right|=1\right\} 
\end{equation} $$

$$ \begin{equation} \label{set:S2}
\mathbb{S}^{2}=\left\{ \mathbf{u}\in\mathbb{R}^{3}:\left|\mathbf{u}\right|=1\right\} 
\end{equation} $$

$ \newcommand{\hcomplex}[1]{\mathbb{C}_{#1}} $
$$ \begin{equation} \label{Cu}
\hcomplex{u}=\left\{ \left(a + b~\mathbf{u}\right) \in\mathbb{H}:a,b~\mathbb{\in R}~; \mathbf{u}\in\mathbb{S}^{2}\right\} 
\end{equation} $$


\\
Returns the current position in the sequence.

<br>

------

SLERP <small>the dreaded circular arc</small>
------

\\
Returns the current position in the sequence.

<br>

------

Basic analysis tools <small></small>
------

\\
Returns the current position in the sequence.

<br>

------

Reference formulation of slerp
------

\\
Returns the current position in the sequence.

<br>

------

xxx
------

#### LERP

\\
Returns the current position in the sequence.


<br>

------

References and Footnotes
------
