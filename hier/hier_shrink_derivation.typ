#set page(paper: "us-letter", margin: 1in)
#set text(size: 11pt)
#set par(justify: true)
#set math.equation(numbering: "(1)")

#align(center)[
  #text(size: 16pt)[Optimal shrinkage in a simple hierarchical model]

  BIOS/BCB 784, companion to `hier/hier_simulation.R`
]

= Setup

For genes $i = 1, dots, m$ and replicates $j = 1, dots, n$:

$ theta_i tilde N(mu_theta, sigma_theta^2), quad
  x_(i j) | theta_i tilde N(theta_i, sigma_x^2). $

Let $macron(x)_i = 1/n sum_j x_(i j)$. We consider the shrinkage estimator

$ hat(theta)_i = (1 - B) macron(x)_i + B mu_theta, quad B in [0, 1], $

where $mu_theta$, $sigma_theta^2$ and $sigma_x^2$ are treated as known.

= Question

Work out the value of $B$ that minimizes the mean squared error
$E[(hat(theta)_i - theta_i)^2]$, averaging over both $theta_i$ and the data.

= Answer

Given $theta_i$, the average of $n$ independent draws has

$ macron(x)_i | theta_i tilde N(theta_i, sigma^2), quad "where" quad
  sigma^2 = sigma_x^2 / n. $

Adding and subtracting $theta_i$,

$ hat(theta)_i - theta_i
  = (1 - B)(macron(x)_i - theta_i) + B(mu_theta - theta_i). $

Let's call the two pieces $e_i = macron(x)_i - theta_i$ (the measurement
error) and $d_i = mu_theta - theta_i$ (how far the gene's true value is from
the population mean). Squaring the expression above gives

$ (hat(theta)_i - theta_i)^2
  = (1 - B)^2 e_i^2 + B^2 d_i^2 + 2 B (1 - B) e_i d_i, $

and we need the expected value of each of the three terms.

- The measurement error $e_i$ is $N(0, sigma^2)$ no matter what value
  $theta_i$ takes. So $E[e_i] = 0$ and $E[e_i^2] = "Var"(e_i) = sigma^2$.
- The gene's deviation $d_i$ comes from the prior, $theta_i tilde
  N(mu_theta, sigma_theta^2)$. So $E[d_i] = 0$ and $E[d_i^2] = sigma_theta^2$.
- The measurement error does not depend on which $theta_i$ the gene has,
  so $e_i$ and $d_i$ are independent. Then
  $E[e_i d_i] = E[e_i] E[d_i] = 0$, and the cross term drops out.

Putting these together:

$ "MSE"(B) = E[(hat(theta)_i - theta_i)^2]
  = (1 - B)^2 sigma^2 + B^2 sigma_theta^2. $ <mse>

The first term is the variance of the estimator and the second term is
its average squared bias. Increasing $B$ lowers the variance but raises the
bias, so the best choice is somewhere in between. This is the same
bias-variance trade-off that appears in ridge regression, smoothing, and
many other settings where we accept a little bias in exchange for a large
reduction in variance. Taking the derivative,

$ (d "MSE") / (d B) = -2 (1 - B) sigma^2 + 2 B sigma_theta^2 = 0, $

and solving for $B$,

$ B^* = sigma^2 / (sigma^2 + sigma_theta^2). $ <bopt>

To interpret $B^*$, look at the variance of the observed averages. We can
write $macron(x)_i = theta_i + e_i$, the true value plus measurement error.
These two pieces are independent, so their variances add:

$ "Var"(macron(x)_i) = "Var"(theta_i) + "Var"(e_i) = sigma_theta^2 + sigma^2. $

The denominator of $B^*$ is exactly this total variance, and the numerator
is the part of it due to measurement error:

$ B^* = sigma^2 / "Var"(macron(x)_i) = "noise variance" / "total variance". $

So $B^*$ is the fraction of the spread in the $macron(x)_i$ that is noise.
If most of the differences between genes' averages are just measurement
error, the $macron(x)_i$ are unreliable and we should shrink them strongly
toward $mu_theta$. If most of the differences reflect real differences in
$theta_i$, we should mostly trust each gene's own average.

With the parameters in the script ($n = 5$, $sigma_x^2 = 4$,
$sigma_theta^2 = 1.5$), we have $sigma^2 = 4 / 5 = 0.8$ and

$ B^* = 0.8 / (0.8 + 1.5) = 0.348, $

which is close to the minimum found by grid search in the simulation. Plugging
$B^*$ back into @mse gives the minimum MSE,

$ "MSE"(B^*) = (sigma_theta^2 sigma^2) / (sigma_theta^2 + sigma^2)
  = 1.2 / 2.3 = 0.522. $

Compare this to the alternatives:

- No shrinkage ($B = 0$): $"MSE"(0) = sigma^2 = 0.8$.
- Complete shrinkage ($B = 1$): $"MSE"(1) = sigma_theta^2 = 1.5$.
- The variance of the data, $"Var"(macron(x)_i) = sigma_theta^2 + sigma^2 = 2.3$.

The optimal MSE is smaller than both $sigma^2$ and $sigma_theta^2$: it is
one half of their harmonic mean, so its inverse is the sum of the two
precisions, $1 \/ "MSE"(B^*) = 1 \/ sigma^2 + 1 \/ sigma_theta^2$.
Equivalently, $"MSE"(B^*) = sigma_theta^2 sigma^2 \/ "Var"(macron(x)_i)$,
the product of the signal and noise variances divided by the total
variance of the $macron(x)_i$.
