# Generating-Random-Numbers
Generation of random numbers via acceptance&rejection and Box-Müller methods.

Here we can find two different programes,one is focused on acceptance&rejection method, and the other one on Box-Müller method.

## Acceptance&Rejection
Generating random numbers with this method distributed according to the fun(x) distribution

$$p(x) = \frac{12}{\pi (2\pi^2-3)} x^2 \sin ^2(x)$$

And according to an exponencial distribution

$$p(x)=\lambda e^{-\lambda x}$$

We estimate the mean, variance, and standard deviation of the variable $x$ in order to discuss its behavior.

## Box-Müller.
Applied case of generating random numbers via Box-Müller.

### Electron inside a double-well potential (1 dimension)
The probability density of finding the electron at a position x is given by

$$p(x)= \frac{5}{324L}(x/L)^2 (9-(x/L)^2)$$

With $L= 3$ nm and $X\in[-3L,3L]$.

Here we generate 40000 possible values of the elctron position and we computate the probability of the electron to be in a region using Simpson method (numericla integration).

### Rb atom confined in 1D region

The probability density is given by:

$$g(x) = \sqrt{e^{-x^2/(2\sigma^2)}}{\sqrt{2\pi \sigma^2}}$$

We study the probability density and we estimate the mean, variance and standard deviation of the variable $x$. 
