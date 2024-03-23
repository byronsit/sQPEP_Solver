# sQPEP_Solver
The corresponding paper is currently under review. Detailed information and the full paper will be shared once the review process is complete.



## Proof the Lemma

Given the description of the theorem, our objective is to demonstrate that for a polynomial function $f(x_1, x_2, \ldots, x_m) = C \prod_{i=1}^m x_i^{p_i}$, the product of its $n^{th}$ order partial derivatives with respect to all variables and $vec(S^{\otimes n})$ is equal to $n!$ times the function itself. Here, $C$ is a constant, $p_i$ are non-negative integers, and $\sum_{i=1}^m p_i = n$.


To establish this, it is essential to consider all possible $n^{th}$ order partial derivatives. For each specific set of $(k_1, k_2, \ldots, k_m)$ that satisfies $\sum_{i=1}^m k_i = n$, the $n^{th}$ order mixed partial derivative of $f$ (provided that all $k_i\! \leq\! p_i$) is given by:


$$
\frac{\partial^n f}{\partial x_1^{k_1} \partial x_2^{k_2} \ldots \partial x_m^{k_m}}\! =\! C \cdot \prod_{i=1}^m \frac{p_i!}{(p_i - k_i)!} \cdot x_i^{p_i - k_i},    
$$




This is due to the fact that for each variable $x_i$, if we take the $k_i^{th}$ partial derivative, the remaining power is $p_i - k_i$, accompanied by a coefficient $\frac{p_i!}{(p_i - k_i)!}$. Should $k_i\! >\! p_i$, the result is zero, in accordance with the definition of partial derivatives.

Next, we turn our attention to $vec(S^{\otimes n})$. This vector comprises all elements of $S^{\otimes n}$, where $S$ is the set formed by the variables of $f$. Consequently, each element within $vec(S^{\otimes n})$ corresponds to the coefficient of a specific $n^{th}$ order partial derivative, represented by the multinomial coefficient $\binom{n}{k_1, k_2, \ldots, k_m}$. This coefficient is the number of ways to distribute $k_i$ derivatives to each $x_i$ from the $n$ derivatives. Here, the multinomial coefficient $\binom{n}{k_1, k_2, \ldots, k_m}$ is defined as $\frac{n!}{k_1! k_2! \ldots k_m!}$.


Now, we can compute the $\mathbf{D}^{(n)} \cdot vec(S^{\otimes n})$:
$$
\sum_{\sum_{i=1}^m k_i = n} \left( C \cdot \prod_{i=1}^m \frac{p_i!}{(p_i - k_i)!} \cdot x_i^{p_i - k_i} \right) \cdot \left( \frac{n!}{k_1! k_2! \ldots k_m!} \right).
$$


Because we have a complete permutation of $n!$, it will be multiplied by each $\frac{p_i!}{(p_i - k_i)!}$. However, to compute the dot product, we must also consider the corresponding $x_i^{k_i}$ terms in $vec(S^{\otimes n})$. It is noted that the product of these terms yields the powers of the original polynomial $f$. This is due to each term $x_i^{p_i - k_i} \cdot x_i^{k_i}$ ultimately resulting in $x_i^{p_i}$. Consequently, we obtain:
$$
\mathbf{D}^{(n)} \cdot vec(S^{\otimes n}) = C \cdot n! \cdot \prod_{i=1}^m x_i^{p_i} = n! \cdot f.    eee
$$
