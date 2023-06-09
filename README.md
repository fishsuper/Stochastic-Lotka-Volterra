# Stochastic-Lotka-Volterra
A type of numerical methods for a class of stochastic Poisson systems with invariant energy. The proposed numerical methods preserve both the energy and the Casimir functions of the systems.
The stochastic Lotka-Volterra systems are written as: 
$$
d \boldsymbol{y}(t)=\boldsymbol{B}(\boldsymbol{y}(t)) \nabla H(\boldsymbol{y}(t))(d t+c \circ d W(t)),
$$
where
$$
\begin{aligned}
\boldsymbol{B}(\boldsymbol{y}) & =\left(b_{i j}^0 y^i y^j\right)=\operatorname{diag}\left(y^1, \ldots, y^m\right) \boldsymbol{B}_0 \operatorname{diag}\left(y^1, \ldots, y^m\right), \\
H(\boldsymbol{y}) & =\sum_{i=1}^m \beta_i y^i-p_i \ln y^i
\end{aligned}
$$
The numerical experiments are prefromed by Matlab. OurMethods.m presents the efficiency of our algorithms, StrongOrder.m tests the strong order of these methods.
Find more details on http://www.math.ualberta.ca/ijnam/Volume-19-2022/No-2-22/2022-02-03.pdf
