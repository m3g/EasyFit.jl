# 
# Least-square function and gradient
#
# f is the inner function and gf is the gradient of the
# function applied to each element of the data vector
#
function sq_residue(p,f,X,Y)
  s = 0.
  for i in 1:length(X)
    s += (Y[i] - f(X[i],p))^2
  end
  s
end
function sq_residue_grad!(p,f,∇f,X,Y,g)
  @. g = 0.
  for i in 1:length(X)
    dsdf = -2*(Y[i] - f(X[i],p))
    grad = ∇f(X[i],p)
    for i in 1:length(g)
      g[i] = g[i] + dsdf*grad[i]
    end
  end
end
