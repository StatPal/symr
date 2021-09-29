# MRF_log_likeli <- function(W, Psi_inv, beta1) {
#  # (-0.5*trace(Psi_inv*W.transpose()*Lambda(beta1)*W)) - 1.5*n*log(2*M_PI) + 1.5*log|Gamma^{-1}| + (n/2)*log|Psi^{-1}|
#
#  likeli_sum = -0.5 * tr(Psi_inv %*% t(W) %*% Lambda(beta1) %*% W)

#  likeli_sum = likeli_sum +
#    ( 3 * sp_log_det_specific(beta1) +
#        n * determinant(Psi_inv, logarithm = T)$modulus - 3*n*log(2*pi) )/2
#
#  return (likeli_sum)
# }
