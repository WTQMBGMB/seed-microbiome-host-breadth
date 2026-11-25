# =============================================================================
# phylosig_fast.R — Pagel's λ / Blomberg's K 的快速实现（仅依赖 ape）
#
# 动机：phytools::phylosig 在 10k tips 上反复 solve()，单机运行过慢。
#   λ 用 vcv 的相关矩阵特征分解 + 一维 optimize（与 phytools 的
#   likelihoodLambda 等价）；K 的置换检验用矩阵化 BLAS 运算。
#
# 输入约定：C = ape::vcv(tree)，rownames 为 tip labels；x 为命名向量。
# =============================================================================

prepare_vcv <- function(tree, x) {
  C <- ape::vcv(tree)
  x <- x[rownames(C)]
  stopifnot(!anyNA(x), length(x) == nrow(C))
  list(C = C, x = x)
}

pagel_lambda_fast <- function(C, x, niter = 10L) {
  x <- x[rownames(C)]
  n <- nrow(C)
  D <- diag(C)
  stopifnot(all(D > 0))
  s <- sqrt(D)
  A <- C / s
  A <- sweep(A, 2, s, "/")
  ev <- eigen(A, symmetric = TRUE, only.values = FALSE)
  U <- ev$vectors
  lam <- ev$values
  lam[lam < 0] <- 0

  z <- x / s
  one <- 1 / s
  logdetD <- sum(log(D))

  lik <- function(lambda) {
    d <- lambda * lam + (1 - lambda)
    if (any(d <= 0)) return(-Inf)
    Minv_z <- U %*% (as.numeric(crossprod(U, z)) / d)
    Minv_one <- U %*% (as.numeric(crossprod(U, one)) / d)
    # GLS mean in z space: mean vector is a * one (not a * 1), so the
    # numerator/denominator must both be weighted by one.  The previous
    # unweighted sums biased the likelihood whenever tip variances differ.
    a <- sum(one * Minv_z) / sum(one * Minv_one)
    r <- z - a * one
    Minv_r <- U %*% (as.numeric(crossprod(U, r)) / d)
    q <- sum(r * Minv_r)
    sig2 <- q / n
    logdet <- logdetD + sum(log(d))
    -q / (2 * sig2) - n * log(2 * pi) / 2 -
      (logdet + n * log(sig2)) / 2
  }

  # phytools 在 10 个子区间分别 optimize 后取最大，避免似然非单峰时
  # 单区间 optimize 落到局部极大。
  intervals <- cbind(
    seq(0, 1 - 1 / niter, length.out = niter),
    seq(1 / niter, 1, length.out = niter))
  fits <- apply(intervals, 1, stats::optimize, f = lik,
                maximum = TRUE, tol = 1e-6)
  liks <- vapply(fits, function(z) z$objective, numeric(1))
  opt <- fits[[which.max(liks)]]
  logL0 <- lik(0)
  P <- stats::pchisq(2 * (opt$objective - logL0), df = 1,
                     lower.tail = FALSE)
  list(lambda = opt$maximum, logL = opt$objective, logL0 = logL0,
       P = P, lik = lik)
}

# ---------------------------------------------------------------------------
# 剪枝（pruning）版 Pagel λ：O(n) / 每次似然，避免 n×n 特征分解。
#
# Pagel λ 等价于把内部枝长 ×λ、末端枝长调整为
#   t'_i = d_i - λ * (d_i - t_i)
# 以保持 root-to-tip 距离不变，然后在该树上拟合标准 BM。phylolm 的 BM
# 拟合就是 Felsenstein 剪枝，单次 22k tips 约 0.2-0.4 s，而特征分解需要
# 数十分钟。已用直接多元正态似然逐点验证（误差 <1e-6）。
# ---------------------------------------------------------------------------
pagel_lambda_pruning <- function(tree, x, tol = 1e-5) {
  if (!requireNamespace("phylolm", quietly = TRUE)) {
    stop("pagel_lambda_pruning() requires the phylolm package")
  }
  stopifnot(ape::is.rooted(tree))
  x <- x[tree$tip.label]
  stopifnot(!anyNA(x))
  ntip <- length(tree$tip.label)
  d <- ape::node.depth.edgelength(tree)[seq_len(ntip)]
  el <- tree$edge.length
  tip_edge <- tree$edge[, 2] <= ntip
  tipnodes <- tree$edge[tip_edge, 2]
  dat <- data.frame(x = as.numeric(x), row.names = tree$tip.label)
  tr2 <- tree

  neg_loglik <- function(lambda) {
    tr2$edge.length[!tip_edge] <- el[!tip_edge] * lambda
    tr2$edge.length[tip_edge] <-
      d[tipnodes] - lambda * (d[tipnodes] - el[tip_edge])
    # phylolm 拒绝“姐妹末端枝同时为 0”的树；lambda=1 时零末端枝会重现，
    # 加一个可忽略的 eps 以保持数值稳定（不影响 lambda < 1 的结果）。
    if (any(tr2$edge.length <= 0)) {
      tr2$edge.length[tr2$edge.length <= 0] <- 1e-12
    }
    -as.numeric(stats::logLik(
      phylolm::phylolm(x ~ 1, data = dat, phy = tr2, model = "BM")))
  }

  opt <- stats::optimize(neg_loglik, c(0, 1), tol = tol)
  logL0 <- -neg_loglik(0)
  P <- stats::pchisq(2 * (-opt$objective - logL0), df = 1,
                     lower.tail = FALSE)
  list(lambda = opt$minimum, logL = -opt$objective, logL0 = logL0, P = P,
       method = "pruning (phylolm BM on lambda-transformed tree)")
}

blomberg_k_fast <- function(C, x, nsim = 1000L, seed = 20261107L) {
  x <- x[rownames(C)]
  n <- nrow(C)

  # FastTree 会把完全相同的序列折叠成零长度末端枝，VCV 因此可能存在
  # 完全相同的行，solve() 报 "system is exactly singular"。仅在此时加
  # 极小 ridge（1e-8 × 平均对角），对正常矩阵不做任何改动。
  ridge <- 0
  invC <- tryCatch(solve(C), error = function(e) NULL)
  if (is.null(invC)) {
    ridge <- 1e-8 * mean(diag(C))
    C <- C + diag(ridge, n)
    invC <- solve(C)
  }
  sInv <- sum(invC)
  a <- sum(invC %*% x) / sInv
  r <- x - a
  K <- (sum(r^2) / sum(r * (invC %*% r))) /
    ((sum(diag(C)) - n / sInv) / (n - 1))

  if (nsim <= 0L) {
    return(list(K = as.numeric(K), P = NA_real_, sim.K = numeric(0),
                ridge = ridge))
  }

  set.seed(seed)
  Xp <- replicate(nsim, sample(x))
  Y <- invC %*% Xp
  a_p <- colSums(Y) / sInv
  Xc <- Xp - rep(a_p, each = n)
  Yc <- Y - rep(a_p * sInv, each = n)
  num <- colSums(Xc^2)
  den <- colSums(Xc * Yc)
  const <- (sum(diag(C)) - n / sInv) / (n - 1)
  Kp <- (num / den) / const
  P <- mean(Kp >= K)
  list(K = as.numeric(K), P = P, sim.K = Kp, ridge = ridge)
}
