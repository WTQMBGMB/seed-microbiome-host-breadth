# =============================================================================
# phylo_d_fast.R — Fritz & Purvis D for binary traits（快速实现）
#
# caper::phylo.d 在 10k+ tips 上极慢：VCV.array 用 R 循环（10,655 tips 需
# 75 s），且 rmvnorm() 默认用特征分解。本实现与 caper 的 D 定义完全一致，
# 但：
#   - VCV 用 ape::vcv()（C 实现，0.3 s）
#   - Brownian 模拟用 Cholesky（chol + crossprod）
#   - 对比计算仍调用 caper:::contrCalc(picMethod="phylo.d")
#
# 依赖：ape, caper（仅内部 contrCalc）
#
# D 的解释（Fritz & Purvis 2010）：
#   D = 0  -> 符合 Brownian 阈值模型（强系统发育信号）
#   D = 1  -> 随机（无系统发育结构）
#   D < 0  -> 比 Brownian 更聚集；D > 1 -> 比随机更分散
#   P_random = P(D = 1)（随机零模型；显著小 => 有信号）
#   P_BM     = P(D = 0)（Brownian 零模型）
# =============================================================================

phylo_d_fast <- function(phy, ds, permut = 1000L, seed = 20261108L,
                         eps = 1e-8, obs_only = FALSE) {
  stopifnot(ape::is.rooted(phy), ape::is.binary(phy))
  n <- length(phy$tip.label)
  ds <- as.numeric(ds)
  stopifnot(length(ds) == n, !anyNA(ds), length(unique(ds)) == 2)
  prop0 <- mean(ds == 0)

  if (any(phy$edge.length <= 0)) {
    if (eps <= 0) stop("tree has zero-length branches; set eps > 0")
    phy$edge.length[phy$edge.length <= 0] <- eps
  }

  phy2 <- ape::reorder.phylo(phy, "pruningwise")
  phy2$node.label <- as.character((n + 1):(n + phy2$Nnode))

  # 与 caper:::contrCalc(picMethod="phylo.d") 等价，但去掉了逐节点
  # data.frame/complete.cases 开销；二元树每个内部节点：
  #   node = (v1/b1 + v2/b2) / (1/b1 + 1/b2)
  #   contrast = |v1 - node| + |v2 - node|
  # 对所有列（观测 + 置换 + Brownian 模拟）同时计算。
  child <- split(phy2$edge[, 2], phy2$edge[, 1])
  # 用 postorder 保证子节点先于父节点（ape 的节点编号不保证 parent < child）
  phy_post <- ape::reorder.phylo(phy2, "postorder")
  nds <- unique(phy_post$edge[, 1])
  edge_child <- phy2$edge[, 2]

  d_metric <- function(M) {
    P <- ncol(M)
    vals <- matrix(0, nrow = n + phy2$Nnode, ncol = P + 1L)
    vals[seq_len(n), 1L] <- ds
    if (P > 0L) vals[seq_len(n), seq_len(P) + 1L] <- M
    total <- numeric(P + 1L)
    el <- phy2$edge.length  # 局部副本；contrCalc 会沿路调整枝长
    for (nd in nds) {
      ch <- child[[as.character(nd)]]
      if (length(ch) != 2L) stop("binary tree required")
      idx <- match(ch, edge_child)
      bl <- el[idx]
      w <- 1 / bl
      v <- (vals[ch[1], ] * w[1] + vals[ch[2], ] * w[2]) / sum(w)
      vals[nd, ] <- v
      total <- total + abs(vals[ch[1], ] - v) + abs(vals[ch[2], ] - v)
      # caper:::contrCalc 的 PIC 枝长调整：父枝 += prod(bl)/sum(bl)
      pidx <- match(nd, edge_child)
      if (!is.na(pidx)) el[pidx] <- el[pidx] + prod(bl) / sum(bl)
    }
    total
  }

  obssocc <- d_metric(matrix(numeric(0), nrow = n))[1]
  if (obs_only) {
    return(list(D = NA_real_, P_random = NA_real_, P_BM = NA_real_,
                observed = unname(obssocc), n = n, permut = 0L,
                eps = eps))
  }

  set.seed(seed)
  ds.ran <- replicate(permut, sample(ds))
  ransocc <- d_metric(ds.ran)[-1]

  C <- ape::vcv(phy2)
  L <- chol(C)
  Z <- matrix(stats::rnorm(n * permut), nrow = n)
  sim <- crossprod(L, Z)
  thr <- apply(sim, 2, stats::quantile, prop0)
  ds.phy <- sweep(sim, 2, thr, "<")
  storage.mode(ds.phy) <- "integer"
  physocc <- d_metric(ds.phy)[-1]

  D <- (obssocc - mean(physocc)) / (mean(ransocc) - mean(physocc))
  list(D = unname(D),
       P_random = mean(ransocc < obssocc),
       P_BM = mean(physocc > obssocc),
       observed = unname(obssocc),
       mean_random = mean(ransocc),
       mean_brownian = mean(physocc),
       n = n, permut = as.integer(permut), eps = eps)
}
