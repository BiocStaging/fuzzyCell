# Parallel loops, run as in spicyR and lisaClust.

# A BiocParallel param for `cores` forked workers (serial on Windows). With a
# `seed`, every element of a loop gets its own random stream, so the results do
# not depend on the number of cores.
parallelParam = function(cores = 1, seed = NULL) {
    BPPARAM = if (cores > 1 && .Platform$OS.type != "windows") BiocParallel::MulticoreParam(cores) else
        BiocParallel::SerialParam()
    if (!is.null(seed)) BiocParallel::bpRNGseed(BPPARAM) = seed
    BPPARAM
}

# The task function and elements of the current parLapply() call, found by
# forked workers in their copy of this environment.
parallelEnv = new.env(parent = emptyenv())
runTask = function(i) parallelEnv$task(parallelEnv$X[[i]])

# lapply(X, FUN) over the workers of BPPARAM (from parallelParam()).
# - The workers are forked once, after the task exists, so they share the
#   session's data without being sent it, and only task numbers and results pass
#   between processes.
# - BLAS runs one thread in each worker during the call. OpenBLAS otherwise
#   starts one thread per core in every worker (128 on a 128-core server), which
#   slows the workers down many times. The setting is changed in this session
#   before the workers start, because changing it inside a worker first restarts
#   every thread.
# - Tasks are handed out in about eight chunks per worker, so that quick and slow
#   tasks balance out.
parLapply = function(X, FUN, BPPARAM) {
    workers = BiocParallel::bpnworkers(BPPARAM)
    # serially, the param still sets each element's random stream
    if (workers <= 1L || length(X) <= 1L) return(BiocParallel::bplapply(X, FUN, BPPARAM = BPPARAM))
    # only when BLAS has several threads: restoring them after forking restarts
    # every thread, so a single-threaded BLAS is left alone
    blasThreads = RhpcBLASctl::blas_get_num_procs()
    if (isTRUE(blasThreads > 1L)) {
        RhpcBLASctl::blas_set_num_threads(1L)
        on.exit(RhpcBLASctl::blas_set_num_threads(blasThreads), add = TRUE)
    }
    # OpenMP likewise: forked workers can hang if this session has used several OpenMP threads
    ompThreads = RhpcBLASctl::omp_get_max_threads()
    if (isTRUE(ompThreads > 1L)) {
        RhpcBLASctl::omp_set_num_threads(1L)
        on.exit(RhpcBLASctl::omp_set_num_threads(ompThreads), add = TRUE)
    }
    BiocParallel::bptasks(BPPARAM) = min(length(X), 8L * workers)
    assign("task", FUN, envir = parallelEnv)
    assign("X", X, envir = parallelEnv)
    on.exit(rm(list = c("task", "X"), envir = parallelEnv), add = TRUE)
    BiocParallel::bpstart(BPPARAM)
    on.exit(BiocParallel::bpstop(BPPARAM), add = TRUE, after = FALSE)
    out = BiocParallel::bplapply(seq_along(X), runTask, BPPARAM = BPPARAM)
    names(out) = names(X)
    out
}
