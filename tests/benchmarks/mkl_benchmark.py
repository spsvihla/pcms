import os
import time
import ctypes
import pcms.haar
import pcms.tree

def get_mkl_threads():
    try:
        mkl_rt = ctypes.CDLL("libmkl_rt.so")
        mkl_get_max_threads = mkl_rt.MKL_Get_Max_Threads
        mkl_get_max_threads.restype = ctypes.c_int
        return mkl_get_max_threads()
    except OSError:
        print("Could not load libmkl_rt.so")
        return None

def benchmark_sparsify():
    DATA = os.environ['DATA']
    GG_DATA = os.path.join(DATA, "greengenes/gg_13_8_otus")
    GG_TREES = os.path.join(GG_DATA, "trees")

    # Print MKL threads reported at runtime
    actual_threads = get_mkl_threads()
    print(f"MKL reports {actual_threads} threads in use.")

    t2 = pcms.tree.nwk2tree(
        os.path.join(GG_TREES, "97_otus_unannotated.nwk"),
        ensure_planted=True
    )
    t0 = time.time()
    cpu0 = time.process_time()
    Q2, S2 = pcms.haar.sparsify(t2)
    t1 = time.time()
    cpu1 = time.process_time()

    print(f"Elapsed wall time: {t1 - t0:.4f}s, CPU time: {cpu1 - cpu0:.4f}s")

if __name__ == "__main__":
    benchmark_sparsify()
