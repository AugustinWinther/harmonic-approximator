"""Functions for testing other functions.

"""

# Standard library
import time
import pstats
import cProfile
from typing import Any, Callable

# Third party
import numpy as np

# Local
import vis
import calc

def profile(func: Callable, **func_kwargs: Any) -> None:
    """Profiles the given function and saves data in the file perf.prof.
    
    """
    with cProfile.Profile() as pr: func(**func_kwargs)
    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.dump_stats("perf.prof")

def timeit(func: Callable, N_runs: int = 1, **func_kwargs: Any) -> float:
    """Times given function and returns exec. time in seconds. 
    
    If N_runs > 1, then the mean of all the exec. times is returned.
    """
    if N_runs == 1:
        start = time.perf_counter()
        func(**func_kwargs)
        end = time.perf_counter()
        return end - start
    else:
        exec_times = np.empty(N_runs)
        for i in range(N_runs):
            exec_times[i] = timeit(func, **func_kwargs)
        return np.mean(exec_times)

if __name__ == "__main__":
    vis.print_ten_first_energies(N=17, L=11)
