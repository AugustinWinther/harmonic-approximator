"""Functions for profiling and testing other functions performance.

"""

# Standard library
import time
import pstats
import cProfile
from typing import Any, Callable

# Third party
import numpy as np

def profile(func: Callable, **func_kwargs: Any) -> None:
    """Profiles the given function and saves data in the file perf.prof.
    
    """
    with cProfile.Profile() as pr: func(**func_kwargs)
    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.dump_stats("perf.prof")

def timeit(func: Callable, **func_kwargs: Any) -> tuple[float, Any]:
    """Times given function and returns exec. time in ms and func return values.
    
    """
    start = time.perf_counter()
    return_values = func(**func_kwargs)
    end = time.perf_counter()
    exec_time = (end - start)*1000

    return exec_time, return_values

def timeits(N_runs: int, func: Callable, **func_kwargs: Any) -> float:
    """Times given function N_runs amount of times and returns mean exec. time 
    
    """
    exec_times = np.empty(N_runs)
    for i in range(N_runs):
        exec_times[i] = timeit(func, **func_kwargs)[0]
    return np.mean(exec_times)
