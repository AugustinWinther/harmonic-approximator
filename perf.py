"""Functions for profiling and testing other functions performance.

"""

# Standard library
import pstats
import cProfile
from typing import Any, Callable

def profile_func(func: Callable, **func_kwargs: Any) -> None:
    """Profiles the given function and saves data in the file perf.prof.
    
    """
    with cProfile.Profile() as pr: func(**func_kwargs)
    stats = pstats.Stats(pr)
    stats.sort_stats(pstats.SortKey.TIME)
    stats.dump_stats("perf.prof")
