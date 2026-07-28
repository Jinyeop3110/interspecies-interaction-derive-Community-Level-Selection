"""Raw ASV tables to processed data tables.

Python port of the MATLAB stage driven by ``Main_synthetic.m`` and
``Main_natural.m``.  Every ported stage reproduces its archived output to
floating-point tolerance; see ``tests/test_pipeline_parity.py``.

Stages
------
``sequences``           ASV tables -> ``processed_Sequences_*.xlsx``
``coalescence_events``  -> ``processed_CoalescenceEvent_*.xlsx``
``communities``         -> ``processed_Communities_*.xlsx``
``summary_metrics``     the metric definitions the two stages above share
"""

from . import coalescence_events, communities, sequences, summary_metrics

__all__ = ["sequences", "coalescence_events", "communities", "summary_metrics"]
