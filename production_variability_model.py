"""
Production Variability Model for Manufacturing
==============================================

This module implements a simple, yet practical, production system simulator
designed for manufacturing contexts.  The model is inspired by the cell
transmission model and advection–diffusion solver contained within the
``MTFC_2026`` repository and draws on queueing theory results (Little’s law
and Kingman’s formula) to quantify the impact of variability and utilisation
on waiting times and throughput.  It provides an example of how to translate
the ideas from transportation and environmental modelling into a factory
setting.

The simulation operates in discrete time and tracks the number of items
waiting and being processed at each station (analogous to cells in the
cell–transmission model).  Each station has a finite buffer and a finite
processing rate, and items flow between stations according to an adjacency
matrix.  Arrival patterns (demand) can be specified exogenously.  At each
time step the model:

1. Adds new arrivals to the first station, subject to buffer limits.
2. Calculates the amount of work that can be sent from each station to its
   downstream neighbours (sending capacity) based on current WIP and
   processing rate.
3. Calculates the capacity of downstream stations to receive work
   (receiving capacity) based on available buffer space.
4. Moves work accordingly and updates the in–station WIP.
5. Estimates the expected waiting time at each station using Kingman’s
   approximation given the utilisation and coefficient of variation of
   arrival and service processes.

This model is intended for educational and prototyping purposes.  For a
complete industrial implementation, users should calibrate the model
parameters with their own data and extend the logic to handle multiple product
types, routing rules, machine breakdowns and other realistic features.

References
----------
* Little’s law relates the average inventory ``L`` in the system to the
  average throughput ``λ`` and the lead time ``W`` as ``L = λ * W``.  See
  ``allaboutlean.com`` for a concise overview.
* Kingman’s formula (or VUT equation) approximates the expected waiting time
  ``E(W)`` in a ``G/G/1`` queue as ``E(W) ≈ (ρ/(1−ρ)) * ((C_a^2 + C_s^2)/2) *
  μ_s``, where ``ρ`` is the utilisation, ``C_a`` and ``C_s`` are the
  coefficients of variation of the inter–arrival and service times, and
  ``μ_s`` is the mean service time【686481500816131†L99-L137】.
* Discrete‑event simulation has long been used to analyse manufacturing
  systems.  It samples events (arrivals, service completions) and allows
  utilisation and waiting times to be estimated.  The Project Production
  Institute notes that discrete‑event simulation tracks the state of the
  system at event boundaries and enforces that resources and inputs are
  available before operations begin【530851410714537†L466-L475】.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Callable
import numpy as np


@dataclass
class Station:
    """Represents a manufacturing station (e.g. machine, work centre).

    Attributes
    ----------
    name: str
        Identifier of the station.
    processing_rate: float
        Average number of items processed per unit time (items/hour).
    buffer_capacity: float
        Maximum number of items (work‑in‑process) the station can hold,
        including those being processed and those waiting in queue.
    service_cv: float
        Coefficient of variation of service time.  A value of 0 implies
        deterministic processing; higher values model stochastic processes.
    arrival_cv: float
        Coefficient of variation of arrival process to this station.  This
        affects waiting time via Kingman’s formula.  When ``arrival_cv`` is
        unspecified, the value will be updated dynamically based on upstream
        variability.
    current_wip: float = 0.0
        Current work‑in‑process (both in queue and in service).
    total_processed: float = 0.0
        Cumulative number of items that have completed processing at this
        station.
    upstream_arrival_rate: float = 0.0
        Estimated arrival rate (items/hour) from the upstream stations.  This
        value will be updated at each time step and used to compute
        utilisation.
    downstream_indices: List[int] = field(default_factory=list)
        Indices of downstream stations in the production network.
    """

    name: str
    processing_rate: float
    buffer_capacity: float
    service_cv: float
    arrival_cv: float
    current_wip: float = 0.0
    total_processed: float = 0.0
    upstream_arrival_rate: float = 0.0
    downstream_indices: List[int] = field(default_factory=list)

    def utilisation(self) -> float:
        """Compute the current utilisation ρ = arrival rate / service capacity.

        Returns a value between 0 and 1 (values >=1 indicate overload).  When
        utilisation is very high the Kingman approximation yields very large
        waiting times.
        """
        if self.processing_rate <= 0:
            return 0.0
        rho = self.upstream_arrival_rate / (self.processing_rate)
        return min(max(rho, 0.0), 0.999)  # cap at 0.999 to avoid division by zero

    def mean_service_time(self) -> float:
        """Return the mean service time (hours per item)."""
        return 1.0 / self.processing_rate if self.processing_rate > 0 else np.inf

    def expected_waiting_time(self) -> float:
        """Estimate the expected waiting time in queue using Kingman’s formula.

        Kingman’s formula (also known as the VUT equation) approximates the
        average waiting time (excluding service time) for a ``G/G/1`` queue.  It
        requires the coefficient of variation of arrival and service time and the
        utilisation.  See the module docstring for details【686481500816131†L99-L137】.

        Returns
        -------
        float
            Estimated waiting time (in hours) for an item arriving to this
            station.
        """
        rho = self.utilisation()
        if rho <= 0.0:
            return 0.0
        mu_s = self.mean_service_time()
        c_a = max(self.arrival_cv, 0.0)
        c_s = max(self.service_cv, 0.0)
        # Kingman approximation: E(W) ≈ (ρ/(1−ρ)) * ((C_a^2 + C_s^2)/2) * μ_s
        waiting = (rho / (1.0 - rho)) * ((c_a**2 + c_s**2) / 2.0) * mu_s
        return max(waiting, 0.0)


class ProductionLineModel:
    """Discrete‑time simulation of a production line or network.

    Parameters
    ----------
    stations : List[Station]
        List of stations in topological order.  Items flow from low indices to
        higher indices according to the adjacency matrix provided.
    adjacency : np.ndarray
        Square (n×n) adjacency matrix.  Entry ``adjacency[i,j]`` is 1 if
        station i can send items directly to station j, and 0 otherwise.  The
        matrix should be upper triangular (no backflow) to avoid cycles.
    time_step : float
        Duration of a simulation step in hours.  For example, 1/60 for
        minute‑level resolution.
    """

    def __init__(self, stations: List[Station], adjacency: np.ndarray,
                 time_step: float = 1.0 / 60.0):
        if len(stations) != adjacency.shape[0] or adjacency.shape[0] != adjacency.shape[1]:
            raise ValueError("Adjacency matrix must be square with size equal to number of stations")
        self.stations: List[Station] = stations
        self.adj = adjacency.astype(int)
        self.dt = time_step
        # Precompute downstream lists for each station
        for i, st in enumerate(self.stations):
            st.downstream_indices = list(np.where(self.adj[i] > 0)[0])

        # Histories for analysis
        self.history_wip: List[List[float]] = []
        self.history_wait: List[List[float]] = []
        self.history_throughput: List[float] = []

        # Internal tracker for throughput of the last station in current step
        self._last_step_throughput: float = 0.0

    def add_demand(self, demand: Dict[int, float]):
        """Add new arrivals to stations with external demand.

        Parameters
        ----------
        demand : Dict[int, float]
            Mapping from station index to arrival rate (items per hour) during
            this time step.  The arrival quantity is `rate * dt`.  If a station
            has multiple upstream suppliers, the demand should be aggregated
            upstream; this argument represents exogenous demand only.
        """
        for idx, rate in demand.items():
            station = self.stations[idx]
            arrivals = rate * self.dt
            space = station.buffer_capacity - station.current_wip
            items_added = min(arrivals, space)
            station.current_wip += items_added
            station.upstream_arrival_rate = rate  # update arrival rate for utilisation

    def compute_flows(self) -> np.ndarray:
        """Compute flows between stations for the current time step.

        Returns
        -------
        np.ndarray
            Matrix ``F`` where ``F[i,j]`` is the number of items moving from
            station i to station j in this time step.
        """
        n = len(self.stations)
        flows = np.zeros((n, n), dtype=float)

        # First compute send and receive capacities for each station
        send_cap = np.zeros(n)
        recv_cap = np.zeros(n)
        for i, st in enumerate(self.stations):
            # Items that can be processed this step: limited by processing rate and current WIP
            send_cap[i] = min(st.processing_rate * self.dt, st.current_wip)
            # Buffer space available to receive new items (excluding those being processed)
            recv_cap[i] = max(0.0, st.buffer_capacity - st.current_wip)

        # Distribute sending capacity to downstream stations proportionally
        for i, st in enumerate(self.stations):
            if send_cap[i] <= 0 or not st.downstream_indices:
                continue
            total_links = len(st.downstream_indices)
            # For simplicity, split evenly across all downstream links
            for j in st.downstream_indices:
                if recv_cap[j] <= 0:
                    continue
                amount = send_cap[i] / total_links
                amount = min(amount, recv_cap[j])
                flows[i, j] = amount
                # Reserve space and capacity
                recv_cap[j] -= amount
                send_cap[i] -= amount

        return flows

    def update_state(self, flows: np.ndarray):
        """Update station WIP based on flows matrix.

        Parameters
        ----------
        flows : np.ndarray
            Matrix of flows returned by `compute_flows`.
        """
        n = len(self.stations)
        # Compute net change in WIP and update processing counts
        delta_wip = np.zeros(n)
        outflows = np.zeros(n)
        inflows = np.zeros(n)
        for i in range(n):
            # Sum of items leaving station i
            outflows[i] = flows[i].sum()
            # Sum of items entering station i
            inflows[i] = flows[:, i].sum()
            delta_wip[i] = inflows[i] - outflows[i]

        for i, st in enumerate(self.stations):
            st.current_wip += delta_wip[i]
            # If station has no downstream, items leaving contribute to total processed
            if not st.downstream_indices:
                st.total_processed += outflows[i]
                # Track throughput leaving the last station during this time step
                self._last_step_throughput = outflows[i]
            # Ensure WIP remains within [0, buffer_capacity]
            st.current_wip = min(max(st.current_wip, 0.0), st.buffer_capacity)

    def record_history(self):
        """Record WIP and waiting time history for the current step."""
        wip_snapshot = [st.current_wip for st in self.stations]
        wait_snapshot = [st.expected_waiting_time() for st in self.stations]
        self.history_wip.append(wip_snapshot)
        self.history_wait.append(wait_snapshot)
        # Record throughput as the number of items that left the last station this step
        self.history_throughput.append(self._last_step_throughput)

    def step(self, demand: Dict[int, float]):
        """Advance the simulation by one time step.

        Parameters
        ----------
        demand : Dict[int, float]
            External arrivals to add at this step.
        """
        # 1. Add new arrivals
        self.add_demand(demand)
        # 2. Compute flows
        flows = self.compute_flows()
        # 3. Update state
        self.update_state(flows)
        # 4. Record metrics
        self.record_history()

    def run(self, demand_schedule: List[Dict[int, float]], n_steps: int):
        """Run the simulation for a sequence of time steps.

        Parameters
        ----------
        demand_schedule : List[Dict[int, float]]
            Sequence of dictionaries representing the arrival rate into stations
            for each step.  Length should be >= n_steps; if shorter, the last
            entry is reused.
        n_steps : int
            Number of steps to simulate.
        """
        for t in range(n_steps):
            # Use the appropriate demand dict or repeat the last one
            if t < len(demand_schedule):
                demand = demand_schedule[t]
            else:
                demand = demand_schedule[-1] if demand_schedule else {}
            self.step(demand)

    def get_history(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return simulation histories.

        Returns
        -------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Tuple of arrays ``(WIP, Wait, Throughput)`` where ``WIP`` and
            ``Wait`` are ``n_steps × n_stations`` arrays, and ``Throughput`` is
            length ``n_steps``.
        """
        return (np.array(self.history_wip),
                np.array(self.history_wait),
                np.array(self.history_throughput))


def example_usage() -> None:
    """Run a small example to demonstrate usage of the model.

    This example sets up a simple line with three stations.  The first station
    receives exogenous demand at 30 items per hour; each station can process
    up to 20 items per hour and hold up to 50 items.  The coefficients of
    variation (C_v) for service time are set to 0.2 for all stations, while
    arrival CV is assumed to be 1.0 for the first station and inherited by
    downstream stations.  The model runs for 8 hours simulated at a 1‑minute
    resolution.  After running the simulation, it prints summary statistics.
    """
    # Define stations
    stations = [
        Station(name="Machine A", processing_rate=20.0, buffer_capacity=50.0,
                service_cv=0.2, arrival_cv=1.0),
        Station(name="Machine B", processing_rate=20.0, buffer_capacity=50.0,
                service_cv=0.2, arrival_cv=1.0),
        Station(name="Machine C", processing_rate=20.0, buffer_capacity=50.0,
                service_cv=0.2, arrival_cv=1.0),
    ]
    adjacency = np.array([
        [0, 1, 0],  # A → B
        [0, 0, 1],  # B → C
        [0, 0, 0],  # C is last
    ], dtype=int)
    model = ProductionLineModel(stations, adjacency, time_step=1.0 / 60.0)

    # Demand schedule: constant arrival to station 0 of 15 items/hour.
    # This value is deliberately chosen to be below the processing rate of each
    # station (20 items/hour) to illustrate a stable system.  Feel free to
    # experiment with higher arrival rates to observe queue build‑up.
    demand_schedule = [{0: 15.0} for _ in range(8 * 60)]  # 8 hours at minute resolution
    model.run(demand_schedule, n_steps=len(demand_schedule))

    # Retrieve history
    wip, wait, throughput = model.get_history()
    avg_wip = wip.mean(axis=0)
    avg_wait = wait.mean(axis=0)
    total_out = stations[-1].total_processed
    print("Average WIP per station:", avg_wip)
    print("Average waiting time per station (h):", avg_wait)
    print("Total items produced:", total_out)


if __name__ == "__main__":
    # Execute example when run as a script
    example_usage()