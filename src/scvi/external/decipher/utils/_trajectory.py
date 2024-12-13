from typing import Any

import numpy as np


class Trajectory:
    """A class that defines a trajectory through a latent space.

    This class creates a trajectory through a latent space by connecting a sequence of points
    (e.g., cluster centroids) with straight lines. The trajectory is parameterized by time,
    where time 0 corresponds to the first point and increases linearly with distance along
    the trajectory.

    Parameters
    ----------
    rep_key : str
        Obsm key for the latent representation being used
    cluster_locations : np.ndarray
        Array of coordinates for points to connect, shape (n_points, n_dimensions)
    cluster_ids : List[str]
        List of cluster IDs corresponding to the locations
    point_density : int, default=50
        Number of interpolated points per unit distance in the latent space

    Attributes
    ----------
    trajectory_latent : np.ndarray
        Coordinates of points along the interpolated trajectory
    trajectory_time : np.ndarray
        Time values corresponding to each point in trajectory_latent
    cumulative_length : np.ndarray
        Cumulative distance along the trajectory at each cluster location
    cluster_locations : np.ndarray
        Original cluster locations used to create the trajectory
    n_points : int
        Total number of points in the interpolated trajectory
    """

    def __init__(
        self,
        rep_key: str,
        cluster_locations: np.ndarray,
        cluster_ids: list[str],
        point_density: int = 50,
    ) -> None:
        self._point_density = point_density
        self.cluster_ids = cluster_ids
        cluster_locations = np.array(cluster_locations)
        distances = []
        for s, e in zip(cluster_locations, cluster_locations[1:], strict=False):
            v = e - s
            d = np.linalg.norm(v)
            distances.append(d)

        cumulative_length = np.cumsum(distances)

        self.cumulative_length = cumulative_length
        self.cluster_locations = cluster_locations
        self.n_points = int(self._point_density * self.cumulative_length[-1])
        self.rep_key = rep_key

        self.trajectory_latent, self.trajectory_time = self._linspace(self.n_points)

    def _linspace(self, num: int = 100) -> tuple[np.ndarray, np.ndarray]:
        total_length = self.cumulative_length[-1]
        times = np.linspace(0, total_length, num)
        res = []
        for t in times:
            res.append(self.at_time(t))
        trajectory_latent = np.array(res).astype(np.float32)
        trajectory_time = times

        return trajectory_latent, trajectory_time

    def at_time(self, t: float) -> np.ndarray:
        """Get the coordinates at a specific time point along the trajectory.

        Parameters
        ----------
        t : float
            Time point along the trajectory.

        Returns
        -------
        np.ndarray
            Coordinates at the specified time point.
        """
        i = 0
        while t > self.cumulative_length[i]:
            i += 1
        if i > 0:
            t = (t - self.cumulative_length[i - 1]) / (
                self.cumulative_length[i] - self.cumulative_length[i - 1]
            )
        else:
            t = t / self.cumulative_length[i]

        return self.cluster_locations[i] * (1 - t) + t * self.cluster_locations[i + 1]

    def to_dict(self) -> dict[str, Any]:
        """Convert trajectory object to dictionary representation.

        Returns
        -------
        Dict[str, Any]
            Dictionary containing trajectory data including cluster locations,
            cluster IDs, trajectory points, times, point density and representation key.
        """
        return {
            "cluster_locations": self.cluster_locations,
            "cluster_ids": self.cluster_ids,
            "points": self.trajectory_latent,
            "times": self.trajectory_time,
            "density": self._point_density,
            "rep_key": self.rep_key,
        }

    @staticmethod
    def from_dict(d: dict[str, Any]) -> "Trajectory":
        """Create a Trajectory object from a dictionary representation.

        Parameters
        ----------
        d : Dict[str, Any]
            Dictionary containing trajectory data.

        Returns
        -------
        Trajectory
            New Trajectory object initialized with data from dictionary.
        """
        trajectory = Trajectory(
            rep_key=d["rep_key"],
            cluster_locations=d["cluster_locations"],
            cluster_ids=d["cluster_ids"],
            point_density=d.get("density", None),
        )
        return trajectory
