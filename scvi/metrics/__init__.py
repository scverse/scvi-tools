from scvi.metrics.clustering import get_latent
from .classification import compute_accuracy
from .log_likelihood import compute_log_likelihood
from .visualization import show_t_sne

__all__ = ['compute_log_likelihood',
           'compute_accuracy',
           'show_t_sne']


class Metric:
    def __init__(self, mode):
        self.mode = mode

class LLMetric(Metric):
    def __init__(self):
        super(LLMetric, self).__init__(mode='min')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        ll = compute_log_likelihood(infer.model, data_loader, *args, **kwargs)
        print("LL for %s is : %.4f"%(name, ll))
        return ll

class AccuracyMetric(Metric):
    def __init__(self):
        super(AccuracyMetric, self).__init__(mode='max')

    def __call__(self, infer, data_loader, name, *args, **kwargs):
        acc = compute_accuracy(infer, data_loader, *args, **kwargs)
        print("Acc for %s is : %.4f"%(name, acc))
        return acc

# class Task:
#     pass
#
# class ClusteringTask(Task):
#     pass
#
# class ImputationTask(Task):
#     pass
#
# class DEStatsTask(Task):
#     pass
