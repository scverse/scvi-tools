import numpy as np
import pandas as pd

from scvi.external.drvi.nn_modules.encodig import FeatureOneHotEncoding


class TestFeatureEmbedding:
    def make_test_data(self):
        sample_df = pd.DataFrame(
            [
                ["g1", "g1", "GE"],
                ["g2", "g2", "GE"],
                ["g3", "g3", "GE"],
                ["p1", "g1", "PK"],
                ["p2", "g1", "PK"],
                ["p3", "g4", "PK"],
            ],
            columns=["id", "nearest_gene", "feature_type"],
        )
        return sample_df

    def test_feature_encoding_simple(self):
        sample_df = self.make_test_data()
        feature_encoding = FeatureOneHotEncoding.from_pandas_dataframe(sample_df)
        print(feature_encoding)
        accessed_data = feature_encoding(
            np.array(
                [
                    ["g3", "g3", "GE"],
                    ["g3", "g3", "PK"],
                    ["g3", "g3", "GE"],
                ]
            ),
            index_cache_key="x",
        )
        assert ((accessed_data[0] - accessed_data[1]) ** 2).sum() >= 1
        assert ((accessed_data[0] - accessed_data[2]) ** 2).sum() < 1e-8

