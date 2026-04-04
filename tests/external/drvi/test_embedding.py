import numpy as np
import pandas as pd
import torch

from scvi.external.drvi.nn_modules.embedding import FeatureEmbedding


class TestFeatureEmbedding:
    def make_test_data(self):
        feature_table = pd.DataFrame(
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
        feature_table_row_extended = pd.DataFrame(
            [
                ["g3", "g3", "GE"],
                ["g4", "g4", "GE"],
                ["p4", "g9", "PK"],
                ["p6", "g3", "PK"],
            ],
            columns=["id", "nearest_gene", "feature_type"],
        )
        feature_table_full_extended = pd.DataFrame(
            [
                ["g3", "g3", "GE", "NEW_FEATURE_1"],
                ["p1", "g1", "KK", "NEW_FEATURE_2"],
                ["p100", "g100", "PK", "NEW_FEATURE_2"],
            ],
            columns=["id", "nearest_gene", "feature_type", "new_feature"],
        )

        return feature_table, feature_table_row_extended, feature_table_full_extended

    def test_feature_embedding_simple(self):
        sample_df = self.make_test_data()[0]
        feature_embedding = FeatureEmbedding.from_pandas_dataframe(sample_df, [2, 2, 1])
        accessed_data = feature_embedding(
            np.array(
                [
                    ["g3", "g3", "GE"],
                    ["p1", "g3", "PK"],
                ]
            )
        )
        assert ((accessed_data[0, 2:4] - accessed_data[1, 2:4]) ** 2).sum() < 1e-8  # g3 is same

    def test_feature_embedding_cache(self):
        sample_df = self.make_test_data()[0]
        feature_embedding = FeatureEmbedding.from_pandas_dataframe(sample_df, [2, 2, 1])
        accessed_data = feature_embedding(
            np.array(
                [
                    ["g3", "g3", "GE"],
                    ["p1", "g3", "PK"],
                ]
            ),
            index_cache_key="x",
        )
        assert ((feature_embedding(None, index_cache_key="x") - accessed_data) ** 2).sum() < 1e-8

    def test_feature_embedding_row_extension(self):
        sample_df, row_extension_df, _ = self.make_test_data()
        union_df = pd.concat([sample_df, row_extension_df])

        feature_embedding = FeatureEmbedding.from_pandas_dataframe(sample_df, [2, 2, 1])
        feature_embedding_extended = FeatureEmbedding.from_pandas_dataframe(union_df, [2, 2, 1])

        key_array = np.array(
            [
                ["g1", "g3", "GE"],
                ["p4", "g9", "PK"],
            ]
        )
        accessed_data_pre_transfer = feature_embedding_extended(key_array)
        feature_embedding_extended.load_weights_from_trained_module(feature_embedding)
        accessed_data_post_transfer = feature_embedding_extended(key_array)

        # First row should come from feature_embedding (overwritten)
        assert (
            (feature_embedding(key_array[0:1, :]) - accessed_data_post_transfer[0:1, :]) ** 2
        ).sum() < 1e-8
        # In second row, 'p4', 'g9' should not change
        assert (
            (accessed_data_pre_transfer[1, :4] - accessed_data_post_transfer[1, :4]) ** 2
        ).sum() < 1e-8

    def test_feature_embedding_full_extension(self):
        sample_df, _, full_extension_df = self.make_test_data()
        union_df = pd.concat([sample_df, full_extension_df])

        feature_embedding = FeatureEmbedding.from_pandas_dataframe(sample_df, [2, 2, 1])
        feature_embedding_extended = FeatureEmbedding.from_pandas_dataframe(union_df, [2, 2, 2, 2])

        key_array = np.array(
            [
                ["g1", "g3", "GE", "NEW_FEATURE_1"],
                ["p100", "g100", "KK", "NEW_FEATURE_2"],
            ]
        )
        accessed_data_pre_transfer = feature_embedding_extended(key_array)
        feature_embedding_extended.load_weights_from_trained_module(feature_embedding)
        accessed_data_post_transfer = feature_embedding_extended(key_array)

        # 'g1', 'g3', 'GE (dim1)' should come from feature_embedding (overwritten)
        assert (
            (feature_embedding(key_array[0:1, :3]) - accessed_data_post_transfer[0:1, :-3]) ** 2
        ).sum() < 1e-8
        # Last 3 dims (1 for 'feature_type' and 2 for 'new_feature') should not change (new)
        assert (
            (accessed_data_pre_transfer[:, -3:] - accessed_data_post_transfer[:, -3:]) ** 2
        ).sum() < 1e-8
        # Second row should not change (new)
        assert (
            (accessed_data_pre_transfer[1, :] - accessed_data_post_transfer[1, :]) ** 2
        ).sum() < 1e-8

    def test_feature_embedding_freeze_on_extension(self):
        sample_df, _, full_extension_df = self.make_test_data()
        union_df = pd.concat([sample_df, full_extension_df])

        feature_embedding = FeatureEmbedding.from_pandas_dataframe(sample_df, [2, 2, 1])
        feature_embedding_extended = FeatureEmbedding.from_pandas_dataframe(union_df, [2, 2, 2, 2])

        key_array = np.array(
            [
                ["g1", "g3", "GE", "NEW_FEATURE_1"],
                ["p100", "g100", "KK", "NEW_FEATURE_2"],
            ]
        )
        feature_embedding_extended.load_weights_from_trained_module(
            feature_embedding, freeze_old=True
        )
        # grad(x^2/2) = 2 -> grad = data if not freeze
        loss = sum([(p**2 / 2).sum() for _, p in feature_embedding_extended.named_parameters()])
        loss.backward()
        with torch.no_grad():
            for _name, p in list(feature_embedding_extended.named_parameters()):
                p.data -= p.grad  # everything that is not freeze should become zero

        accessed_data_post_transfer = feature_embedding_extended(key_array)
        # 'g1', 'g3', 'GE (dim1)' should come from feature_embedding (overwritten)
        assert (
            (feature_embedding(key_array[0:1, :3]) - accessed_data_post_transfer[0:1, :-3]) ** 2
        ).sum() < 1e-8
        # Last 3 dims should become zero
        assert ((accessed_data_post_transfer[:, -3:]) ** 2).sum() < 1e-8
        # Second row should become zero
        assert ((accessed_data_post_transfer[1, :]) ** 2).sum() < 1e-8
