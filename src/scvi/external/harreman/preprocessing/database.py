import time
from itertools import zip_longest
from re import compile, match
from typing import Literal

import pandas as pd
import pooch
from anndata import AnnData

from scvi.external.harreman._data import harreman_data_hash, harreman_data_url

IMPORT_METAB_KEY = "IMPORT"
EXPORT_METAB_KEY = "EXPORT"
BOTH_METAB_KEY = "IMPORT_EXPORT"
CELLCHAT_DB_DIR = "CellChatDB"
HARREMAN_DB_DIR = "HarremanDB"


def extract_interaction_db(
    adata: AnnData,
    use_raw: bool = False,
    species: Literal["mouse", "human"] | None = None,
    database: Literal["transporter", "LR", "both"] | None = None,
    extracellular_only: bool | None = True,
    verbose: bool | None = False,
) -> None:
    """Extract the metabolite transporter or ligand-receptor (LR) database from .csv files.

    Parameters
    ----------
    adata
        AnnData object to compute database for.
    use_raw
        Whether to use adata.raw.X for database computation.
    species
        Species identity to select the LR database from CellChatDB.
    database
        Whether to use the transporter database, the LR database, or both.
    extracellular_only
        Whether to restrict the communication inference to extracellular metabolites.
    verbose
        Whether to print progress and status messages.

    Returns
    -------
    Genes by metabolites (or LRs) dataframe. Index is aligned to genes from adata.

    """
    start = time.time()
    if verbose:
        print("Extracting interaction database...")

    if species not in {"mouse", "human"}:
        raise ValueError(f"Unsupported species: {species}. Choose 'mouse' or 'human'.")
    if database not in {"transporter", "LR", "both"}:
        raise ValueError("Choose one of: 'transporter', 'LR', or 'both'.")

    adata.uns["species"] = species
    index = adata.raw.var.index if use_raw else adata.var_names
    df_list = []

    if database in ["LR", "both"]:
        extract_lr_pairs(adata, species)
        lr_data = build_LR_matrix(
            index, adata.uns["LR_database"], adata.uns["ligand"], adata.uns["receptor"]
        )
        df_list.append(lr_data)

    if database in ["transporter", "both"]:
        metab_dict = extract_transporter_info(adata, species, extracellular_only)
        metab_data = build_transporter_matrix(index, metab_dict)
        df_list.append(metab_data)

    database_df = pd.concat(df_list, axis=1).fillna(0)

    adata.uns["database_varm_key"] = "database"
    adata.uns["database"] = database
    adata.varm["database"] = database_df

    if verbose:
        print("Finished extracting interaction database in %.3f seconds" % (time.time() - start))

    return


def build_LR_matrix(index, database, ligands, receptors):
    """Build a gene-by-ligand-receptor interaction matrix."""
    matrix = pd.DataFrame(0, index=index, columns=database.index)
    matrix.index = matrix.index.str.lower()
    for col in matrix.columns:
        for _key, df, sign in [("ligand", ligands, 1), ("receptor", receptors, -1)]:
            genes = df.loc[col].dropna().astype(str).str.lower()
            genes = genes[genes.isin(matrix.index)]
            matrix.loc[genes, col] = sign
    matrix.index = index
    return matrix.loc[:, matrix.any()]


def build_transporter_matrix(index, metab_dict):
    """Build a gene-by-metabolite transporter direction matrix."""
    matrix = pd.DataFrame(0, index=index, columns=metab_dict.keys())
    matrix.index = matrix.index.str.lower()
    for metab, gene_dirs in metab_dict.items():
        for direction, sign in [
            (IMPORT_METAB_KEY, -1),
            (EXPORT_METAB_KEY, 1),
            (BOTH_METAB_KEY, 2),
        ]:
            genes = pd.Index(gene_dirs.get(direction, [])).str.lower()
            genes = genes.intersection(matrix.index)
            matrix.loc[genes, metab] = sign
    matrix.index = index
    return matrix.loc[:, matrix.any()]


def extract_transporter_info(
    adata: AnnData,
    species: Literal["mouse", "human"],
    extracellular_only: bool = True,
    export_suffix: str = "(_exp|_export)",
    import_suffix: str = "(_imp|_import)",
    verbose: bool = False,
) -> dict[str, dict[str, list[str]]]:
    """Read csv file to extract the metabolite database."""
    cache = pooch.os_cache("scvi_harreman")

    filenames = {
        "extracellular": pooch.retrieve(
            url=harreman_data_url(HARREMAN_DB_DIR, f"HarremanDB_{species}_extracellular.csv"),
            known_hash=harreman_data_hash(
                HARREMAN_DB_DIR, f"HarremanDB_{species}_extracellular.csv"
            ),
            fname=f"HarremanDB_{species}_extracellular.csv",
            path=cache,
            progressbar=False,
        ),
        "all": pooch.retrieve(
            url=harreman_data_url(HARREMAN_DB_DIR, f"HarremanDB_{species}.csv"),
            known_hash=harreman_data_hash(HARREMAN_DB_DIR, f"HarremanDB_{species}.csv"),
            fname=f"HarremanDB_{species}.csv",
            path=cache,
            progressbar=False,
        ),
        "heterodimer": pooch.retrieve(
            url=harreman_data_url(HARREMAN_DB_DIR, f"Heterodimer_info_{species}.csv"),
            known_hash=harreman_data_hash(HARREMAN_DB_DIR, f"Heterodimer_info_{species}.csv"),
            fname=f"Heterodimer_info_{species}.csv",
            path=cache,
            progressbar=False,
        ),
    }
    df = pd.read_csv(filenames["extracellular" if extracellular_only else "all"], index_col=0)
    heterodimer_info = pd.read_csv(filenames["heterodimer"], index_col=0)
    df["Metabolite"] = df["Metabolite"].str.replace("/", "_", regex=False)

    pattern_import = compile(r"(\S+)" + import_suffix + "$")
    pattern_export = compile(r"(\S+)" + export_suffix + "$")

    metab_dict = {}
    for _, row in df.iterrows():
        metabolite, genes_str = row["Metabolite"], row["Gene"]
        if not genes_str:
            continue
        genes = [g.strip() for g in genes_str.split("/") if g.strip()]
        genes = sorted(set(genes))
        m = metabolite.lower()
        match_import = match(pattern_import, m)
        match_export = match(pattern_export, m)
        if match_import:
            name, direction = match_import.group(1), IMPORT_METAB_KEY
        elif match_export:
            name, direction = match_export.group(1), EXPORT_METAB_KEY
        else:
            name, direction = metabolite, BOTH_METAB_KEY
        metab_dict.setdefault(
            name, {IMPORT_METAB_KEY: [], EXPORT_METAB_KEY: [], BOTH_METAB_KEY: []}
        )
        genes_in_var = pd.Series(genes).isin(adata.var_names)
        metab_dict[name][direction] = pd.Series(genes)[genes_in_var].tolist()

    adata.uns["metabolite_database"] = df
    adata.uns["heterodimer_info"] = heterodimer_info
    adata.uns["num_metabolites"] = df.shape[0]
    adata.uns["importer"] = build_metabolite_df(metab_dict, IMPORT_METAB_KEY)
    adata.uns["exporter"] = build_metabolite_df(metab_dict, EXPORT_METAB_KEY)
    adata.uns["import_export"] = build_metabolite_df(metab_dict, BOTH_METAB_KEY)

    return metab_dict


def build_metabolite_df(metab_dict, key):
    """Build a metabolite annotation dataframe for one transport direction."""
    df = pd.DataFrame.from_dict({k: v[key] for k, v in metab_dict.items()}, orient="index").T
    df.columns = [f"{key}{i}" for i in range(df.shape[1])]
    return df


def extract_lr_pairs(adata, species):
    """Extracting LR pairs from CellChatDB."""
    cache = pooch.os_cache("scvi_harreman")

    interaction_fname = f"interaction_input_CellChatDB_v2_{species}.csv"
    complex_fname = f"complex_input_CellChatDB_v2_{species}.csv"

    interaction_path = pooch.retrieve(
        url=harreman_data_url(CELLCHAT_DB_DIR, interaction_fname),
        known_hash=harreman_data_hash(CELLCHAT_DB_DIR, interaction_fname),
        fname=interaction_fname,
        path=cache,
        progressbar=False,
    )
    complex_path = pooch.retrieve(
        url=harreman_data_url(CELLCHAT_DB_DIR, complex_fname),
        known_hash=harreman_data_hash(CELLCHAT_DB_DIR, complex_fname),
        fname=complex_fname,
        path=cache,
        progressbar=False,
    )

    interaction = pd.read_csv(interaction_path, index_col=0).sort_values("annotation")
    complex = pd.read_csv(complex_path, index_col=0)

    ligands, receptors = interaction.pop("ligand").values, interaction.pop("receptor").values

    for i in range(len(ligands)):
        for n in [ligands, receptors]:
            l = n[i]
            if l in complex.index:
                n[i] = (
                    complex.loc[l]
                    .dropna()
                    .values[pd.Series(complex.loc[l].dropna().values).isin(adata.var_names)]
                )
            else:
                n[i] = pd.Series(l).values[pd.Series(l).isin(adata.var_names)]

    lig_df = pd.DataFrame.from_records(zip_longest(*pd.Series(ligands).values)).transpose()
    rec_df = pd.DataFrame.from_records(zip_longest(*pd.Series(receptors).values)).transpose()

    lig_df.columns = [f"Ligand{i}" for i in range(lig_df.shape[1])]
    rec_df.columns = [f"Receptor{i}" for i in range(rec_df.shape[1])]

    lig_df.index = rec_df.index = interaction.index

    adata.uns["ligand"] = lig_df
    adata.uns["receptor"] = rec_df
    adata.uns["LR_database"] = interaction

    return
